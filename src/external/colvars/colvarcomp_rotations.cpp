// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"
#include "colvar_rotation_derivative.h"


struct colvar::orientation::rotation_derivative_impl_: public rotation_derivative<cvm::atom_pos, cvm::atom_pos> {
public:
  rotation_derivative_impl_(colvar::orientation* orientation_cvc):
   rotation_derivative<cvm::atom_pos, cvm::atom_pos>(
    orientation_cvc->rot, orientation_cvc->ref_pos, orientation_cvc->shifted_pos) {}
};


colvar::orientation::orientation(std::string const &conf)
  : cvc()
{
  set_function_type("orientation");
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_quaternion);
  colvar::orientation::init(conf);
}


colvar::orientation::~orientation() {}


int colvar::orientation::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  atoms = parse_group(conf, "atoms");
  ref_pos.reserve(atoms->size());

  if (get_keyval(conf, "refPositions", ref_pos, ref_pos)) {
    cvm::log("Using reference positions from input file.\n");
    if (ref_pos.size() != atoms->size()) {
      return cvm::error("Error: reference positions do not "
                        "match the number of requested atoms.\n", COLVARS_INPUT_ERROR);
    }
  }

  {
    std::string file_name;
    if (get_keyval(conf, "refPositionsFile", file_name)) {

      std::string file_col;
      double file_col_value=0.0;
      if (get_keyval(conf, "refPositionsCol", file_col, std::string(""))) {
        // use PDB flags if column is provided
        bool found = get_keyval(conf, "refPositionsColValue", file_col_value, 0.0);
        if (found && file_col_value==0.0) {
          return cvm::error("Error: refPositionsColValue, "
                            "if provided, must be non-zero.\n", COLVARS_INPUT_ERROR);
        }
      }

      ref_pos.resize(atoms->size());
      cvm::load_coords(file_name.c_str(), &ref_pos, atoms,
                       file_col, file_col_value);
    }
  }

  if (!ref_pos.size()) {
    return cvm::error("Error: must define a set of "
                      "reference coordinates.\n", COLVARS_INPUT_ERROR);
  }


  cvm::rvector ref_cog(0.0, 0.0, 0.0);
  size_t i;
  for (i = 0; i < ref_pos.size(); i++) {
    ref_cog += ref_pos[i];
  }
  ref_cog /= cvm::real(ref_pos.size());
  cvm::log("Centering the reference coordinates on the origin by subtracting "
           "the center of geometry at "+
           cvm::to_str(-1.0 * ref_cog)+"; it is "
           "assumed that each atom is the closest "
           "periodic image to the center of geometry.\n");
  for (i = 0; i < ref_pos.size(); i++) {
    ref_pos[i] -= ref_cog;
  }

  get_keyval(conf, "closestToQuaternion", ref_quat, cvm::quaternion(1.0, 0.0, 0.0, 0.0));

  rot_deriv_impl = std::unique_ptr<rotation_derivative_impl_>(new rotation_derivative_impl_(this));

  // If the debug gradients feature is active, debug the rotation gradients
  // (note that this won't be active for the orientation CVC itself, because
  // colvardeps prevents the flag's activation)
  rot.b_debug_gradients = is_enabled(f_cvc_debug_gradient);

  return error_code;
}


colvar::orientation::orientation()
  : cvc()
{
  set_function_type("orientation");
  disable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_quaternion);
}


void colvar::orientation::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  if ((rot.q).inner(ref_quat) >= 0.0) {
    x.quaternion_value = rot.q;
  } else {
    x.quaternion_value = -1.0 * rot.q;
  }
}


void colvar::orientation::calc_gradients()
{
  // gradients have already been calculated and stored within the
  // member object "rot"; we're not using the "grad" member of each
  // atom object, because it only can represent the gradient of a
  // scalar colvar
}


void colvar::orientation::apply_force(colvarvalue const &force)
{
  cvm::quaternion const &FQ = force.quaternion_value;

  if (!atoms->noforce) {
    rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
    cvm::vector1d<cvm::rvector> dq0_2;
    auto ag_force = atoms->get_group_force_object();
    for (size_t ia = 0; ia < atoms->size(); ia++) {
      rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
      const auto f_ia = FQ[0] * dq0_2[0] +
                        FQ[1] * dq0_2[1] +
                        FQ[2] * dq0_2[2] +
                        FQ[3] * dq0_2[3];
      ag_force.add_atom_force(ia, f_ia);
    }
  }
}


cvm::real colvar::orientation::dist2(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2(x2);
}


colvarvalue colvar::orientation::dist2_lgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x1.quaternion_value.dist2_grad(x2);
}


colvarvalue colvar::orientation::dist2_rgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  return x2.quaternion_value.dist2_grad(x1);
}



colvar::orientation_angle::orientation_angle(std::string const &conf)
  : orientation()
{
  set_function_type("orientationAngle");
  init_as_angle();
  enable(f_cvc_explicit_gradient);
  orientation_angle::init(conf);
}


int colvar::orientation_angle::init(std::string const &conf)
{
  return orientation::init(conf);
}


void colvar::orientation_angle::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  if ((rot.q).q0 >= 0.0) {
    x.real_value = (180.0/PI) * 2.0 * cvm::acos((rot.q).q0);
  } else {
    x.real_value = (180.0/PI) * 2.0 * cvm::acos(-1.0 * (rot.q).q0);
  }
}


void colvar::orientation_angle::calc_gradients()
{
  cvm::real const dxdq0 =
    ( ((rot.q).q0 * (rot.q).q0 < 1.0) ?
      ((180.0 / PI) * (-2.0) / cvm::sqrt(1.0 - ((rot.q).q0 * (rot.q).q0))) :
      0.0 );

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    (*atoms)[ia].grad = (dxdq0 * dq0_2[0]);
  }
}


void colvar::orientation_angle::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;
  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


simple_scalar_dist_functions(orientation_angle)



colvar::orientation_proj::orientation_proj(std::string const &conf)
  : orientation()
{
  set_function_type("orientationProj");
  enable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
  orientation_proj::init(conf);
}


int colvar::orientation_proj::init(std::string const &conf)
{
  return orientation::init(conf);
}


void colvar::orientation_proj::calc_value()
{
  atoms_cog = atoms->center_of_geometry();
  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);
  x.real_value = 2.0 * (rot.q).q0 * (rot.q).q0 - 1.0;
}


void colvar::orientation_proj::calc_gradients()
{
  cvm::real const dxdq0 = 2.0 * 2.0 * (rot.q).q0;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    (*atoms)[ia].grad = (dxdq0 * dq0_2[0]);
  }
}


void colvar::orientation_proj::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


simple_scalar_dist_functions(orientation_proj)



colvar::tilt::tilt(std::string const &conf)
  : orientation()
{
  set_function_type("tilt");
  x.type(colvarvalue::type_scalar);
  enable(f_cvc_explicit_gradient);
  init_scalar_boundaries(-1.0, 1.0);
  tilt::init(conf);
}


int colvar::tilt::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  error_code |= orientation::init(conf);

  get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0));
  if (axis.norm2() != 1.0) {
    axis /= axis.norm();
    cvm::log("Normalizing rotation axis to "+cvm::to_str(axis)+".\n");
  }

  return error_code;
}


void colvar::tilt::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  x.real_value = rot.cos_theta(axis);
}


void colvar::tilt::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dcos_theta_dq(axis);

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    (*atoms)[ia].grad = cvm::rvector(0.0, 0.0, 0.0);
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    for (size_t iq = 0; iq < 4; iq++) {
      (*atoms)[ia].grad += (dxdq[iq] * dq0_2[iq]);
    }
  }
}


void colvar::tilt::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


simple_scalar_dist_functions(tilt)



colvar::spin_angle::spin_angle(std::string const &conf)
  : orientation()
{
  set_function_type("spinAngle");
  init_as_periodic_angle();
  enable(f_cvc_periodic);
  enable(f_cvc_explicit_gradient);
  spin_angle::init(conf);
}


int colvar::spin_angle::init(std::string const &conf)
{
  int error_code = COLVARS_OK;

  error_code |= orientation::init(conf);

  get_keyval(conf, "axis", axis, cvm::rvector(0.0, 0.0, 1.0));
  if (axis.norm2() != 1.0) {
    axis /= axis.norm();
    cvm::log("Normalizing rotation axis to "+cvm::to_str(axis)+".\n");
  }

  return error_code;
}


colvar::spin_angle::spin_angle()
  : orientation()
{
  set_function_type("spinAngle");
  period = 360.0;
  enable(f_cvc_periodic);
  enable(f_cvc_explicit_gradient);
  x.type(colvarvalue::type_scalar);
}


void colvar::spin_angle::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  x.real_value = rot.spin_angle(axis);
  this->wrap(x);
}


void colvar::spin_angle::calc_gradients()
{
  cvm::quaternion const dxdq = rot.dspin_angle_dq(axis);

  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    (*atoms)[ia].grad = cvm::rvector(0.0, 0.0, 0.0);
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    for (size_t iq = 0; iq < 4; iq++) {
      (*atoms)[ia].grad += (dxdq[iq] * dq0_2[iq]);
    }
  }
}


void colvar::spin_angle::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;

  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


cvm::real colvar::spin_angle::dist2(colvarvalue const &x1,
                                    colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return diff * diff;
}


colvarvalue colvar::spin_angle::dist2_lgrad(colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return 2.0 * diff;
}


colvarvalue colvar::spin_angle::dist2_rgrad(colvarvalue const &x1,
                                            colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return (-2.0) * diff;
}


void colvar::spin_angle::wrap(colvarvalue &x_unwrapped) const
{
  if ((x_unwrapped.real_value - wrap_center) >= 180.0) {
    x_unwrapped.real_value -= 360.0;
    return;
  }

  if ((x_unwrapped.real_value - wrap_center) < -180.0) {
    x_unwrapped.real_value += 360.0;
    return;
  }

  return;
}


colvar::euler_phi::euler_phi(std::string const &conf)
  : orientation()
{
  set_function_type("eulerPhi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
  euler_phi::init(conf);
}


colvar::euler_phi::euler_phi()
  : orientation()
{
  set_function_type("eulerPhi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
}


int colvar::euler_phi::init(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= orientation::init(conf);
  return error_code;
}


void colvar::euler_phi::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real tmp_y = 2 * (q0 * q1 + q2 * q3);
  const cvm::real tmp_x = 1 - 2 * (q1 * q1 + q2 * q2);
  x.real_value = cvm::atan2(tmp_y, tmp_x) * (180.0/PI);
}


void colvar::euler_phi::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = (2 * q0 * q1 + 2 * q2 * q3) * (2 * q0 * q1 + 2 * q2 * q3) + (-2 * q1 * q1 - 2 * q2 * q2 + 1) * (-2 * q1 * q1 - 2 * q2 * q2 + 1);
  const cvm::real dxdq0 = (180.0/PI) * 2 * q1 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) / denominator;
  const cvm::real dxdq1 = (180.0/PI) * (2 * q0 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) - 4 * q1 * (-2 * q0 * q1 - 2 * q2 * q3)) / denominator;
  const cvm::real dxdq2 = (180.0/PI) * (-4 * q2 * (-2 * q0 * q1 - 2 * q2 * q3) + 2 * q3 * (-2 * q1 * q1 - 2 * q2 * q2 + 1)) / denominator;
  const cvm::real dxdq3 = (180.0/PI) * 2 * q2 * (-2 * q1 * q1 - 2 * q2 * q2 + 1) / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    (*atoms)[ia].grad = (dxdq0 * dq0_2[0]) +
                        (dxdq1 * dq0_2[1]) +
                        (dxdq2 * dq0_2[2]) +
                        (dxdq3 * dq0_2[3]);
  }
}


void colvar::euler_phi::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;
  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


cvm::real colvar::euler_phi::dist2(colvarvalue const &x1,
                                  colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return diff * diff;
}


colvarvalue colvar::euler_phi::dist2_lgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return 2.0 * diff;
}


colvarvalue colvar::euler_phi::dist2_rgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return (-2.0) * diff;
}


void colvar::euler_phi::wrap(colvarvalue &x_unwrapped) const
{
  if ((x_unwrapped.real_value - wrap_center) >= 180.0) {
    x_unwrapped.real_value -= 360.0;
    return;
  }

  if ((x_unwrapped.real_value - wrap_center) < -180.0) {
    x_unwrapped.real_value += 360.0;
    return;
  }

  return;
}


colvar::euler_psi::euler_psi(std::string const &conf)
  : orientation()
{
  set_function_type("eulerPsi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
  euler_psi::init(conf);
}


colvar::euler_psi::euler_psi()
  : orientation()
{
  set_function_type("eulerPsi");
  init_as_periodic_angle();
  enable(f_cvc_explicit_gradient);
}


int colvar::euler_psi::init(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= orientation::init(conf);
  return error_code;
}


void colvar::euler_psi::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real tmp_y = 2 * (q0 * q3 + q1 * q2);
  const cvm::real tmp_x = 1 - 2 * (q2 * q2 + q3 * q3);
  x.real_value = cvm::atan2(tmp_y, tmp_x) * (180.0/PI);
}


void colvar::euler_psi::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = (2 * q0 * q3 + 2 * q1 * q2) * (2 * q0 * q3 + 2 * q1 * q2) + (-2 * q2 * q2 - 2 * q3 * q3 + 1) * (-2 * q2 * q2 - 2 * q3 * q3 + 1);
  const cvm::real dxdq0 = (180.0/PI) * 2 * q3 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) / denominator;
  const cvm::real dxdq1 = (180.0/PI) * 2 * q2 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) / denominator;
  const cvm::real dxdq2 = (180.0/PI) * (2 * q1 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) - 4 * q2 * (-2 * q0 * q3 - 2 * q1 * q2)) / denominator;
  const cvm::real dxdq3 = (180.0/PI) * (2 * q0 * (-2 * q2 * q2 - 2 * q3 * q3 + 1) - 4 * q3 * (-2 * q0 * q3 - 2 * q1 * q2)) / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    (*atoms)[ia].grad = (dxdq0 * dq0_2[0]) +
                        (dxdq1 * dq0_2[1]) +
                        (dxdq2 * dq0_2[2]) +
                        (dxdq3 * dq0_2[3]);
  }
}


void colvar::euler_psi::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;
  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


cvm::real colvar::euler_psi::dist2(colvarvalue const &x1,
                                  colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return diff * diff;
}


colvarvalue colvar::euler_psi::dist2_lgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return 2.0 * diff;
}


colvarvalue colvar::euler_psi::dist2_rgrad(colvarvalue const &x1,
                                          colvarvalue const &x2) const
{
  cvm::real diff = x1.real_value - x2.real_value;
  diff = (diff < -180.0 ? diff + 360.0 : (diff > 180.0 ? diff - 360.0 : diff));
  return (-2.0) * diff;
}


void colvar::euler_psi::wrap(colvarvalue &x_unwrapped) const
{
  if ((x_unwrapped.real_value - wrap_center) >= 180.0) {
    x_unwrapped.real_value -= 360.0;
    return;
  }

  if ((x_unwrapped.real_value - wrap_center) < -180.0) {
    x_unwrapped.real_value += 360.0;
    return;
  }

  return;
}


colvar::euler_theta::euler_theta(std::string const &conf)
  : orientation()
{
  set_function_type("eulerTheta");
  init_as_angle();
  enable(f_cvc_explicit_gradient);
  euler_theta::init(conf);
}


colvar::euler_theta::euler_theta()
  : orientation()
{
  set_function_type("eulerTheta");
  init_as_angle();
  enable(f_cvc_explicit_gradient);
}


int colvar::euler_theta::init(std::string const &conf)
{
  int error_code = COLVARS_OK;
  error_code |= orientation::init(conf);
  return error_code;
}


void colvar::euler_theta::calc_value()
{
  atoms_cog = atoms->center_of_geometry();

  shifted_pos = atoms->positions_shifted(-1.0 * atoms_cog);
  rot.calc_optimal_rotation(ref_pos, shifted_pos);

  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  x.real_value = cvm::asin(2 * (q0 * q2 - q3 * q1)) * (180.0/PI);
}


void colvar::euler_theta::calc_gradients()
{
  const cvm::real& q0 = rot.q.q0;
  const cvm::real& q1 = rot.q.q1;
  const cvm::real& q2 = rot.q.q2;
  const cvm::real& q3 = rot.q.q3;
  const cvm::real denominator = cvm::sqrt(1 - (2 * q0 * q2 - 2 * q1 * q3) * (2 * q0 * q2 - 2 * q1 * q3));
  const cvm::real dxdq0 = (180.0/PI) * 2 * q2 / denominator;
  const cvm::real dxdq1 = (180.0/PI) * -2 * q3 / denominator;
  const cvm::real dxdq2 = (180.0/PI) * 2 * q0 / denominator;
  const cvm::real dxdq3 = (180.0/PI) * -2 * q1 / denominator;
  rot_deriv_impl->prepare_derivative(rotation_derivative_dldq::use_dq);
  cvm::vector1d<cvm::rvector> dq0_2;
  for (size_t ia = 0; ia < atoms->size(); ia++) {
    rot_deriv_impl->calc_derivative_wrt_group2(ia, nullptr, &dq0_2);
    (*atoms)[ia].grad = (dxdq0 * dq0_2[0]) +
                        (dxdq1 * dq0_2[1]) +
                        (dxdq2 * dq0_2[2]) +
                        (dxdq3 * dq0_2[3]);
  }
}


void colvar::euler_theta::apply_force(colvarvalue const &force)
{
  cvm::real const &fw = force.real_value;
  if (!atoms->noforce) {
    atoms->apply_colvar_force(fw);
  }
}


cvm::real colvar::euler_theta::dist2(colvarvalue const &x1,
                                     colvarvalue const &x2) const
{
  // theta angle is not periodic
  return cvc::dist2(x1, x2);
}


colvarvalue colvar::euler_theta::dist2_lgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  // theta angle is not periodic
  return cvc::dist2_lgrad(x1, x2);
}


colvarvalue colvar::euler_theta::dist2_rgrad(colvarvalue const &x1,
                                             colvarvalue const &x2) const
{
  // theta angle is not periodic
  return cvc::dist2_rgrad(x1, x2);
}
