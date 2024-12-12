import os
import sys
import json
import shutil
import datetime
import subprocess
from rich import print
from pathlib import Path
    

def create_dir(path):
    if not os.path.exists(path): 
        os.makedirs(path) 
        
def run(cmd : str):
    subprocess.run(
        cmd, 
        shell=True)

if __name__ == "__main__":
    config_filepath = sys.argv[1]
    f = open(config_filepath)
    config = json.load(open(config_filepath))
    
    print(json.dumps(config, indent=2))
    
    config_name = Path(config_filepath).stem
    gmx_path = config['gmx_dir'] + "/gmx"
    setting_dir = config['setting_dir']
    nb = config['nb']
    
    dt_now = datetime.datetime.now().strftime('%Y%m%d%H%M%S')
    result_dir = f"exp/result/{dt_now}@{config_name}"
    
    os.makedirs(result_dir, exist_ok=True)
    os.makedirs(f"{result_dir}/{config_name}", exist_ok=True)
    shutil.copyfile(config_filepath, f"{result_dir}/{config_name}/config.json")
    shutil.copytree(setting_dir, f"{result_dir}/{config_name}/setting")
    
    print("[orange1]start[/orange1]")
    
    cmd1 = f"{gmx_path} grompp -f {setting_dir}/grompp.mdp -p {setting_dir}/topol.top -c {setting_dir}/conf.gro -o {result_dir}/run.tpr -po {result_dir}/mdout.mdp"
    print(f"{cmd1=}")
    run(cmd1)
    
    cmd2 = f"{gmx_path} mdrun -v -s {result_dir}/run.tpr -e {result_dir}/result.edr -g {result_dir}/log.log -c {result_dir}/confout.gro -cpo {result_dir}/state.cpt -nb {nb} "
    print(f"{cmd2=}")
    run(cmd2)
    
    print("[orange1]completed[/orange1]")
    
