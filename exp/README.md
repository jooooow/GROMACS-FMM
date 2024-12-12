### Example

```sh
python3 exp/rungromacs.py exp/noLJ.json 
```

This python script will
+ create result directory named `exp/result/<date>@<config>/` 
+ copy `<config>.json` and settings into the result directory
+ save all results in the result directory