### Example

```sh
python3 exp/rungromacs.py exp/noLJ.json 
```

This python script will
+ creates result directory named `exp/result/<date>@<config>/` 
+ copys `<config>.json` and settings into the result directory
+ saves all results in the result directory