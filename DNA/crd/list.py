import os, sys

f = open("list")
lines = f.readlines()
f.close()

for line in lines:
    term = line[:-1]
    os.chdir(term)
    os.system(f"rm *arc")
    os.system(f"rm *dyn")
    os.system(f"rm *sh")
    os.system(f"rm *log")
    os.system(f"rm r_MD-{term}.key")
    os.system(f"rm read_config.py")
    os.system(f"rm {term}-results.log")
    os.system(f"rm *dat")
    os.system("rm *traj*")
    os.chdir("..")
