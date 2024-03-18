import os, sys

f = open("list")
lines = f.readlines()
f.close()

for i, line in enumerate(lines):
    term = line[:-1]
    os.chdir(term)
    os.system("cp %s.xyz ../../anchors/%d.xyz" %(term, i+1))
    os.chdir("..")
