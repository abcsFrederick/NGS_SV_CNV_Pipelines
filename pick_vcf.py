# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 12:15:49 2018

@author: Jack
"""
import sys

filepath = sys.argv[1]
filepath1 = sys.argv[2]
lw = sys.argv[3]
out = ""

fl=open(filepath, 'r')
ll = fl.readlines()
fl.close()

for ln in ll:
    if ln.startswith("#"):
        out = out + ln
    else:
        if lw == "lumpy":
            if ln.strip().split("\t")[-1] == "./.:0:0:0" and ln.strip().split("\t")[-2] != "./.:0:0:0":
                out = out + ln
        else:
            if ln.strip().split("\t")[-1] == ".:.:0" and ln.strip().split("\t")[-2] != ".:.:0":
                out = out + ln
        
fout = open(filepath1, 'w')
fout.write(out)
fout.close()