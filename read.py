#!/usr/bin/env python
from lecroy import *
from sys import argv


DataDirectory = "testdata/"
OutputRootFile = "outdata.root"
for i in range(len(argv)):
    if argv[i] in ['-i','--indir']:
        DataDirectory=argv[i+1]
    if argv[i] in ['-o','--outfile']:
        OutputRootFile=argv[i+1]

if DataDirectory[-1] != "/":
    DataDirectory = DataDirectory + "/"
print "reading from : " + DataDirectory
print "output file: " + OutputRootFile

LeCroy2Root(DataDirectory, OutputRootFile)
