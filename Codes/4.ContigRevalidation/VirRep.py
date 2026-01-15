#!/usr/bin/env python3

import os, sys, pickle, subprocess, time
import pandas as pd

### This code is for running VirRep in GPU server on deduplicated viral contigs

df = sys.argv[1]
name = df.split("/")[-1].split(".fa")[0]
svdir = sys.argv[2]
threads = int(sys.argv[3])
gpu = sys.argv[4]
model = sys.argv[5]  ## Virrep model path
svpath ="/".join([svdir, name])
if not os.path.isdir(svdir):
    os.makedirs(svpath, exist_ok=True)

cmd = "VirRep.py -m %s -i %s -o %s --provirus-off --use-amp -w %s -k 2000 --batch-size 4096 --gpu-device %s -l 2000 --label %s" % (model, df, svpath, threads, gpu, name)
st = time.time()
subprocess.run([cmd], shell=True)
et = time.time() -st
print("done in %s min" % (et/60))
