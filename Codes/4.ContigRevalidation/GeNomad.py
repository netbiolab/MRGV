#!/usr/bin/env python3

### This code is for running genomad on deduplicated viral contigs
import time
import os, sys, subprocess
svdir = sys.argv[1]
ctg_list = sys.argv[2]
threads = sys.argv[3]
genomad_db = sys.argv[4] 
outdir = "/".join([svdir, "%s" % ctg_list.split("/")[-1].split(".fa")[0]])
if not os.path.isdir(outdir):
    os.makedirs(outdir, exist_ok=True)
cmd = "genomad end-to-end --enable-score-calibration --threads %s %s %s %s" % (threads, ctg_list, outdir, genomad_db)
st = time.time()
subprocess.run([cmd], shell=True)
#print(cmd)
et = time.time() - st
print("Done in %s min" % (et/60))
