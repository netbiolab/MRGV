#!/usr/bin/env python3

import os, sys, pickle
import pandas as pd
from multiprocessing import Process
import subprocess
from subprocess import PIPE

## This code for viral binning with Metabat2
def run_MB2(svdir, sub_list):
    for idx in sub_list.index:
        print("Starts on %s" % idx)
        cohort = sub_list.loc[idx, "Cohort"]
        sample = sub_list.loc[idx, "Sample"]
        vctg = sub_list.loc[idx, "vContigs_path"]
        depth = sub_list.loc[idx, "MB2_depth"]
        svpath1 = "/".join([svdir, cohort, sample])
        fname = idx
        svpath2 = "/".join([svpath1, fname])
        if not os.path.isdir(svpath1):
            os.makedirs(svpath1, exist_ok=True)
        mb2_cmd = "metabat2 -i %s -m 2000 -s 5000 -t 1 -a %s -o %s" % (vctg, depth, svpath2)
        out = subprocess.run([mb2_cmd], shell=True, stdout=PIPE, stderr=PIPE)
        #os.system(mb2_cmd)
    

if __name__ == "__main__":
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    svdir = sys.argv[2]
    threads = int(sys.argv[3])
    procs = []
    for t in range(threads):
        sub_list = meta.iloc[t::threads, :]
        proc = Process(target=run_MB2, args=(svdir, sub_list))
        proc.start()
        procs.append(proc)
    for proc in procs:
        proc.join()


