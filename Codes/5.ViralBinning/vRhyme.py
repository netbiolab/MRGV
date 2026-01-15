#!/usr/bin/env python3

import os, sys
import subprocess
import pandas as pd
from multiprocessing import Process

### This code is for viral binning with vRhyme 
def run_vRhyme(svdir, sub_df):
    for idx in sub_df.index:
        print("start on %s" % idx)
        cohort = sub_df.loc[idx, "Cohort"]
        sample = sub_df.loc[idx, "Sample"]
        vctg = sub_df.loc[idx, "vContigs_path"]
        vrcov = sub_df.loc[idx, "VR_cov"]
        svpath1 = "/".join([svdir, cohort])
        
        svpath2 = "/".join([svpath1, sample])
        cmd = "vRhyme -t 1 -i %s -c %s -o %s" % (vctg, vrcov, svpath2)
        if not os.path.isdir(svpath1):
            os.makedirs(svpath1, exist_ok=True)
        subprocess.run([cmd], shell=True)


if __name__ == "__main__":
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    svdir = sys.argv[2]
    threads = int(sys.argv[3])
    procs = []
    for t in range(threads):
        sub_df = meta.iloc[t::threads, :]
        proc = Process(target=run_vRhyme, args=(svdir, sub_df))
        proc.start()
        procs.append(proc)
    for proc in procs:
        proc.join()
