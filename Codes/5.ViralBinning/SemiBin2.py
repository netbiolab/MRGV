#!/usr/bin/env python3

import os, sys
import pandas as pd
from multiprocessing import Process
import subprocess

### This code is for viral binning with Semibin2
### To apply to viral contigs, checkM step of Sembin2 was removed
def run_SM(svdir, sub_df):
    for idx in sub_df.index:
        print("Start on %s" % idx)
        cohort = sub_df.loc[idx, "Cohort"]
        sample = sub_df.loc[idx, "Sample"]
        vctg = sub_df.loc[idx, "vContigs_path"]
        sbam = sub_df.loc[idx, "SBam"]
        depth = sub_df.loc[idx, "MB2_depth"]
        svpath1 = "/".join([svdir, cohort, sample])
        if not os.path.isdir(svpath1):
            os.makedirs(svpath1, exist_ok=True)
        ## Bypass checkm step in Semibin2 due to viral contigs
        cmd = "SemiBin2 single_easy_bin --quiet -m 2000 --compression none --minfasta-kbs 5 --engine cpu --self-supervised -i %s  -b %s -o %s --threads %s" % (vctg, sbam, svpath1, 1)
        out = subprocess.run([cmd], shell=True, capture_output=True)
        #print(cmd)

if __name__ == "__main__":
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    svdir = sys.argv[2]
    threads = int(sys.argv[3])
    from_ = int(sys.argv[4])
    interval_ = int(sys.argv[5])
    procs = []
    for t in range(threads):
        sub_df = meta.iloc[t::threads, :]
        proc = Process(target=run_SM, args=(svdir, sub_df))
        proc.start()
        procs.append(proc)
    for proc in procs:
        proc.join()
