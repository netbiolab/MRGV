#!/usr/bin/env python3

import os, sys
import pandas as pd
from multiprocessing import Process
import subprocess

### This code is for generating Metabat2 depth file for viral binning prerequisite

def map_cmd(sub_meta, svdir):
    task= []
    for i in sub_meta.index:
        cohort = sub_meta.loc[i, "Cohort"]
        sample = sub_meta.loc[i, "Sample"]
        vctg = sub_meta.loc[i, "FF_fa_path"]
        sbam = sub_meta.loc[i, "sbam_path"]
        svpath1 = "/".join([svdir, cohort, sample])
        print("Starts on %s" % i)
        if not os.path.isdir(svpath1):
            os.makedirs(svpath1, exist_ok=True)
        fname = ".".join([i, "depth"])
        out_depth = "/".join([svpath1, fname])
        cmd = "jgi_summarize_bam_contig_depths --minContigLength 2000 --referenceFasta {} {} --outputDepth {}".format(vctg, sbam, out_depth)
        out = subprocess.run([cmd], shell=True, capture_output=True)
        print("Done on %s" % i)

if __name__ == "__main__":
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    svdir = sys.argv[2]
    threads = int(sys.argv[3])
    procs = []    
    for t in range(threads):
        sub_meta = meta.iloc[t::threads, :]
        proc = Process(target=map_cmd, args=(sub_meta, svdir))
        proc.start()
        procs.append(proc)
    for proc in procs:
        proc.join()
        




