#! /usr/bin/env python3

import pandas as pd
import os, sys, pickle, subprocess
import time
import numpy as np

def align(sub_meta, mm_idx, svdir, threads):
    print("Start sub %s samples" % len(sub_meta))
    times = []
    for n, i in enumerate(sub_meta.index):
        cohort = sub_meta.loc[i, "Cohort"]
        sample = sub_meta.loc[i, "Sample"]
        name_base = "-".join([cohort, sample])
        svdir2 = "/".join([svdir, cohort, sample])
        if not os.path.isdir(svdir2):
            os.makedirs(svdir2, exist_ok=True)
        svpath = "/".join([svdir2, name_base+".sam"])
        stype = sub_meta.loc[i, "stype"]
        if stype == "paired":
            fq1 = sub_meta.loc[i, "fq1"]
            fq2 = sub_meta.loc[i, "fq2"]
            if not os.path.isfile(fq1):
                fq1 = fq1 + ".gz"
            if not os.path.isfile(fq2):
                fq2 = fq2 + ".gz"

            cmd = "minimap2 -a sr -t %s %s %s %s >> %s" % (threads, mm_idx, fq1, fq2, svpath)
        elif stype == "single":
            fq = sub_meta.loc[i, "fq1"]
            if not os.path.isfile(fq):
                fq = fq + ".gz"
            cmd = "minimap2 -a sr -t %s %s %s >> %s" % (threads, mm_idx, fq, svpath)
        
        cmd2 = "samtools view -bS -@ %s %s | samtools sort -@ %s --write-index -o %s" % (threads, svpath, threads, svpath.replace(".sam", ".bam.sorted"))
        cmd3 = "rm %s" % svpath
        
        cmd_all = " && ".join([cmd, cmd2, cmd3])
        st = time.time()
        if not os.path.isfile(svpath.replace(".sam", ".bam.sorted")):
            if os.path.isfile(svpath):
                os.remove(svpath)
        #subprocess.run([cmd_all], shell=True)
        print(cmd_all)
        break
        et = (time.time() - st)/60
        times.append(et)
        if (n+1) % 10 == 0:
            print("Finished %s samples" % (n+1))
            print("Average time: %s" % np.mean(times))
    print("Finished all %s samples" % len(sub_meta))
    
    return None
        



if __name__ == "__main__":
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    mm_idx = sys.argv[2]
    svdir = sys.argv[3]
    from_ = int(sys.argv[4])
    interval_ = int(sys.argv[5])
    threads = int(sys.argv[6])
    if not os.path.isdir(svdir):
        os.makedirs(svdir, exist_ok=True)
    align(meta.iloc[from_::interval_], mm_idx, svdir, threads)
    
