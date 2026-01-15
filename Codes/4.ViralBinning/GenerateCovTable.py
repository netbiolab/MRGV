#! /usr/bin/env python3
# vRhyme
# Author: Kristopher Kieft
# University of Wisconsin-Madison


import argparse
import os
import sys
import pandas as pd
import subprocess
from multiprocessing import Process

### This code is for converting Metabat2 depth file into vRhyme coverage file for running vRhyme
def convert(svdir, sub_df):
    for idx in sub_df.index:
        cohort = sub_df.loc[idx, "Cohort"]
        sample = sub_df.loc[idx, "Sample"]
        depth= sub_df.loc[idx, "MB2_depth"]
        svpath = "/".join([svdir, cohort, sample])
        if not os.path.isdir(svpath):
            os.makedirs(svpath, exist_ok=True)
        fname = idx+"_VR_cov.tsv"
        svname = "/".join([svpath, fname])

        with open(depth, 'r') as in_table, open(svname, 'w') as out_table:
            first = in_table.readline().strip("\n").split("\t")[3:]
            total = len(first)+3
            out_table.write("scaffold")
            for item in first:
                if item[-4:] == '-var':
                    out_table.write("\tstdev_" + str(item.rsplit("-",1)[0]))
                else:
                    out_table.write("\tavg_" + str(item))
            for line in in_table:
                line = line.strip("\n").split("\t")
                out_table.write("\n" + line[0])
                for n in range(3,total,2):
                    out_table.write("\t" + line[n])
                    sd = float(line[n+1])**0.5
                    out_table.write("\t" + str(sd))


if __name__ == '__main__':
    meta = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    svdir = sys.argv[2]
    threads = int(sys.argv[3])
    procs = []
    for t in range(threads):
        sub_df = meta.iloc[t::threads, :]
        proc = Process(target=convert, args=(svdir, sub_df))
        proc.start()
        procs.append(proc)
    for proc in procs:
        proc.join()
