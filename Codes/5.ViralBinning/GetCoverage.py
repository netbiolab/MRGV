#!/usr/bin/env python3

import os, sys, pickle
import pandas as pd
import subprocess

### This code is for computing sample-wise contig coverage 

## Buile sample-wise bowtie2 index
def get_btidx(vctg, sv_dir, idx, threads):
    btidx_dir = "/".join([sv_dir, "btidx"])
    
    if not os.path.isdir(btidx_dir):
        os.makedirs(btidx_dir)
    btidx_name = "/".join([btidx_dir, idx])
     
    cmd = "bowtie2-build -q --threads %s %s %s" % (threads, vctg, btidx_name)
    print("Start to create bowtie2 index on %s " % idx)
    subprocess.run([cmd], shell=True)
    #print(cmd)
    print("Bowtie2 index of %s has been created\n" % idx)
    return btidx_name

## Run Bowtie2 for each sample
def run_bowtie2(sv_dir, idx, threads, btidx, f1, f2, stype):
    bam_path = "/".join([sv_dir, "%s.bam" % idx])
    sbam_path = "/".join([sv_dir, "%s_sorted.bam" % idx])
    
    if os.path.isfile(sbam_path):
        print("file %s exists, so skipped" % sbam_path)
    else:
        if stype == "paired":
            cmd = "bowtie2 --quiet --threads %s -x %s -1 %s -2 %s|samtools view -bS - > %s" % (threads, btidx, f1, f2, bam_path)
        else:
            cmd = "bowtie2 --quiet --threads %s -x %s -U %s|samtools view -bS - > %s" % (threads,  btidx, f1, bam_path)
        print("Start to align with bowtie2 on %s" % idx)
        subprocess.run([cmd], shell=True)
        #print(cmd)
        print("Align on %s has been done\n" % idx)
        print("Start to sort and index")
        cmd = "samtools sort %s -o %s" % (bam_path, sbam_path)
        subprocess.run([cmd], shell=True)
        #print(cmd)
        print("Bam file of %s has been sorted\n" % idx)
        print("Start to index of sorted bam file of %s" % idx)
        cmd = "samtools index %s" % (sbam_path)
        subprocess.run([cmd], shell=True)
        #print(cmd)
        print("Bam alignment profile on %s has been created!\n" % idx)
        print("Start to remove original bamfiles")
        cmd = "rm %s" % bam_path
        subprocess.run([cmd], shell=True)
        print("Successfully removed\n")


def main(sub_df, sv_dir, threads):
    for idx in sub_df.index:
        cohort = sub_df.loc[idx, "Cohort"]
        sample = sub_df.loc[idx, "Sample"]
        f1 = sub_df.loc[idx, "fq1"]
        f2 = sub_df.loc[idx, "fq2"]
        stype = sub_df.loc[idx, "stype"]
        vctg = sub_df.loc[idx, "FF_fa_path"]
        if vctg == "None":
            continue
        sv_dir1 = "/".join([sv_dir, cohort, sample])
        if not os.path.isdir(sv_dir1):
            os.makedirs(sv_dir1)
        sbam_path = "/".join([sv_dir, "%s_sorted.bam" % idx])
        if os.path.isfile(sbam_path):
            print("file %s exist, so skipped" % sbam_path)
        else:            
            btidx_name = get_btidx(vctg, sv_dir1, idx, threads)
            sbam_path = run_bowtie2(sv_dir1, idx, threads, btidx_name, f1, f2, stype)

if __name__ == "__main__":
    meta_df = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    sv_dir = sys.argv[2]
    from_ = int(sys.argv[3]) 
    interval_ = int(sys.argv[4])
    threads = int(sys.argv[5])
    sub_df = meta_df.iloc[from_::interval_, :]
    main(sub_df, sv_dir, threads)


