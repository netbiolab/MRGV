#!/usr/bin/env python3
import os, sys, pickle
import pandas as pd

### This code for running DeepVirFinder and filtering confident viral contigs

def run_dvf(ctg_name, cohort, ctg_path, threads, sv_dir):
    dvf_dir = "/".join([sv_dir, cohort, ctg_name, "dvf"])
    fname = "/".join([dvf_dir, ctg_path.split("/")[-1]+"_gt2000bp_dvfpred.txt"])
    filtered = "/".join([dvf_dir, ctg_name+"_filtered.fa"])
    if os.path.isfile(filtered):
        print("file %s exists, skipped" % filtered)
        return
    if not os.path.isdir(dvf_dir):
        os.makedirs(dvf_dir, exist_ok=True)
      
    ## DeepVirFinder parameter settings
    cmd = "DeepVirFinder/dvf.py -i %s -o %s -l 2000 -c %s" % (ctg_path, dvf_dir, threads)
    os.system(cmd)
    ## Check the output file 
    if not os.path.isfile(fname):
        raise FileNotFoundError(fname)
    else:
        ## Filter contigs with viral score >= 0.95 and pvalue <= 0.01
        filtered = "/".join([dvf_dir, ctg_name+"_filtered.fa"])
        predtxt = pd.read_csv(fname, sep="\t", index_col = 0)
        predtxt_filtered = predtxt.loc[(predtxt["score"] >= 0.95) & (predtxt["pvalue"] <= 0.01)]
        indices = predtxt_filtered.index.tolist()
        new_ctg_dict = refine_ctg(ctg_path)

        with open(filtered, "w") as f:
            for idx in indices:
                idx = "_".join(idx.split())
                header = ">"+idx
                seq = new_ctg_dict[header]
                f.write("%s\n%s\n" % (header, seq))
        f.close()

## Concatenate multiple sequences lines of a single contig into one line of sequences
def refine_ctg(ctg_path):
    new_ctg = []
    new_ctg_dict = {}
    idx = 0
    ctg = [i.strip() for i in open(ctg_path)]
    for n, l in enumerate(ctg):
        if n ==0 and l.startswith(">"):
            l = "_".join(l.split())
            new_ctg.append(l)
        elif n != 0 and l.startswith(">"):
            l = "_".join(l.split())
            seq = "".join(ctg[idx+1:n])
            new_ctg.append(seq)
            new_ctg.append(l)
            idx = n
        elif n == (len(ctg)-1):
            seq = "".join(ctg[idx+1:])
            new_ctg.append(seq)
    for h, s in zip(new_ctg[::2], new_ctg[1::2]):
        if not h.startswith(">"):
            raise ValueError(h)
        elif ">" in s:
            raise ValueError(s)
        else:
            new_ctg_dict[h] = s
    return new_ctg_dict

def main(sub_ctg_tsv, sv_dir, threads):
    ## Unifying assembly name format
    for idx in sub_ctg_tsv.index:
        if len(idx.split("/")) == 1:
            ctg_name = idx.rsplit(".", 1)[0]
        else:
            ctg_name = idx.split("/")[0]
    
        cohort = sub_ctg_tsv.loc[idx, "Cohort"]
        ctg_path = sub_ctg_tsv.loc[idx, "Contigs_path"]
        print("Start dvf on %s %s" % (cohort, ctg_name))
        run_dvf(ctg_name, cohort, ctg_path, threads, sv_dir)
        print("Done dvf")

if __name__ == "__main__":
    ## ctg_tsv contains name of cohort and samples, and absolute path for assemblies from MetaSPAdes or MEGAHIT
    ctg_tsv = pd.read_csv(sys.argv[1], sep="\t", index_col = 0)
    sv_dir = sys.argv[2]
    threads = sys.argv[3]
    from_ = int(sys.argv[4])
    interval_ = int(sys.argv[5])
    ## Split queries into chuncks
    sub_ctg_tsv = ctg_tsv.iloc[from_::interval_, :]
    main(sub_ctg_tsv, sv_dir, threads)
