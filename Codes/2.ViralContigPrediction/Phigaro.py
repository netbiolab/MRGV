#!/usr/bin/env python3
import os, sys, pickle
import pandas as pd

### This code is for running Phigaro 

def run_phigaro(ctg_name, cohort, ctg_path, threads, sv_dir):
    
    phigaro_dir = "/".join([sv_dir, cohort, ctg_name, "phigaro"])
    if not os.path.isdir(phigaro_dir):
        os.makedirs(phigaro_dir, exist_ok=True)
    new_ctg = refine_ctg(ctg_path)
    new_ctg_name = "/".join([phigaro_dir, "%s_header_revised.fa" % ctg_name])
    out_name = new_ctg_name.rsplit(".", 1)[0]+".phigaro.fasta"
    final_name = out_name.rsplit(".", 1)[0] + ".refined.fa"

    if os.path.isfile(final_name):
        print("File exists %s skipped" % final_name.split("/")[-1])
        return
    if not os.path.isfile(new_ctg_name):
        with open(new_ctg_name, "w") as f:
            for h, s in new_ctg.items():
                f.write("%s\n%s\n" % (h, s))
        f.close()
    if not os.path.isfile(out_name):
        cmd = "phigaro -f %s -p -e gff -o %s --not-open -t %s --no-cleanup --save-fasta -d" % (new_ctg_name, phigaro_dir, threads)
        os.system(cmd)
    if not os.path.isfile(out_name):
        print("No prophages %s" % out_name.split("/")[-1])
        return
    if os.path.isfile(out_name) and not os.path.isfile(final_name):
        rctg = refine_ctg(out_name)
        with open(final_name, "w") as f:
            for h, s in rctg.items():
                f.write("%s\n%s\n" % (h, s))
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
   
    for idx in sub_ctg_tsv.index:
        ## Unifying assembly name format
        if len(idx.split("/")) == 1:
            ctg_name = idx.rsplit(".", 1)[0]
        else:
            ctg_name = idx.split("/")[0]
        ## Run VIBRANT
        cohort = sub_ctg_tsv.loc[idx, "Cohort"]
        ctg_path = sub_ctg_tsv.loc[idx, "Contigs_path"]
        print("Start Phigaro on %s %s" % (cohort, ctg_name))
        run_phigaro(ctg_name, cohort, ctg_path, threads, sv_dir)
        print("Done Phigaro")

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
