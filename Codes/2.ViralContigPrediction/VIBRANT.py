#!/usr/bin/env python3
import os, sys, pickle
import pandas as pd

### This code for running VIBRANT and annotating lifestyles in contig headers

def run_vibrant(ctg_name, cohort, ctg_path, threads, sv_dir):
    vibrant_dir = "/".join([sv_dir, cohort, ctg_name, "vibrant"])
    if not os.path.isdir(vibrant_dir):
        os.makedirs(vibrant_dir, exist_ok=True)
    cmd = "python3 /home/hjkim/inetbiolab/Phage/HRGV_2nd/VIBRANT/VIBRANT_run.py -i %s -folder %s -t %s" % (ctg_path, vibrant_dir, threads)
    out_path = "/".join([vibrant_dir, "VIBRANT_%s" % ctg_path.split("/")[-1].rsplit(".", 1)[0], "VIBRANT_phages_%s" % ctg_path.split("/")[-1].rsplit(".", 1)[0]])
    out_lytic = "/".join([out_path, ctg_path.split("/")[-1].rsplit(".", 1)[0]+".phages_lytic.fna"])
    out_lysogenic = "/".join([out_path, ctg_path.split("/")[-1].rsplit(".", 1)[0]+".phages_lysogenic.fna"])
    out_combine = "/".join([out_path, ctg_path.split("/")[-1].rsplit(".", 1)[0]+".phages_combined.fna"])
    if not os.path.isdir(out_path):
        print("output dir %s not exist" % out_path)
        
        os.system(cmd)

    else:
        print("Dir %s exist so, skipped" % out_path)

    if os.path.isfile("/".join([out_path, "%s_annotated.fa" % ctg_name])):
        print("Final file %s exist, so skipped" % "/".join([out_path, "%s_annotated.fa" % ctg_name]))
        return
    if os.path.isfile(out_combine) and not os.path.isfile("/".join([out_path, "%s_annotated.fa" % ctg_name])):
        print("ctg %s will be annotated" % ctg_name)
        combine = []
        if os.path.isfile(out_lytic) and os.path.getsize(out_lytic) != 0:
            print("starts process lytic ctgs")
            new_lytic = refine_ctg(out_lytic)
            for h, s in new_lytic.items():
                nh = "@".join([h, "lytic"])
                combine.append(nh)
                combine.append(s)
        if os.path.isfile(out_lysogenic) and os.path.getsize(out_lysogenic) != 0:
            print("starts process lysogenic ctgs")
            new_lysogenic = refine_ctg(out_lysogenic)
            for h, s in new_lysogenic.items():
                nh = "@".join([h, "lysogenic"])
                combine.append(nh)
                combine.append(s)
        if not len(combine) == 0:
            with open("/".join([out_path, "%s_annotated.fa" % ctg_name]), "w") as f:
                for i in combine:
                    f.write("%s\n" % i)
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
        print("Start VIBRANT on %s %s" % (cohort, ctg_name))
        run_vibrant(ctg_name, cohort, ctg_path, threads, sv_dir)
        print("Done VIBRANT")

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
