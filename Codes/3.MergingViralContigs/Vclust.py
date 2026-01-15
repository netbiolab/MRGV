#!/usr/bin/env python3

import os, sys, subprocess, pickle
from multiprocessing import Process

### This code is for running Vclust to sample-wise deduplicate viral contigs from DeepVirFinder, VIBRANT, and Phigaro
def run_vclust(t, svdir, sub_keys, contig_dict):
  for n, sub_k in enumerate(sub_keys):
    sample_idx = sub_k.split(".f")[0]
    input = contig_dict[sub_k]
    out1 = "/".join([svdir, sample_idx])
    out2 = "/".join([svdir+"_clu_opts", sample_idx])
    ids = ".".join([sample_idx.rsplit(".", 1)[0], "ids", sample_idx.rsplit(".", 1)[-1]])
    out3 = "/".join([svdir, ids])
    ## Vclust alignment with 0.999 ANI
    align_cmd = "vclust.py align -i %s -o %s --out-ani 0.999" % (input, out1)
    ## Vclust UCLUST with ANI 100% and COV 100%
    cluster_cmd = "vclust.py cluster -i %s -o %s.clus --ids %s --algorithm uclust --metric ani --ani 1.0, --cov 1.0 --out-repr" % (out1, out2, out3)
    _ = subprocess.run([align_cmd], shell=True, capture_output=True)
    _ = subprocess.run([cluster_cmd], shell=True, capture_output=True)
    
if __name__ == "__main__":
  procs = []
  ## contig_dict = {Samples : concatanated contigs from outputs of DeepVirFinder, VIBRANT and Phigaro, ...}
  contig_dict = pickle.load(open(sys.argv[1], "rb"))
  svidr = sys.argv[2]
  threads = int(sys.argv[3])
  from_ = int(sys.argv[4])
  interval_ = int(sys.argv[5])
  contig_keys = list(contig_dict.keys())
  ## Split samples into chunk
  pre_sub_keys = keys[from_::interval_]
  for t in range(threads):
    ## Split samples in chunks into second chunks for multiprocessing
    sub_keys = pre_sub_keys[t::threads]
    proc = Process(target=run_vclust, args=(t, svdir, sub_keys, contig_dict))
    proc.start()
    procs.append(proc)
  for proc in procs:
    proc.join()
  
