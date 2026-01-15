#!/usr/bin/env python3

import os, sys, pickle, subprocess

gbk_in = sys.argv[1]
out = sys.argv[2]
phold_db = sys.argv[3]
threads = sys.argv[4]
pred_dir = sys.argv[5]
phold_cmd = "phold compare -i %s --predictions_dir %s -o %s -t %s -d %s --keep_tmp_files" % (gbk_in, pred_dir, out, threads, phold_db)
subprocess.run([phld_cmd], shell=True)
