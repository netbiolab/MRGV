#!/usr/bin/env python3

import os, sys, pickle, subprocess

gbk_in = sys.argv[1]
out = sys.argv[2]
phold_db = sys.argv[3]
threads = sys.argv[4]
phold_cmd = "phold predict -i %s -o %s -t %s -d %s" % (gbk_in, out, threads, phold_db)
subprocess.run([phld_cmd], shell=True)
