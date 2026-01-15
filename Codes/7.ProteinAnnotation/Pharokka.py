#!/usr/bin/env python3

import os, sys, subprocess

query_fasta = sys.argv[1]
out = sys.argv[2]
pharokka_db = sys.argv[3]
threads = sys.argv[4]
pharokka_cmd = "pharokka.py -i %s -o %s -d %s -t %s -g phanotate" % (query, out, pharokka_db, threads)
subprocess.run([pharokka_cmd], shell=True)
