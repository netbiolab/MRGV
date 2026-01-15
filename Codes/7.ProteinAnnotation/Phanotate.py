#!/usr/bin/env python3

import os, sys, subprocess

query_fasta = sys.argv[1]
out = sys.argv[2]
phanotate_cmd = "phanotate.py -o %s -f genbank"
subprocess.run([phanotate_cmd], shell=True)
