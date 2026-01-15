#!/usr/bin/env python3

import os, sys, subprocess

query = sys.argv[1]
out = sys.argv[2]
min_id = sys.argv[3]
mmseq_cmd = "mmseqs easy-linclust %s %s tmp --min-seq-id %s" % (query, out, min_id)
