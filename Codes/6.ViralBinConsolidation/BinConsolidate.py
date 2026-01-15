#! /usr/bin/env python3

import os, sys, pickle, fastcluster, subprocess
import pandas as pd
import numpy as np
from itertools import combinations
from scipy.cluster.hierarchy import cut_tree
from multiprocessing import Process

### This code is for sample-wise bin-consolidation using UPGMA clustering

def comb2(N):
    return int(N*(N-1)/2)

def get_darray(ani_info, glist):
    g2idx = {i : n for n, i in enumerate(glist)}
    darray = [-1] * comb2(len(glist))
    for i in combinations(glist, 2):
        g1, g2 = i
        dist = 1 - get_ani(ani_info, g1, g2)
        i = g2idx[g1]
        j = g2idx[g2]
        if not (i < j):
            continue
        else:
            da_i = comb2(len(glist)) - comb2((len(glist)-i)) + (j-i-1)
            darray[da_i] = dist
    if -1 in darray:
        raise ValueError("negative value in darray")
    else:
        return darray

def get_cluster(darray, glist, cutoff):
    z = fastcluster.linkage(darray, method = "average", preserve_input=False)
    cutree = cut_tree(z, height = cutoff)
    clu_dict = dict()
    glist = list(glist)
    for idx in range(len(glist)):
        clu = cutree[idx][0]
        if not clu in clu_dict.keys():
            clu_dict[clu] = set([glist[idx]])
        else:
            clu_dict[clu].add(glist[idx])
    return clu_dict

def compute_ani(g1, g2, g1_ln, g2_ln):
    if g1 == g2:
        return 1
    else:
        ovl_ln = len("".join(g1 & g2))
        ani = max(ovl_ln/g1_ln, ovl_ln/g2_ln)
        return ani

def get_ani(ani_info, g1, g2):
    if g1 == g2:
        return 1.0
    try:
        ani = ani_info[g1][g2]
    except KeyError:
        try:
            ani = ani_info[g2][g1]
        except KeyError:
            ani = 0.0
    return ani


def update_ani(cs, bidx1, bidx2, ani_info, ani):
    if not bidx1 in ani_info.keys():
        if not bidx2 in ani_info.keys():
            ani_info[bidx1] = {bidx2 : ani}
        else:
            if not bidx1 in ani_info[bidx2].keys():
                ani_info[bidx2][bidx1] = ani
    else:
        if not bidx2 in ani_info[bidx1].keys():
            ani_info[bidx1][bidx2] = ani
    return ani_info

def cherry_pick(score_d, clu_dict):
    rep_dict = {}
    for clu, members in clu_dict.items():
        members = list(members)
        if len(members) == 1:
            rep_dict[members[0]] = set()
            continue
        else:
            scores = [score_d[i] for i in members]
            max_idx = scores.index(max(scores))
            rep_g = members[max_idx]
            rep_dict[rep_g] = set(members) - set([rep_g])
    return rep_dict       
            
def check_redunancy(rep_dict, sub_hsd):
    for rep1 in rep_dict.keys():
        for rep2 in rep_dict.keys():
            if rep1 == rep2:
                continue
            else:
                ovl = set(sub_hsd[rep1].values()) & set(sub_hsd[rep2].values())
                if len(ovl) > 0:
                    return False
    return True


def bin_consolidation(t, score_d, sub_keys, hsd):
    t_clu_dict = {}
    t_rep_dict = {}
    t_check_idx = True
    for cs in sub_keys:
        v = hsd[cs]
        ani_info = {}
        glist = list(v.keys())
        if len(glist) == 1:
            if not cs in t_rep_dict.keys():
                t_rep_dict[cs] = {glist[0] : set()}
            if not cs in t_clu_dict.keys():
                t_clu_dict[cs] = {glist[0] : set([glist[0]])}
            continue
        for bidx1, v2 in v.items():
            seqs1 = set(v2.values())
            seq1_ln = len("".join(seqs1))
            for bidx2, v2 in v.items():
                if bidx1 == bidx2:
                    continue
                seqs2 = set(v2.values())
                seq2_ln = len("".join(seqs2))
                ani = compute_ani(seqs1, seqs2, seq1_ln, seq2_ln)
                ani_info = update_ani(cs, bidx1, bidx2, ani_info, ani)
        da = get_darray(ani_info, glist)

        clu_dict = get_cluster(da, glist, 1.0)
        rep_dict = cherry_pick(score_d, clu_dict)
        check_idx = check_redunancy(rep_dict, v)
        if not check_idx:
            t_check_idx = False

        if not cs in t_rep_dict.keys():
            t_rep_dict[cs] = rep_dict
        else:
            for k3, v3 in rep_dict.items():
                t_rep_dict[cs][k3] = v3
        if not cs in t_clu_dict.keys():
            t_clu_dict[cs] = clu_dict
        else:
            for k3, v3 in clu_dict.items():
                t_clu_dict[cs][k3] = v3
    if t_check_idx:
        return t_clu_dict, t_rep_dict
    else:
        print(t_check_idx)
        raise ValueError("Redundancy found")
        
