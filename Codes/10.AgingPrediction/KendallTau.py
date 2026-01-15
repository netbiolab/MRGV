import os, sys, pickle, subprocess
from multiprocessing import Process
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from itertools import combinations_with_replacement
import warnings
from scipy.stats import kendalltau
warnings.filterwarnings("ignore")

def get_kendall(t, sub_comb, df, svdir):
    print("start sub %s in process %s" % (len(sub_comb), t))
    taud = {}
    pvald = {}
    for n, (i, j) in enumerate(sub_comb):
        x = df.loc[i]
        x_nz = x.loc[x!=0]
        y = df.loc[j]
        y_nz = y.loc[y!=0]
        nz_union = list(set(x_nz.index.tolist() + y_nz.index.tolist()))
        x_nz_u = x.loc[nz_union]
        y_nz_u = y.loc[nz_union]
        tau, p = kendalltau(x_nz_u, y_nz_u)

        if not i in taud.keys():
            if not j in taud.keys():
                taud[i] = {j : tau}
                pvald[i] = {j: p}
            else:
                if not i in taud[j].keys():
                    taud[j][i] = tau
                    pvald[j][i] = p
        else:
            if not j in taud[i].keys():
                taud[i][j] = tau
                pvald[i][j] = p
        if (n+1) % 10000 == 0:
            print("Done sub %s in process %s" % ((n+1), t))
    dd = [taud, pvald]
    pickle.dump(dd, open("%s/tau_pval_d_part%s.pkl" % (svdir, t), "wb"))
    print("Done all in process %s" % t)

def get_tau_pval_d(dir, svdir, prefix):
    taud = {}
    pvald = {}
    for i in os.listdir(dir):
        _1st = "/".join([dir, i])
        tmp = pickle.load(open(_1st, "rb"))
        for n, d in enumerate(tmp):
            if n == 0:
                tod = taud
            else:
                tod = pvald
            for k1, v1 in d.items():
                for k2, v2 in v1.items():
                    if not k1 in tod.keys():
                        if not k2 in tod.keys():
                            tod[k1] = {k2:v2}
                        else:
                            if not k1 in tod[k2].keys():
                                tod[k2][k1] = v2
                    else:
                        if not k2 in tod[k1].keys():
                            tod[k1][k2] = v2
    print("Done to collect tau and pvalue dictionary")

    tau_list = []
    pval_list = []
    for k, v in taud.items():
        for k2, v2 in v.items():
            tau_list.append([k, k2, v2])
            tau_list.append([k2, k, v2])
    for k, v in pvald.items():
        for k2, v2 in v.items():
            pval_list.append([k, k2, v2])
            pval_list.append([k2, k, v2])
    tau_df = pd.DataFrame(tau_list)
    pval_list = pd.DataFrame(pval_list)
    tau_df = tau_df.set_index([0,1])
    pval_df = pval_list.set_index([0, 1])
    tau_df = tau_df.loc[~tau_df.index.duplicated(keep="first")]
    pval_df = pval_df.loc[~pval_df.index.duplicated(keep="first")]
    tau_df = tau_df.unstack()
    tau_df.columns= [i[-1] for i in tau_df.columns]
    tau_df = tau_df[tau_df.index]
    pval_df = pval_df.unstack()
    pval_df.columns= [i[-1] for i in pval_df.columns]
    pval_df = pval_df[pval_df.index]
    tau_dm = tau_df.apply(lambda x: (1-x)/2)
    pval_df.to_csv("%s/%s_kendall_pvaldf.tsv" % (svdir, prefix), sep="\t")
    tau_df.to_csv("%s/%s_kendall_tau_coef.tsv" % (svdir, prefix), sep="\t")
    tau_dm.to_csv("%s/%s_kendall_tau_dm.tsv" % (svdir, prefix), sep="\t")
    return tau_df, tau_dm


if __name__ == "__main__":
    svdir_tmp = sys.argv[1]
    svdir_main = sys.argv[2]
    threads = sys.argv[3]
    abundance_df = sys.argv[4]
    indicator = sys.argv[5]
    tau_dm_dict = {}
    for idx, df in zip([indicator],  [abundance_df]):
        do_sample_comb_wr = list(combinations_with_replacement(df.index.tolist(), 2))
        print("Start %s kendall tau dm calculation" % idx)
        procs = []
        for t in range(n_threads):
            sub_comb = do_sample_comb_wr[t::n_threads]
            proc = Process(target=get_kendall, args=(t, sub_comb, df, svdir_tmp))
            procs.append(proc)
            proc.start()
        for proc in procs:
            proc.join()
        print("Start to collect kendall tau dm on %s" % idx)
        tau_dm_dict[idx] = get_tau_pval_d(svdir_tmp, svdir, idx)
        subprocess.run(["rm %s/*" % svdir_tmp], shell=True)
        print("Delete tmp files in %s" % svdir_tmp)
    pickle.dump(tau_dm_dict, open("/".join([svdir_main, "tau_dm.pkl", "wb"))
    print("Done all kendall tau dm calculations")

