import os, sys
import pandas as pd

def get_uniqueness_min(meta, dm, rank):
    uni_list = []
    for i in dm.index:
        target = dm.loc[dm.index.isin([i])]
        cage = meta.loc[i, "Cage"]
        cage_mates = meta.loc[meta["Cage"] == cage].index.tolist()
        others_all = list(set(dm.columns.tolist()) - set([i]))
        others_nocage = list(set(others_all) - set(cage_mates))
        others_cage = list((set(dm.columns.tolist()) & set(cage_mates)) - set([i]))
        target_w_others = target.loc[:, others_all]
        target_wo_cagemates = target.loc[:, others_nocage]
        target_w_cagemates = target.loc[:, others_cage]
        min_w_others = target_w_others.min(axis=1).values[0]
        min_wo_cagemates = target_wo_cagemates.min(axis=1).values[0]
        min_w_cagemates = target_w_cagemates.min(axis=1).values[0]
        uni_list.append([i, rank, min_w_others, min_wo_cagemates, min_w_cagemates])
    uni_df = pd.DataFrame(uni_list, columns=["SampleID", "Rank", "With CageMates", "Without CageMates", "Within CageMates"])
    uni_df.set_index("SampleID", inplace=True)
    return uni_df

if __name__ == "__main__":
  meta = pd.rad_csv(sys.argv[1], sep="\t", index_col =0)
  dist_matrix = pd.read_csv(sys.argv[2], sep="\t", index_col = 0)
  viral_tax_rank = sys.argv[3]
  uni_df = get_uniqueness_min(meta, dist_matrix, vrial_tax_rank)
  uni_df.to_csv("%s_uni_df.tsv" % viral_tax_rank, sep="\t")
