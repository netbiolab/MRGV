import os, sys, pickle, subprocess

def run_coverm(bam_list, sub_list, threads):
    for n, i in enumerate(sub_list):
        print("Start on %s" % (i.split("/")[-1]))
        svpath = i.rsplit("/", 1)[0]
        svname1 = i.split("/")[-1].replace(".bam.sorted", ".covered_bases")
        svname2 = i.split("/")[-1].replace(".bam.sorted", ".aligned_reads_counts")
        svpath1 = "/".join([svpath, svname1])
        svpath2 = "/".join([svpath, svname2])
        cmd1 = "coverm contig -b %s -m covered_bases -t %s -o %s" % (i, threads, svpath1)
        cmd2 = "coverm contig -b %s -m count -t %s -o %s" % (i, threads, svpath2)
        cmd_all = " && ".join([cmd1, cmd2])
        subprocess.run([cmd_all], shell=True)
        #print(cmd_all)
    return None

if __name__ == "__main__":
    bam_list = [i.strip() for i in open(sys.argv[1])]
    from_ = int(sys.argv[2])
    interval_ = int(sys.argv[3])
    threads = int(sys.argv[4])
    run_coverm(bam_list, bam_list[from_::interval_], threads)   
