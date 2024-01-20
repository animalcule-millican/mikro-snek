#!/usr/bin/env python3
import snakemake
import os
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="Run snakemake workflow")
    parser.add_argument('-s', "--snakefile", type=str, help="path to snakefile", default = "/home/glbrc.org/millican/repos/mikro-snek/workflow/Snakefile")
    #parser.add_argument("--configfile", type=str, help="path to configfile", default = "/home/glbrc.org/millican/repos/mikro-snek/workflow/config.yml")
    parser.add_argument('-i', "--sequence_directory", type=str, help="Path to directory with sequence files", required=True)
    parser.add_argument('-c', "--chunks", type=int, help="Number of chunks to split sequence files into", default=1)
    parser.add_argument("--cores", help="number of cores", default='all')
    parser.add_argument("--dryrun", help="Testing workflow with a dry run. Will not execute any rules.", action="store_true")
    parser.add_argument("--printrulegraph", help="printrulegraph", action="store_true")
    parser.add_argument("--keepgoing", help="keepgoing", action="store_true", default=True)
    parser.add_argument("--latency_wait", type=str, help="latency_wait", default=120)
    parser.add_argument("--use_conda", help="use_conda", action="store_true")
    parser.add_argument("--conda_prefix", type=str, help="Path to location of Conda/Mamba", default=os.environ['CONDA_PREFIX'])
    parser.add_argument("--conda_cleanup_pkgs", help="conda_cleanup_pkgs", action="store_true")
    parser.add_argument("--profile", type=str, help="Profile for running on cluster. 'slurm' or 'HTCondor'. Default is none.", default = "HTCondor")
    args = parser.parse_args()
    return args

def build_config(args):
    config = {}
    config["sequence_files"] = args.sequence_directory
    config["split_chunks"] = args.chunks
    config['project_directory'] = os.path.abspath(os.path.join(os.path.dirname(args.snakefile), os.pardir))
    return config

def main():
    args = parse_args()
    config = build_config(args)
    workdir = os.path.dirname(args.snakefile)

    os.system("~/mambaforge/bin/snakemake --snakefile {0} --configfile {1} --directory {2} --cores {3} --jobs {8} --dryrun --stats {4} --printrulegraph --keepgoing --latency-wait {5} --use-conda --conda-prefix {6} --conda-cleanup-pkgs --profile {7}".format(args.snakefile, args.configfile, workdir, args.cores, os.path.join(workdir, "stats.txt"), args.latency_wait, args.conda_prefix, args.profile))
    #snakemake.snakemake(snakefile = args.snakefile, config = config, workdir = workdir, cores = args.cores, dryrun = args.dryrun, stats = os.path.join(workdir, "stats.txt"), printrulegraph = args.printrulegraph, keepgoing = args.keepgoing, latency_wait = args.latency_wait, use_conda = args.use_conda, conda_prefix = args.conda_prefix, conda_cleanup_pkgs = args.conda_cleanup_pkgs)

path = subprocess.check_output(args="which snakemake", shell=True, text=True).strip()

--snakefile
--configfile
--config
--cores
--jobs
--directory
--conda-base-path
--rulegraph
--unlock 
--latency-wait
--use-conda
--conda-prefix
--conda-cleanup-pkgs
if __name__ == "__main__":
    main()