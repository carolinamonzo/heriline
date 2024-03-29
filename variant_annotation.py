#! usr/bin/python3.5

# 2018-04-12, Carolina Monzo

# Variant annotation

import os
import argparse
import json
import datetime
import glob
import subprocess

def parseArguments():
    '''
    Function to parse arguments
    Input: path to project
    Output: parsed arguments from command line
    '''

    # Create argument parser class
    parser = argparse.ArgumentParser(description = "Annotate variants in VCF split by chromosome")

    # Define arguments to parse
    parser.add_argument("--project_path", "-p", required = False, type = str, help = "Argument to set the path to the project")

    # Call for arguments
    args = parser.parse_args()

    return(args)

def get_config(project_path):
    """
    Load config dictionary
    """

    config_list = []

    os.chdir(project_path)
    for file in glob.glob("config_*.json"):
        config_list.append(file)

    config_pwd = project_path + sorted(config_list)[-1]

    with open(config_pwd, "r") as jconfig:
        config = json.load(jconfig)

    return(config)


def read_input_fof(config):
    """

    """

    # Get fof file
    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + "decompose_chr_*.fof"):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    # Read fof file
    with open(fof_pwd, "r") as fi:

        fof = fi.read().splitlines()

    return(fof)


def cmd_annotate_vcf(config, fof):
    """

    """

    cmd_sh = datetime.datetime.now().strftime("cmd_annotate_variants_%Y%m%d_%H-%M-%S.sh")

    with open(config["paths"]["cmd_files"] + cmd_sh, "a") as cmd_file:
        for vcf in fof:
            vcf_name = vcf.split("/")[-1].split(".")[0]
            cmd_str = "vep -i {} --format vcf --vcf -o {}{}_vep.vcf --cache_version 92 --dir_cache /nfs/qnapugdg8tb3a/nfs_databases/vep/vep --assembly GRCh37 --offline --force_overwrite --everything --sift s --polyphen s".format(vcf, config["paths"]["variant_annotation"], vcf_name)


            cmd_file.write(cmd_str + "\n")

    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))

    return(cmd_sh)


def run_parallel(config, cmd_sh):
    """

    """

    log_str = datetime.datetime.now().strftime("annotate_variants_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 18 :::: {}{}".format(config["paths"]["variant_annotation"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD]: " + cmd)

    subprocess.call(cmd + " 2> /dev/null", shell = True)

def write_output_fof(config):
    """

    """

    fof = datetime.datetime.now().strftime("variant_annotation_chr_%Y%m%d_%H-%M-%S.fof")

    cmd_fof = 'find {} -name "all_chr*_merged-decomp-norm_vep.vcf" > {}{}'.format(config["paths"]["variant_annotation"], config["paths"]["fof_files"], fof)

    print("[CMD] " + cmd_fof)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))

    subprocess.call(cmd_fof, shell = True)

def main():
    """

    """

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config)

    cmd_sh = cmd_annotate_vcf(config, fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

if __name__ == "__main__":
    main()









