#! usr/bin/python3.6

# 2018-03-28, Carolina Monzo

# Merge bams

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
    parser = argparse.ArgumentParser(description = "Merge fastq files")

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
    '''
    Load bam files
    '''

    # Get fof file
    fof_list = []

    for file in glob.glob(config["paths"]["fof_files"] + "splitted_bams_*.fof"):
        fof_list.append(file)

    fof_pwd = os.path.normpath(sorted(fof_list)[-1])

    # Read fof file

    with open(fof_pwd, "r") as fi:
        fof = fi.read().splitlines()
    return(fof)

def cmd_merge_bam(config, fof):
    '''
    Create cmd file
    '''

    cmd_sh = datetime.datetime.now().strftime("cmd_merge_bam_%Y%m%d_%H-%M-%S.sh")
    chrom = []

    chr_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]

    with open(config["paths"]["cmd_files"] + cmd_sh, "a") as cmd_file:

        for i in chr_list:
            
            for fi in fof:
                if ("_chr" + i + ".bam") in fi:
                    chrom.append(fi)
            files = " ".join(sorted(chrom))
            cmd_str = "samtools merge -c {}all_chr{}_merged.bam {}; samtools index {}all_chr{}_merged.bam".format(config["paths"]["mapping_merged"], i, files, config["paths"]["mapping_merged"], i)

            cmd_file.write(cmd_str + "\n")
            chrom = []

    print("[INFO]: CMD_FILE - {}{}".format(config["paths"]["cmd_files"], cmd_sh))
    
    return(cmd_sh)

def run_parallel(config, cmd_sh):
    '''

    '''

    log_str = datetime.datetime.now().strftime("merge_bam_%Y%m%d_%H-%M-%S.log")

    cmd = "parallel --joblog {}{} -j 18 :::: {}{}".format(config["paths"]["mapping_merged"], log_str, config["paths"]["cmd_files"], cmd_sh)

    print("[CMD] " + cmd)

    # Execute command
    subprocess.call(cmd + " 2> /dev/null", shell = True)

def write_output_fof(config):
    '''

    '''

    fof = datetime.datetime.now().strftime("merged_bams_chr_%Y%m%d_%H-%M-%S.fof")

    # Create fof file
    cmd_fof = 'find {} -name "all_chr*_merged.bam" > {}{}'.format(config["paths"]["mapping_merged"], config["paths"]["fof_files"], fof)

    # Write command on the command line

    print("[CMD] " + cmd_fof)
    print("[INFO]: FOF_FILE - {}{}".format(config["paths"]["fof_files"], fof))

    subprocess.call(cmd_fof, shell = True)

def main():

    args = parseArguments()

    config = get_config(args.project_path)

    fof = read_input_fof(config)

    cmd_sh = cmd_merge_bam(config, fof)

    run_parallel(config, cmd_sh)

    write_output_fof(config)

if __name__ == '__main__':

    main()
