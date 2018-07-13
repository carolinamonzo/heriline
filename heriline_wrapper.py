#! /usr/bin/python3.5

# Carolina Monzo - 2018-07-13

# Final wrapper for automatic running of heriline

import argparse
import os
import re
import pandas as pd
import glob
import logging
import coloredlogs
import json
import pprint
import datetime
import subprocess

coloredlogs.install()


## RUN UTILS.PY

from utils import *

args = parseArguments()

logger = setLogger(setMode(args), "__main__")

config = load_config(args)

check_args(args, config, logger)

config = create_directories(config, logger)

get_software_versions(config, logger)

# This must be the last thing to do each step
export_config(config, logger)

## RUN MERGE FASTQ_FILES.PY

import merge_fastq
# Running only functions of interest

logger.info("MERGING FASTQ FILES using zcat\n")

fof = merge_fastq.read_input_fof(config)

df_fastq_sorted = merge_fastq.fastq_dataframe(config, fof)

cmd_sh = merge_fastq.cmd_zcat_fastq(config, df_fastq_sorted)

merge_fastq.run_parallel(config, cmd_sh)

merge_fastq.write_output_fof(config)

# Document on config file

config["steps"]["step1_merge_fastq"] = "DONE_{}".format(output_fof)

# This must be the last thing to do each step
export_config(config, logger)


## RUN TRIMM FASTQ_FILES.PY

import trimm_fastq

# Running only functions of interest
logger.info("TRIMMING FASTQ FILES using seqtk\n")

fof = trimm_fastq.read_input_fof(config)

cmd_sh = trimm_fastq.cmd_trimm_fastq(config, fof)

trimm_fastq.run_parallel(config, cmd_sh)

trimm_fastq.write_output_fof(config)

config["steps"]["step2_trimm_fastq"] = "DONE_{}".format(output_fof)

# This must be the last thing to do each step
export_config(config, logger)


## RUN MAP.PY

import map_fastq

logger.info("MAPPING using BWA-MEM\n")

fof = map_fastq.read_input_fof(config)

df_fastq_sorted = map_fastq.merged_fastq_dataframe(config, fof)

cmd_sh = map_fastq.cmd_map_fastq(config, df_fastq_sorted)

map_fastq.run_parallel(config, cmd_sh)

map_fastq.write_output_fof(config)

config["steps"]["step3_map"] = "DONE_{}".format(output_fof)

# This must be the last thing to do each step

export_config(config, logger)


## RUN MARK DUPLICATES.PY

import mark_duplicates

logger.info("MARKING DUPLICATES using picard\n")

fof = mark_duplicates.read_input_fof(config)

cmd_sh = mark_duplicates.cmd_mark_duplicates(config, fof)

mark_duplicates.run_parallel(config, cmd_sh)

mark_duplicates.write_output_fof(config)

config["steps"]["step4_mark_duplicates"] = "DONE_{}".format(output_fof)

export_config(config, logger)


## RUN SPLIT CHR.PY

import split_chr

logger.info("SPLITTING BAM FILES PER CHROMOSOME using samtools\n")

fof = split_chr.read_input_fof(config)

cmd_sh = split_chr.cmd_split_bam(config, fof)

split_chr.run_parallel(config, cmd_sh)

split_chr.write_output_fof(config)

config["steps"]["step5_split_chr"] = "DONE_{}".format(output_fof)

export_config(config, logger)

## RUN MERGE BAM FILES PER SAMPLE AND CHROMOSOME

import merge_bams

logger.info("MERGING BAM FILES PER SAMPLE AND CHROMOSOME using samtools\n")

fof = merge_bams.read_input_fof(config)

cmd_sh = merge_bams.cmd_merge_bam(config, fof)

merge_bams.run_parallel(config, cmd_sh)

merge_bams.write_output_fof(config)

config["steps"]["step6_merge_bam"] = "DONE_{}".format(output_fof)

export_config(config, logger)

## RUN VARIANT CALLING

import variant_calling

fof = variant_calling.read_input_fof(config)

cmd_sh = variant_calling.cmd_call_variants(config, fof)

variant_calling.run_parallel(config, cmd_sh)

variant_calling.write_output_fof(config)

config["steps"]["step7_call_variants"] = "DONE_{}".format(output_fof)

export_config(config, logger)

pprint.pprint(config)





