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

# This must be the last thing to do here
#export_config(config, logger)

## RUN MERGE FASTQ_FILES.PY

import merge_fastq
# Running only functions of interest

logger.info("MERGING FASTQ FILES\n")

fof = merge_fastq.read_input_fof(config)

df_fastq_sorted = merge_fastq.fastq_dataframe(config, fof)

cmd_sh = merge_fastq.cmd_zcat_fastq(config, df_fastq_sorted)

merge_fastq.run_parallel(config, cmd_sh)

merge_fastq.write_output_fof(config)

# Document on config file

config["steps"]["step1_merge_fastq"] = "DONE_{}".format(output_fof)


pprint.pprint(config)
