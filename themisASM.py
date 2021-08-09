#!/usr/bin/env python3

import sys
import argparse
import os
import logging
import subprocess as sp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
LOG = logging.getLogger()

version = "0.0.1"

def parse_user_input():
    parser = argparse.ArgumentParser(
            description = "Themis-ASM: Compare multiple genome assemblies with short-read alignment data" + version
            )
    parser.add_argument('-a', '--assembly',
                        help="An assembly fasta file to compare. Provide the full path to the file",
                        action="append", default=[]
                        )
    parser.add_argument('-n', '--name',
                        help="A name for the assembly file. Add a name in the order in which the -a option is specified.",
                        action="append", default=[]
                        )
    parser.add_argument('-b', '--busco',
                        help="The name of the BUSCO database to use for assessment",
                        type=str, required=True
                        )
    parser.add_argument('-f', '--fastq',
                        help="Paired short read fastqs separated by commas.",
                        action="append", default=[]
                        )
    parser.add_argument('-s', '--sample',
                        help="A sample name for the fastqs provided. Add in the order in which the fastq files are specified",
                        action="append", default=[]
                        )
    parser.add_argument('-c', '--cluster_string',
                        help="A short formatted string for instructions on how to submit to your cluster [Default: do not submit to cluster]",
                        type=str, default=None
                        )
    parser.add_argument('-j', '--jobs',
                        help="Maximum concurrent jobs to run [Default: 20]",
                        type=int, default=20
                        )
    parser.add_argument('-r', '--resume',
                        help="Resume a previously terminated job [Flag]",
                        action="store_true", default=False
                        )
    parser.add_argument('-u', '--unlock',
                        help="Unlock a previously terminated job directory [Flag]",
                        action="store_true", default=False
                        )

    return parser.parse_args(), parser

def main(args, parser):
    # Sanity checks
    try:
        validate(args)
    except RuntimeError as inst:
        LOG.error(inst.args)
        sys.exit()

    # Create the config file
    try:
        createConfig(args)
    except RuntimeError as inst:
        LOG.error(inst.args)
        sys.exit()

    # Generate the command
    cmd = createSNKCmd(args)

    # Now attempt to run it!
    try:
        LOG.info("Trying to run the snakemake job in this scrpt... wish me luck!")
        sp.run(cmd, shell=True, check=True)
    except Exception as inst:
        LOG.info("We ran into an error! Please try to resume using this script or run the snake command manually")
        LOG.error(inst.args)
        sys.exit()

    LOG.info("I think we're finished! Check that zip file!")

def createSNKCmd(args):
    cmd = f'snakemake -s {SCRIPT_DIR}/themisSnakefile --jobs {args.jobs} -p --use-conda'

    if args.cluster_string != None:
        LOG.info("Detected cluster submission information! Running cluster jobs")
        cmd += f' --cluster-config {SCRIPT_DIR}/cluster.json --cluster "{args.cluster_string}"'

    LOG.info("Created snakemake submission string\nIf you want to run this on your own, please use the following command:\n" + cmd)
    return cmd

def validate(args):
    if len(args.name) != len(args.assembly):
        raise RuntimeError(f'Please enter the same count of assembly fastas {len(args.assembly)} as your names {len(args.name)}')
    if len(args.sample) != len(args.fastq):
        raise RuntimeError(f'Please enter the same count of read fastqs {len(args.fastq)} as your sample names {len(args.sample)}')

    rdir = os.path.join(os.getcwd(), '.snakemake')
    if os.path.isdir(rdir) and not args.resume:
        raise RuntimeError("It looks like you've already run snakemake in this directory -- please remove or resume your job!")
    ldir = os.path.join(rdir, 'locks')
    if os.path.isdir(ldir) and os.listdir(ldir):
        raise RuntimeError("It looks like the directory is locked by Snakemake. Please run an unlock job first with --unlock!")

def createConfig(args):
    config = f'{{\n  "buscoLineage" : "{args.busco}",\n  "assembly" : {{\n';
    for aname, afile in zip(args.name, args.assembly):
        config += f'    "{aname}" : "{afile}",\n'
    config += "  },\n  \"samples\" : {\n"

    for sname, fqstr in zip(args.sample, args.fastq):
        fqseg = fqstr.split(',')
        if len(fqseg) != 2:
            raise RuntimeError(f'It looks like you entered an unpaired fastq set! Please add both reads as a comma delmited entry!\n')
        config += f'    "{sname}" : [\n      "{fqseg[0]}",\n      "{fqseg[1]}"\n    ],\n'

    config += "  }\n}\n"
    with open("default.json", 'w') as co:
        co.write(config)
    LOG.info("Created new configuration file: default.json")

if __name__ == "__main__":
    args, parser = parse_user_input()
    main(args, parser)
