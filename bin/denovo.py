#!/usr/bin/env python3

import argparse
import glob
import os
import sys
import subprocess
from datetime import datetime
import time


SINGULARITY_LIB_PATH = 'library://lakinsm/default/slss-denovo:alpha1'
NANOPORE_TEMP_FILENAME = 'nanopore_concat.fastq'
SLURM_CONFIG = 'slss_hpc_slurm.config'


parser = argparse.ArgumentParser('denovo.py')
parser.add_argument('-r1', '--forward', required=True, type=str,
                    help='Forward, single- or paired-end Illumina FASTQ or Nanopore fastq_pass directory')
parser.add_argument('-r2', '--reverse', default=None, required=False, type=str,
                    help='Reverse paired-end Illumina FASTQ (optional)')
parser.add_argument('-l', '--long_reads', default=None, required=False, type=str,
                    help='Nanopore fastq_pass directory if hybrid assembly (must have -r1, -r2, and -l specified, optional)')
parser.add_argument('--host_reference', default=None, required=False, type=str,
                    help='Host reference sequence FASTA for host depletion (optional)')
parser.add_argument('-w', '--work_dir', default='/tmp', type=str,
                    help='Working directory for temporary files (must be shared if using SLURM/HPC)')
parser.add_argument('-m', '--metagenomic', default=False, action='store_true',
                    help='Flag for metagenomic input file (cannot be used with hybrid assembly)')
parser.add_argument('-s', '--singularity', type=str, default=None,
                    help='Path to Singularity container if other than default (pulls from cloud if this argument isn\'t used)')
parser.add_argument('-t', '--threads', default=1, type=int,
                    help='Number of threads to use')
parser.add_argument('-f', '--force', action='store_true', default=False,
                    help='Overwrite intermediate Nanopore concatenated files.  If not used, will cache Nanopore concatenation.')
parser.add_argument('--slurm', action='store_true', default=False,
                    help='Flag for use of SLURM for HPC clusters.  Modify {} to change cluster options.'.format(
						SLURM_CONFIG
                    ))


def concatenate_nanopore_reads(input_dir, temp_dir):
	with open(temp_dir + '/' + NANOPORE_TEMP_FILENAME, 'w') as out:
		fastq_files = glob.glob(input_dir + '/*.fastq')
		sys.stdout.write('\n{} Nanopore FASTQ files detected in directory {}, concatenating...\n\n'.format(
			len(fastq_files),
			input_dir
		))
		for nanopore_fastq in fastq_files:
			with open(nanopore_fastq, 'r') as f:
				line = f.readline()
				while line:
					out.write(line)
					line = f.readline()


def validate_input_args(arg_obj):
	if arg_obj.reverse:
		if not arg_obj.forward:
			sys.stderr.write('\nError: use of paired-end Illumina reads requires -r1 and -r2 to be specified.\n')
			raise ValueError

	if arg_obj.long_reads:
		if not (arg_obj.forward and arg_obj.reverse):
			sys.stderr.write('\nError: use of hybrid Illumina/Nanopore assembly requires -r1, -r2, and -l to be specified.\n')
			raise ValueError

	if arg_obj.metagenomic:
		if arg_obj.forward and arg_obj.reverse and arg_obj.long_reads:
			sys.stderr.write('\nError: metagenomic assembly cannot be used with hybrid assembly.  Please use paired-end '
			                 'Illumina or exclusively long-read Nanopore assembly options instead.\n')
			raise ValueError

		elif not arg_obj.reverse:
			sys.stderr.write('\nError: metagenomic assembly with Illumina data requires -r1 and -r2 to be specified (paired-end)')
			raise ValueError

	if not arg_obj.reverse and not os.path.isdir(arg_obj.forward):
		sys.stderr.write('\nError: Using only -r1 requires the -r1 path to be a Nanopore fastq_pass directory.  Single-end '
		                 'assembly for Illumina reads is not supported.')
		raise ValueError

	if arg_obj.reverse and not arg_obj.long_reads:
		if os.path.isdir(arg_obj.forward):
			sys.stderr.write('\nUse of forward and reverse reads requires -r1 and -r2 to be Illumina FASTQ reads (paired-end).\n')
			raise ValueError

	if arg_obj.reverse and arg_obj.long_reads:
		if not os.path.isdir(arg_obj.long_reads):
			sys.stderr.write('\nUse of hybrid assembly requires -l to be a Nanopore fastq_pass directory and for '
			                 '-r1 and -r2 to be Illumina paired-end read FASTQ files.\n')
			raise ValueError


def get_real_dir_from_file(infile):
	return os.path.dirname(os.path.realpath(infile))


def get_real_dir_from_dir(indir):
	return os.path.realpath(indir)


if __name__ == '__main__':
	args = parser.parse_args()
	validate_input_args(args)
	temp_dir = get_real_dir_from_dir(args.work_dir)
	nextflow_path = '/'.join(get_real_dir_from_file(sys.argv[0]).split('/')[:-1]) + '/denovo.nf'

	if args.slurm:
		nextflow_config = nextflow_path.replace('denovo.nf', SLURM_CONFIG)
	else:
		nextflow_config = nextflow_path.replace('denovo.nf', 'nextflow.config')

	nextflow_arglist = [
		nextflow_path,
		'-config',
		nextflow_config,
		'-w',
		temp_dir,
		'-resume',
		'--threads',
		args.threads
	]

	if args.host_reference:
		nextflow_arglist += ['--reference', args.host_reference]
	else:
		nextflow_arglist += ['--reference', 'NONE_REF']

	if args.singularity:
		nextflow_arglist += ['-with-singularity', args.singularity]
	else:
		nextflow_arglist += ['-with-singularity', SINGULARITY_LIB_PATH]

	sys.stdout.write('\nBeginning SLSS-DENOVO Pipeline at {:%Y-%m-%d %H:%M:%S}\n\n'.format(
		datetime.now()
	))

	if not args.reverse and os.path.isdir(args.forward):
		# Long read assembly using Flye
		if args.force:
			concatenate_nanopore_reads(args.forward, temp_dir)
		elif not os.path.exists(temp_dir + '/' + NANOPORE_TEMP_FILENAME):
			concatenate_nanopore_reads(args.forward, temp_dir)

		nextflow_arglist += [
			'--forward',
			temp_dir + '/' + NANOPORE_TEMP_FILENAME,
			'--nanopore'
		]

		if args.metagenomic:
			sys.stdout.write('Performing Nanopore long read metagenomic de novo assembly using Flye...\n\n')
			nextflow_arglist += ['--metagenomic']
		else:
			sys.stdout.write('Performing Nanopore long read de novo assembly using Flye...\n\n')
		p = subprocess.Popen([str(x) for x in nextflow_arglist])

	elif args.reverse and not args.long_reads:
		# Paired-end Illumina with SPAdes

		nextflow_arglist += [
			'--forward',
			args.forward,
			'--reverse',
			args.reverse
		]

		if args.metagenomic:
			sys.stdout.write('Performing paired-end Illumina de novo assembly using metaSPAdes...\n\n')
			nextflow_arglist += ['--metagenomic']
		else:
			sys.stdout.write('Performing paired-end Illumina de novo assembly using SPAdes...\n\n')
		p = subprocess.Popen([str(x) for x in nextflow_arglist])

	elif args.reverse and args.long_reads:
		# Hybrid assembly with SPAdes

		nextflow_arglist += [
			'--forward',
			args.forward,
			'--reverse',
			args.reverse,
			'--hybrid',
			temp_dir + '/' + NANOPORE_TEMP_FILENAME
		]

		if args.force:
			concatenate_nanopore_reads(args.long_reads, temp_dir)
		elif not os.path.exists(temp_dir + '/' + NANOPORE_TEMP_FILENAME):
			concatenate_nanopore_reads(args.long_reads, temp_dir)
		sys.stdout.write('Performing hybrid Illumina paired-end with Nanopore long read de novo assembly using SPAdes...\n\n')
		p = subprocess.Popen([str(x) for x in nextflow_arglist])
	else:
		sys.stderr.write('\nThe combination of arguments passed to -r1, -r2, and -l are not a valid pipeline.'
		                 ' See -h for usage.\n')
		raise ValueError
	exit_code = p.wait()
	sys.stdout.write('\nSLSS-DENOVO exited with code: {} at {:%Y-%m-%d %H:%M:%S}\n\n'.format(
		exit_code,
		datetime.now()
	))
