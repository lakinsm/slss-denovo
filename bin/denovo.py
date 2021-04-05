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
parser.add_argument('-r1', '--forward', default=None, required=False, type=str,
                    help='Forward paired-end Illumina FASTQ')
parser.add_argument('-r2', '--reverse', default=None, required=False, type=str,
                    help='Reverse paired-end Illumina FASTQ')
parser.add_argument('-l', '--long_reads', default=None, required=False, type=str,
                    help='Nanopore fastq_pass directory, default Flye assembly, for hybrid must have -r1, -r2, and -l specified')
parser.add_argument('-o', '--output', default=None, required=True, type=str,
                    help='Path to output directory')
parser.add_argument('--host_reference', default=None, required=False, type=str,
                    help='Host reference sequence FASTA for host depletion (optional)')
parser.add_argument('--target_reference', default=None, required=False, type=str,
                    help='Target reference sequence FASTA for assembly metric calculation (optional)')
parser.add_argument('-d', '--diamond_db', default=None, required=True, type=str,
                    help='Location of Diamond indexed BLAST nr database file')
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
parser.add_argument('--flye_overrides', required=False, type=str,
                    help='Custom arguments for Flye, overrides -m')


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

	if arg_obj.reverse and not arg_obj.long_reads:
		if os.path.isdir(arg_obj.forward):
			sys.stderr.write('\nUse of forward and reverse reads requires -r1 and -r2 to be Illumina FASTQ reads (paired-end).\n')
			raise ValueError

	if arg_obj.long_reads:
		if not os.path.isdir(arg_obj.long_reads):
			sys.stderr.write('\nUse of Nanopore data requires -l to be a Nanopore fastq_pass directory.\n')
			raise ValueError


def get_real_dir_from_file(infile):
	return os.path.dirname(os.path.realpath(infile))


def get_real_dir_from_dir(indir):
	return os.path.realpath(indir)


def bwa_index_not_present_message(host):
	sys.stdout.write(
		'\nWARNING: No BWA index detected for host reference genome: {}.  If this genome is large, '
		'it may take up to several hours to index with BWA.  It is recommended to provide a '
		'pre-indexed host genome to prevent wasted computing time rebuilding the index.\n'.format(
			host
		))


if __name__ == '__main__':
	args = parser.parse_args()
	validate_input_args(args)
	temp_dir = get_real_dir_from_dir(args.work_dir)
	nextflow_path = '/'.join(get_real_dir_from_file(sys.argv[0]).split('/')[:-1]) + '/denovo.nf'
	diamond_path = os.path.realpath(args.diamond_db)

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
		args.threads,
		'--output',
		args.output,
		'--scratch_dir',
		temp_dir,
		'--diamond_db',
		diamond_path
	]

	if args.host_reference:
		nextflow_arglist += ['--host', args.host_reference]
		host_files = glob.glob(args.host_reference + '.*')
		print(host_files)
		if len(host_files) > 1:
			use_existing = True
			for ext in ('.amb', '.ann', '.bwt', '.pac', '.sa'):
				local_exists = False
				for f in host_files:
					local_exists |= ext in f
				use_existing &= local_exists
			if use_existing:
				nextflow_arglist += ['--use_index']
				sys.stdout.write('\nExisting BWA indices found for host genome: {}'.format(
					args.host_reference
				))
			else:
				bwa_index_not_present_message(args.host_reference)
		else:
			bwa_index_not_present_message(args.host_reference)
	else:
		nextflow_arglist += ['--host', 'NONE_REF']

	if args.target_reference:
		nextflow_arglist += ['--reference', args.target_reference]
	else:
		nextflow_arglist += ['--reference', 'NONE_REF']

	if args.singularity:
		nextflow_arglist += ['-with-singularity', args.singularity]
	else:
		nextflow_arglist += ['-with-singularity', SINGULARITY_LIB_PATH]

	sys.stdout.write('\nBeginning SLSS-DENOVO Pipeline at {:%Y-%m-%d %H:%M:%S}\n\n'.format(
		datetime.now()
	))

	if not args.forward and os.path.isdir(args.long_reads):
		# Long read assembly using Flye
		if args.force:
			concatenate_nanopore_reads(args.long_reads, temp_dir)
		elif not os.path.exists(temp_dir + '/' + NANOPORE_TEMP_FILENAME):
			concatenate_nanopore_reads(args.long_reads, temp_dir)

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

	elif args.forward and args.reverse and not args.long_reads:
		# Paired-end Illumina with SPAdes

		nextflow_arglist += [
			'--forward',
			args.forward,
			'--reverse',
			args.reverse
		]

		if args.metagenomic:
			sys.stdout.write('Performing paired-end Illumina de novo assembly using metaviralSPAdes...\n\n')
			nextflow_arglist += ['--metagenomic']
		else:
			sys.stdout.write('Performing paired-end Illumina de novo assembly using SPAdes...\n\n')
		p = subprocess.Popen([str(x) for x in nextflow_arglist])

	elif args.forward and args.reverse and args.long_reads:
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
		if args.metagenomic:
			sys.stdout.write('Performing hybrid Illumina/Nanopore de novo assembly using metaviralSPAdes...\n\n')
			nextflow_arglist += ['--metagenomic']
		else:
			sys.stdout.write('Performing hybrid Illumina/Nanopore de novo assembly using SPAdes...\n\n')
		p = subprocess.Popen([str(x) for x in nextflow_arglist])
	else:
		sys.stderr.write('\nThe combination of arguments passed to -r1, -r2, and -l are not a valid pipeline.'
		                 ' Use either -l alone for Flye assembly, -r1 and -r2 together for paired-end assembly, '
		                 'or -r1 -r2 and -l for hybrid assembly.'
		                 ' See -h for usage.\n')
		raise ValueError
	exit_code = p.wait()
	sys.stdout.write('\nSLSS-DENOVO exited with code: {} at {:%Y-%m-%d %H:%M:%S}\n\n'.format(
		exit_code,
		datetime.now()
	))
