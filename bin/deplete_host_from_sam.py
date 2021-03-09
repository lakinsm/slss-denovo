#!/usr/bin/env python3

import argparse
import sys
import os


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files. Outputs non-aligned reads into one or
	two FASTQ files depending on input arguments (single-end or long read vs. paired-end), resulting in a depletion of
	aligned reads (host).  Paired reads where one mate maps but the other does not are both output.
	"""

	def __init__(self, sam_path, forward_outpath, reverse_outpath=None, filter_ns=10):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param forward_outpath: STR, path to forward read output file in FASTQ format
		:param reverse_outpath: STR, path to reverse read output file in FASTQ format (optional if single-end or long)
		:param filter_ns: FLOAT, filter out reads containing more than this % Ns in the nucleotide string
		"""
		self.total_reads_processed = 0
		self.reads_pass_filter = 0
		self.n_filter_threshold = float(filter_ns)
		self.pass_filter = False
		self.prev_header = ''
		self.prev_qual = ''
		self.prev_seq = ''
		self.prev_rev = False
		self.prev_pass_filter = False
		self.sam_path = sam_path
		self.forward_fp = forward_outpath
		self.forward_handle = None
		self.paired = False
		if reverse_outpath:
			self.reverse_fp = reverse_outpath
			self.paired = True
		else:
			self.reverse_fp = None
		self.reverse_handle = None

		# Skip SAM header information
		self._open()
		self.line = self.handle.readline()
		if self.line:
			while self.line[0] == '@':
				self.line = self.handle.readline()

	def __iter__(self):
		return self

	def __next__(self):
		"""
		Iterator for SAM file handle, outputs non-aligned read entries.
		"""
		if self.paired:
			if not self.line:
				self._close()
				raise StopIteration
			entries = self.line.split('\t')
			sam_flag = int(entries[1])
			# Skip alternative alignments
			if sam_flag & 2048 == 0:
				self.total_reads_processed += 1
				if sam_flag & 4 != 0:
					query_header = entries[0]
					query_seq = entries[9]
					query_qual = entries[10]
					rev = (sam_flag & 128) != 0
					if (100*float(query_seq.count('N')) / float(len(query_seq))) <= self.n_filter_threshold:
						self.pass_filter = True
						self.reads_pass_filter += 1
						if query_header == self.prev_header and (self.prev_pass_filter or self.pass_filter):
							if rev and not self.prev_rev:
								self.reverse_handle.write('@{}\n{}\n+\n{}\n'.format(
									query_header,
									query_seq,
									query_qual
								))
								self.forward_handle.write('@{}\n{}\n+\n{}\n'.format(
									self.prev_header,
									self.prev_seq,
									self.prev_qual
								))
							elif not rev and self.prev_rev:
								self.reverse_handle.write('@{}\n{}\n+\n{}\n'.format(
									self.prev_header,
									self.prev_seq,
									self.prev_qual
								))
								self.forward_handle.write('@{}\n{}\n+\n{}\n'.format(
									query_header,
									query_seq,
									query_qual
								))
							else:
								sys.stderr.write('Reads with the same name have matching orientations: {}'.format(
									query_header
								))
								raise ValueError
				else:
					self.pass_filter = False
				self.prev_header = entries[0]
				self.prev_seq = entries[9]
				self.prev_qual = entries[10]
				self.prev_rev = (sam_flag & 128) != 0
				self.prev_pass_filter = self.pass_filter
		else:
			if not self.line:
				self._close()
				raise StopIteration
			entries = self.line.split('\t')
			sam_flag = int(entries[1])
			# Skip alternative alignments
			if sam_flag & 2048 == 0:
				self.total_reads_processed += 1
				if sam_flag & 4 != 0:
					query_header = entries[0]
					query_seq = entries[9]
					query_qual = entries[10]
					if (100*float(query_seq.count('N')) / float(len(query_seq))) <= self.n_filter_threshold:
						self.reads_pass_filter += 1
						self.forward_handle.write('@{}\n{}\n+\n{}\n'.format(
							query_header,
							query_seq,
							query_qual
						))
		self.line = self.handle.readline()

	def _open(self):
		self.handle = open(self.sam_path, 'r')
		self.forward_handle = open(self.forward_fp, 'w')
		if self.reverse_fp:
			self.reverse_handle = open(self.reverse_fp, 'w')

	def _close(self):
		self.handle.close()
		self.forward_handle.close()
		if self.reverse_fp:
			self.reverse_handle.close()


def get_real_file_from_file(infile):
	return os.path.realpath(infile)


parser = argparse.ArgumentParser('deplete_host_from_sam.py')
parser.add_argument('sam_file', type=str, help='Input SAM filepath')
parser.add_argument('-r1', '--forward', type=str, default=None, required=True,
                    help='Output forward single- or paired-end read FASTQ file or Nanopore long read FASTQ file')
parser.add_argument('-r2', '--reverse', type=str, default=None, required=False,
                    help='Output reverse paired-end read FASTQ file')


if __name__ == '__main__':
	args = parser.parse_args()
	sam_parser = SamParser(args.sam_file, args.forward, args.reverse)
	for _ in sam_parser:
		continue
	forward_path = get_real_file_from_file(args.forward)
	fname = forward_path.split('/')[-1].split('.')[0]
	logpath = '/'.join(forward_path.split('/')[:-1]) + '/' + fname + '_depletion_stats.log'
	with open(logpath, 'w') as log:
		log.write('total_reads: {}\nreads_pass_filter: {}\npercent_remaining_after_depletion: {} %\n'.format(
			sam_parser.total_reads_processed,
			sam_parser.reads_pass_filter,
			100 * float(sam_parser.reads_pass_filter) / float(sam_parser.total_reads_processed)
		))
