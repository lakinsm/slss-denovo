#!/usr/bin/env python3

import argparse


class SamParser(object):
	"""
	Line-by-line parsing of Sequence Alignment Map (SAM) UTF-8 encoded files. Outputs non-aligned reads into one or
	two FASTQ files depending on input arguments (single-end or long read vs. paired-end), resulting in a depletion of
	aligned reads (host).
	"""

	def __init__(self, sam_path, forward_outpath, reverse_outpath=None):
		"""
		Initialize object and data structures for storing coverage and coverage over time, if specified.
		:param sam_path: STR, path to input SAM file
		:param forward_outpath: STR, path to forward read output file in FASTQ format
		:param reverse_outpath: STR, path to reverse read output file in FASTQ format (optional if single-end or long)
		"""
		self.sam_path = sam_path
		self.forward_fp = forward_outpath
		self.forward_handle = None
		if reverse_outpath:
			self.reverse_fp = reverse_outpath
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
		if not self.line:
			self._close()
			raise StopIteration
		entries = self.line.split('\t')
		sam_flag = int(entries[1])
		# Skip aligned reads
		while sam_flag & 4 == 0:
			self.line = self.handle.readline()
			if not self.line:
				self._close()
				raise StopIteration
			entries = self.line.split('\t')
			sam_flag = int(entries[1])
		query_header = entries[0]
		query_seq = entries[9]
		query_qual = entries[10]
		rev = (sam_flag & 128) != 0
		if rev:
			self.reverse_handle.write('@{}\n{}\n+\n{}\n'.format(
				query_header,
				query_seq,
				query_qual
			))
		else:
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


parser = argparse.ArgumentParser('deplete_host_from_sam.py')
parser.add_argument('sam_file', type=str, default=None, required=True,
                    help='Input SAM filepath')
parser.add_argument('-r1', '--forward', type=str, default=None, required=True,
                    help='Output forward single- or paired-end read FASTQ file or Nanopore long read FASTQ file')
parser.add_argument('-r2', '--reverse', type=str, default=None, required=False,
                    help='Output reverse paired-end read FASTQ file')


if __name__ == '__main__':
	args = parser.parse_args()
	sam_parser = SamParser(args.sam_file, args.forward, args.reverse)
	for _ in sam_parser:
		continue
