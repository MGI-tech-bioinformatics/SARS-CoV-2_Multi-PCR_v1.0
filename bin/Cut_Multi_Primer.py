#!/usr/bin/env python3

import sys
import os
import pysam
import time
from subprocess import check_call
from optparse import OptionParser

def get_primer_info(file_primer):
	R_set = set()
	with open(file_primer,'r') as fa:
		for line in fa:
			chrn, primer1, start1, end1, primer2, start2, end2 = line.split()
			if int(start1) > int(start2):
				R1_start = start2
				R1_end = end2
				R2_start = start1
				R2_end = end1
			elif int(start1) < int(start2):
				R1_start = start1
				R1_end = end1
				R2_start = start2
				R2_end = end2
			R_set.add('%s\t%s\t%s\t%s\t%s'%(chrn,R1_start,R1_end,R2_start,R2_end))
	return R_set

def cut_primer(chrn,pos,pos_end,reverse,seq,qua,R_set):
	result = 'NA'
	for i in R_set:
		p_chrn, p_start1, p_end1, p_start2, p_end2 = i.split()
		if chrn == p_chrn and int(pos) - 10 <= int(p_start1) <= int(pos) + 10 and reverse == 'mate_reverse':
			result = 'PASS'
			cut_seq = seq
			cut_qua = '!'*(int(p_end1) - int(p_start1) + 1) + qua[(int(p_end1) - int(p_start1) + 1):]
			break
		elif chrn == p_chrn and int(pos_end) - 10 <= int(p_end2) <= int(pos_end) + 10 and reverse == 'read_reverse':
			result = 'PASS'
			cut_seq = seq
			cut_qua = qua[0:(len(seq) - int(p_end2) + int(p_start2) - 1)] + '!'*(int(p_end2) - int(p_start2) + 1)
			break
	if result == 'NA':
		cut_seq = seq
		cut_qua = qua
	return cut_seq, cut_qua, result

def revcom(seq):
	return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

def gzip_fq(fastq):
	check_call('gzip %s'%(fastq),shell=True)
	return

def usage():
	"""
Primer adapter processing program.
----------------------------------
Version: 1.1.0

Usage:
	Cut_Multi_Primer.py -p <STR> -b <STR> -s <STR> -o <STR>

Param list:
	-h		Help information
	-p		Primer position file.
	-b		Input bam file
	-s		Sample name
	-o		Output directory
	"""
	print(usage.__doc__)
	sys.exit(1)
	return

if __name__ == '__main__':
	if len(sys.argv) < 2:
		usage()
	parser = OptionParser()
	parser.add_option('-p', dest = 'opt_p', help = 'Primer position file', type = 'string')
	parser.add_option('-b', dest = 'opt_b', help = 'Input bam file', type = 'string')
	parser.add_option('-s', dest = 'opt_s', help = 'Sample name', type = 'string')
	parser.add_option('-o', dest = 'opt_o', help = 'Output directory', type = 'string')
	optlist, args = parser.parse_args()

	file_primer = optlist.opt_p
	inbam = optlist.opt_b
	sample = optlist.opt_s
	outdir = optlist.opt_o
	outbam = '%s/%s.bam'%(outdir,sample)

	print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Started processing ')

	R_set = get_primer_info(file_primer)
	file_inbam = pysam.AlignmentFile(inbam,'rb',check_sq=False)
	header = file_inbam.header
	file_outbam = pysam.AlignmentFile(outbam,'wb',header=header,check_sq=False)
	log = '%s/%s.cutprimer.log'%(outdir,sample)

	fq1_dict = {}
	fq2_dict = {}
	#abnormal_set = set()

	n = 0
	a = 0

	for r in file_inbam:
		readid = r.query_name
		seq = r.seq
		qua = r.qual
		try:
			chrn = r.reference_name
			pos = int(r.pos) + 1
			pos_end = pos
			cigar = r.cigartuples
			for key, value in cigar:
				if key == 0 or key == 1 or key == 2:
					pos_end += value
		except:
			continue

		if r.is_unmapped:
			continue
		if r.is_reverse:
			reverse = 'read_reverse'
		else:
			reverse = 'mate_reverse'

		if r.is_secondary:
			continue
		#if chrn != 'nCoV' and n > 1:
		if chrn != 'MN908947.3' and n > 1:
			continue
		else:
			n += 1
		cut_seq, cut_qua, result = cut_primer(chrn,pos,pos_end,reverse,seq,qua,R_set)
		if result == 'NA':
			a += 1
			continue
			#abnormal_set.add(readid)
		
		r.seq = cut_seq
		r.qual = cut_qua
		file_outbam.write(r)
	file_inbam.close()
	file_outbam.close()
	print('[',time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),'] Finished processing ')
