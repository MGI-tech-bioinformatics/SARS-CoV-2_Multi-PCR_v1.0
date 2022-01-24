#!/usr/bin/env python3

import sys
import os
from optparse import OptionParser

def stat_qc_iden(sample_l,resultdir,outdir):
	qc_sum = open('%s/Summary.QC.xls'%(outdir),'w')
	iden_sum = open('%s/Summary.Identification.xls'%(outdir),'w')
	n = 0
	for sample in sample_l:
		qc_file = '%s/%s/05.Stat/QC.txt'%(resultdir,sample)
		iden_file = '%s/%s/05.Stat/Identification.txt'%(resultdir,sample)
		fq = open(qc_file,'r')
		fi = open(iden_file,'r')
		qc_head = fq.readline()
		iden_head = fi.readline()
		if n == 0:
			qc_sum.write(qc_head+fq.readline())
			iden_sum.write(iden_head+fi.readline())
		else:
			qc_sum.write(fq.readline())
			iden_sum.write(fi.readline())
		n += 1
	qc_sum.close()
	iden_sum.close()
	return

def stat_variant(sample_l,resultdir,outdir):
	var_sum = open('%s/Summary.Mutation.xls'%(outdir),'w')
	n = 0
	for sample in sample_l:
		var_file = '%s/%s/05.Stat/%s.snpEff.anno.txt'%(resultdir,sample,sample)
		with open(var_file,'r') as fv:
			head = fv.readline()
			for line in fv:
				if n == 0:
					var_sum.write(head+line)
				else:
					var_sum.write(line)
		n += 1
	var_sum.close()
	return

def stat_fa(resultdir,outdir):
	os.system('cat %s/*/05.Stat/*.Consensus.fa > %s/Summary.Consensus.fa'%(resultdir,outdir))
	return

def main():
	if len(sys.argv) < 2:
		print('ERROR: Invalid input, use -h/--help to output help information.')
		sys.exit(1)
	parser = OptionParser()
	parser.add_option('-l', dest = 'opt_l', help = 'Sample list', type = 'string')
	parser.add_option('-r', dest = 'opt_r', help = 'Result directory', type = 'string')
	parser.add_option('-o', dest = 'opt_o', help = 'Output directory', type = 'string')
	optlist, args = parser.parse_args()

	sample_list = optlist.opt_l
	resultdir = optlist.opt_r
	outdir = optlist.opt_o

	sample_l = []
	with open(sample_list,'r') as fs:
		for line in fs:
			sample = line.split()[0]
			if sample not in sample_l:
				sample_l.append(sample)

	stat_qc_iden(sample_l,resultdir,outdir)
	stat_variant(sample_l,resultdir,outdir)
	stat_fa(resultdir,outdir)
	return

if __name__ == '__main__':
	main()
