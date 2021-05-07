#!/usr/bin/env python3

import sys
import re
import os
import getopt
import gzip
import pandas as pd
from statistics import mean

def usage():
	"""
Statistics of FQ and mapping information.
	-d	result directory
	-l	sample list
	-t	FASTQ file type [PE/SE]
	"""
	print(usage.__doc__)
	return

if len(sys.argv) < 2:
	usage()
	sys.exit(1)

try:
	opts,args = getopt.getopt(sys.argv[1:],"hd:l:w:t:")
except:
	sys.exit(1)

for key, value in opts:
	if key == '-h':
		usage()
		sys.exit()
	if key == '-d':
		resultdir = value
	if key == '-l':
		sample_list = value
	if key == '-t':
		fqtype = value

def stat_fq(fqstat,outdict,fqtype):
	with open(fqstat,'r') as fs:
		for line in fs:
			info = line.split()
			if 'Total number of reads' in line:
				if fqtype == 'PE':
					raw_reads = str(int(info[4])*2)
					clean_reads = str(int(info[6])*2)
					outdict['raw_reads'] = raw_reads
					outdict['clean_reads'] = clean_reads
				elif fqtype == 'SE':
					raw_reads = info[4]
					clean_reads = info[6]
					outdict['raw_reads'] = raw_reads
					outdict['clean_reads'] = clean_reads
			if 'Total number of bases' in line:
				total_bases = int(info[4])*2
			if 'Number of base C (%)' in line:
				if fqtype == 'PE':
					C1 = info[5]
					C2 = line.split('(')[3].split(')')[1].strip()
					C_number = int(C1) + int(C2)
				elif fqtype == 'SE':
					C_number = int(info[5])
			if 'Number of base G (%)' in line:
				if fqtype == 'PE':
					G1 = info[5]
					G2 = line.split('(')[3].split(')')[1].strip()
					G_number = int(C1) + int(C2)
				elif fqtype == 'SE':
					G_number = int(info[5])
			if 'Number of filtered reads' in line:
				filter_rate = line.split('(')[2].split(')')[0].strip()
				outdict['filter_rate'] = filter_rate
				outdict['clean_rate'] = str('%.2f'%(100 - float(filter_rate.strip('%'))))+'%'
			if 'Number of base calls with quality value of 20 or higher (Q20+)' in line:
				if fqtype == 'PE':
					raw_Q20_fq1 = line.split('(')[3].split(')')[0].strip('%')
					clean_Q20_fq1 = line.split('(')[4].split(')')[0].strip('%')
					raw_Q20_fq2 = line.split('(')[5].split(')')[0].strip('%')
					clean_Q20_fq2 = line.split('(')[6].split(')')[0].strip('%')
					raw_Q20 = '%.2f'%((float(raw_Q20_fq1) + float(raw_Q20_fq2))/2)
					clean_Q20 = '%.2f'%((float(clean_Q20_fq1) + float(clean_Q20_fq2))/2)
					outdict['raw_Q20'] = raw_Q20 + '%'
					outdict['clean_Q20'] = clean_Q20 + '%'
				elif fqtype == 'SE':
					raw_Q20 = line.split('(')[3].split(')')[0].strip('%')
					clean_Q20 = line.split('(')[4].split(')')[0].strip('%')
					outdict['raw_Q20'] = raw_Q20
					outdict['clean_Q20'] = clean_Q20
			if 'Number of base calls with quality value of 30 or higher (Q30+)' in line:
				if fqtype == 'PE':
					raw_Q30_fq1 = line.split('(')[3].split(')')[0].strip('%')
					clean_Q30_fq1 = line.split('(')[4].split(')')[0].strip('%')
					raw_Q30_fq2 = line.split('(')[5].split(')')[0].strip('%')
					clean_Q30_fq2 = line.split('(')[6].split(')')[0].strip('%')
					raw_Q30 = '%.2f'%((float(raw_Q30_fq1) + float(raw_Q30_fq2))/2)
					clean_Q30 = '%.2f'%((float(clean_Q30_fq1) + float(clean_Q30_fq2))/2)
					outdict['raw_Q30'] = raw_Q30 + '%'
					outdict['clean_Q30'] = clean_Q30 + '%'
				elif fqtype == 'SE':
					raw_Q30 = line.split('(')[3].split(')')[0].strip('%')
					clean_Q30 = line.split('(')[4].split(')')[0].strip('%')
					outdict['raw_Q30'] = raw_Q30
					outdict['clean_Q30'] = clean_Q30
	try:
		outdict['GC_content'] = str('%.2f'%(100*(C_number + G_number)/total_bases))+'%'
	except:
		outdict['GC_content'] = '0%'
	return

#def stat_map(mapcov,bamstat,uniqbamstat,primer_reads,outdict):
def stat_map(mapcov,outdict):
	with open(mapcov,'r') as fm:
		for line in fm:
			info = line.split('\t')
			if '[Total] Fraction of Mapped Data(Mb)' in line:
				Mapping_rate = info[1].strip()
				outdict['Mapping_rate'] = Mapping_rate
			if '[Target] Target Reads' in line:
				outdict['Target_reads'] = info[1].strip()
			if '[Target] Fraction of Target Reads in all reads' in line:
				capture_rate = info[1].strip()
				outdict['capture_rate'] = capture_rate
			if '[Target] Coverage (>0x)' in line:
				outdict['Coverage_1X'] = info[1].strip()
			if '[Target] Coverage (>=100x)' in line:
				Coverage_100X = info[1].strip()
				outdict['Coverage_100X'] = Coverage_100X
	return

def stat_depth(depth_tsv,outdict):
	n_total = 1
	n_1X = 1
	n_100X = 0
	n_depth = 0
	with gzip.open(depth_tsv,'rt') as fd:
		for line in fd:
			if '#' in line:
				continue
			depth = line.split()[2]
			n_total += 1
			if int(depth) == 0:
				continue
			if int(depth) >= 100:
				n_100X += 1
			n_1X += 1
			n_depth += int(depth)
	average_depth_total = '%.2f'%(n_depth/n_total)
	average_depth_1X = '%.2f'%(n_depth/n_1X)
	outdict['Genome_Average_Depth'] = average_depth_total
	outdict['1X_Region_Average_Depth'] = average_depth_1X
	outdict['1X_size'] = n_1X
	outdict['100X_size'] = n_100X
	return

def stat_sample(resultdir,sample,out_dict,lambda_dict,GAPDH_dict):
	file_fq_stat = resultdir + '/' + sample + '/01.Clean/Basic_Statistics_of_Sequencing_Quality.txt'
	file_mapcov = resultdir + '/' + sample + '/03.covdep/coverage.report'
	file_depth_tsv = resultdir + '/' + sample + '/03.covdep/depth.tsv.gz'
	file_lambdacov = resultdir + '/' + sample + '/03.covdep/lambda_cov/coverage.report'
	file_GAPDHcov = resultdir + '/' + sample + '/03.covdep/GAPDH_cov/coverage.report'
	stat_fq(file_fq_stat,out_dict,fqtype)
	stat_map(file_mapcov,out_dict)
	stat_map(file_lambdacov,lambda_dict)
	stat_map(file_GAPDHcov,GAPDH_dict)
	stat_depth(file_depth_tsv,out_dict)
	return

def main():
	with open(sample_list,'r') as fs:
		for line in fs:
			out_dict = {}
			lambda_dict = {}
			GAPDH_dict = {}
			sample = line.split()[0]
			if not os.path.exists('%s/%s/05.Stat'%(resultdir,sample)):
				os.makedirs('%s/%s/05.Stat'%(resultdir,sample))
			file_QC_stat = open('%s/%s/05.Stat/QC.txt'%(resultdir,sample),'w')
			file_iden_out = open('%s/%s/05.Stat/Identification.txt'%(resultdir,sample),'w',encoding='utf-8')
			stat_sample(resultdir,sample,out_dict,lambda_dict,GAPDH_dict)
			lambda_reads = lambda_dict['Target_reads']
			GAPDH_reads = GAPDH_dict['Target_reads']
			try:
				lambda_rate = str('%.2f'%(100*int(lambda_reads)/int(out_dict['clean_reads'])))+'%'
			except:
				lambda_rate = '0%'
			try:
				GAPDH_rate = str('%.2f'%(100*int(GAPDH_reads)/int(out_dict['clean_reads'])))+'%'
			except:
				GAPDH_rate = '0%'
			try:
				PCT = str('%.2f'%(100*int(out_dict['Target_reads'])/(int(lambda_reads)+int(out_dict['Target_reads']))))+'%'
			except:
				PCT = '0%'
			Target_reads = out_dict['Target_reads']
			Coverage_1X = float(out_dict['Coverage_1X'].strip('%'))
			iden_number = float(PCT.strip('%'))
			if iden_number >= 0.1 and Coverage_1X >= 1:
				iden = 'Positive'
			elif iden_number < 0.05:
				iden = 'Negative'
			else:
				iden = 'Indetermination'
			#if Coverage_1X < 30:
			#	iden = 'Negative'
			file_QC_stat.write('Sample\tRaw_Q30\tGC_Content\tRaw_Reads\tClean_Reads\tClean_Rate\tMapping_Rate\n')
			file_QC_stat.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(sample,out_dict['raw_Q30'],out_dict['GC_content'],out_dict['raw_reads'],out_dict['clean_reads'],out_dict['clean_rate'],out_dict['Mapping_rate']))
			file_iden_out.write('Sample\tClean_Reads\tGAPDH_Reads\tLambda_Reads\tSARS-CoV-2_Reads\tGAPDH_Rate\tLambda_Rate\tSARS-CoV-2_Reads_Pct\tGenome_Average_Depth\t≥1X_Region_Average_Depth\t≥1X_Cov\t≥1X_size\t≥100X_Cov\t≥100X_size\tIdentification_Result\n')
			file_iden_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n'%(sample,out_dict['clean_reads'],GAPDH_reads,lambda_reads,Target_reads,GAPDH_rate,lambda_rate,PCT,out_dict['Genome_Average_Depth'],out_dict['1X_Region_Average_Depth'],out_dict['Coverage_1X'],out_dict['1X_size'],out_dict['Coverage_100X'],out_dict['100X_size'],iden))
			file_QC_stat.close()
			file_iden_out.close()
			df_QC = pd.read_csv('%s/%s/05.Stat/QC.txt'%(resultdir,sample), sep='\t')
			df_iden = pd.read_csv('%s/%s/05.Stat/Identification.txt'%(resultdir,sample), sep='\t')
			df_QC.to_excel('%s/%s/05.Stat/QC.xlsx'%(resultdir,sample), 'Sheet1', index=False)
			df_iden.to_excel('%s/%s/05.Stat/Identification.xlsx'%(resultdir,sample), 'Sheet1', index=False)
	return

if __name__ == '__main__':
	main()
