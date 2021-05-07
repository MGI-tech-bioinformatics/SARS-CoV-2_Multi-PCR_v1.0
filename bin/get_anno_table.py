#!/usr/bin/env python3

import sys
import re
import pandas as pd

file_in = sys.argv[1]
resultdir = sys.argv[2]
sample = sys.argv[3]

anno_list = [1,3,5,7,9,10,2]
#DIS 14

with open(resultdir+'/Identification.txt','r') as fi:
	head = fi.readline()
	line = fi.readline()
	iden = line.split()[-1]

file_out = open('%s/%s.snpEff.anno.txt'%(resultdir,sample),'w')

with open(file_in,'r') as fa:
	for line in fa:
		if '##INFO=<ID=ANN' in line:
			header = line.split("'")[1].split('|')
			header_out = 'Sample\tReference\tPosition\tRefBase\tAltBase'
			for index in anno_list:
				if 'HGVS.c' in header[index]:
					header_out += '\tGene_Position_Change'
				elif 'HGVS.p' in header[index]:
					header_out += '\tAA_Change'
				else:
					header_out += '\t'+header[index]
			header_out += '\tReference_Depth\tAlternate_Depth\tAllele_Frequency\n'
			#print(header_out)
			file_out.write(header_out)
			continue
		if '#' in line:
			continue
		i = line.split()
		pos,ref,alt,INFO,vcf_tag,vcf_form = i[1],i[3],i[4],i[7],i[-2],i[-1]
		ANN = re.match(r'.+ANN=(.+)',INFO).group(1).split(',')[0].split('|')
		anno_out = '%s\t%s\t%s\t%s\t%s'%(sample,'2019-nCoV',pos,ref,alt)
		for index in anno_list:
			if not ANN[index]:
				anno_out += '\t.'
			else:
				anno_out += '\t'+ANN[index]
		if ANN[14]:
			anno_out += ';DISTANCE='+ANN[14]

		vcf_dict = {}
		tag_list = vcf_tag.split(':')
		form_list = vcf_form.split(':')
		for n in range(0,len(tag_list)):
			vcf_dict[tag_list[n]] = form_list[n]
		RO = vcf_dict['RO']
		AO = vcf_dict['AO']
		Freq = str('%.2f'%(100*int(AO)/(int(RO)+int(AO))))+'%'
		anno_out += '\t%s\t%s\t%s\n'%(RO,AO,Freq)
		#print(anno_out)
		if iden == 'Positive':
			file_out.write(anno_out)

file_out.close()

df_anno = pd.read_csv('%s/%s.snpEff.anno.txt'%(resultdir,sample), sep='\t')
df_anno.to_excel('%s/%s.snpEff.anno.xlsx'%(resultdir,sample), 'Sheet1', index=False)
