#!/usr/bin/env python3

import sys
import re


file_in = sys.argv[1]

anno_list = [1,3,5,7,9,10,2]
#DIS 14

with open(file_in,'r') as fa:
	for line in fa:
		if '##INFO=<ID=ANN' in line:
			header = line.split("'")[1].split('|')
			header_out = 'Reference\tPosition\tRefBase\tAltBase'
			for index in anno_list:
				if 'HGVS.c' in header[index]:
					header_out += '\tGene_Position_Change'
				elif 'HGVS.p' in header[index]:
					header_out += '\tAA_Change'
				else:
					header_out += '\t'+header[index]
			print(header_out)
			continue
		if '#' in line:
			continue
		i = line.split()
		pos,ref,alt,INFO = i[1],i[3],i[4],i[7]
		ANN = re.match(r'.+ANN=(.+)',INFO).group(1).split(',')[0].split('|')
		anno_out = '%s\t%s\t%s\t%s'%('2019-nCoV',pos,ref,alt)
		for index in anno_list:
			if not ANN[index]:
				anno_out += '\t.'
			else:
				anno_out += '\t'+ANN[index]
		if ANN[14]:
			anno_out += ';DISTANCE='+ANN[14]
		print(anno_out)
