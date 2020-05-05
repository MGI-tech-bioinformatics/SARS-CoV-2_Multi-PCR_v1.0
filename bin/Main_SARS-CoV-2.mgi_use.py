#!/usr/bin/env python3

import sys
import os
import re
import getopt
import configparser
import json
from subprocess import check_call

def usage():
	"""Main program of SARS-CoV-2 panel analysis pipeline.
	Usage: Main_SARS-CoV-2.py -i <config>
	-i	Config json file
	-h	Help information
	-s	Submit model [local/qsubsge]
	"""
	print(usage.__doc__)
	return

def create_dirs(*dirname):
	for each in dirname:
		if not os.path.exists(each):
			os.system("mkdir -p %s" % each)
	return

def fqtype_error():
	print('FqType is not valid,use SE/PE.')
	sys.exit(1)
	return

def CleanData(script,sample,barcode,fqtype,SplitData):
	barcode_dir = result_dir + '/' + sample
	Clean_dir = barcode_dir + '/01.Clean'
	create_dirs(Clean_dir)
	if fqtype == 'SE':
		rawfq = raw_data_path + '/*' + barcode + '.fq.gz'
		cleanfq = Clean_dir + '/Clean_' + sample + '.fq.gz'

		if SplitData:
			t_dict = {'G':10**9,'M':10**6,'K':10**3}
			SplitData_n = float(SplitData.strip()[0:-1])*t_dict[SplitData.strip()[-1]]
			splitfq = Clean_dir + '/Split_' + sample + '.fq'
			script.write("%(seqtk)s sample -s100 %(rawfq)s %(SplitData_n)s > %(splitfq)s && gzip %(splitfq)s && "\
				%{'seqtk':seqtk,'rawfq':rawfq,'splitfq':splitfq,'SplitData_n':SplitData_n})
			script.write("%(SOAPnuke)s filter -l 10 -q 0.2 -n 0.05 -Q 2 -G -T 1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -1 %(splitfq)s.gz -C %(cleanfq)s -o %(Clean_dir)s \n"\
				%{'SOAPnuke':SOAPnuke,'splitfq':splitfq,'cleanfq':cleanfq,'Clean_dir':Clean_dir})
		else:
			script.write("%(SOAPnuke)s filter -l 10 -q 0.2 -n 0.05 -Q 2 -G -T 1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -1 %(rawfq)s -C %(cleanfq)s -o %(Clean_dir)s \n"\
				%{'SOAPnuke':SOAPnuke,'rawfq':rawfq,'cleanfq':cleanfq,'Clean_dir':Clean_dir})
	elif fqtype == 'PE':
		rawfq1 = raw_data_path + '/*' + barcode + '_1.fq.gz'
		rawfq2 = raw_data_path + '/*' + barcode + '_2.fq.gz'
		cleanfq1 = Clean_dir + '/Clean_' + sample + '_1.fq.gz'
		cleanfq2 = Clean_dir + '/Clean_' + sample + '_2.fq.gz'

		if SplitData:
			t_dict = {'G':10**9,'M':10**6,'K':10**3}
			SplitData_n = float(SplitData.strip()[0:-1])*t_dict[SplitData.strip()[-1]]/2
			splitfq1 = Clean_dir + '/Split_' + sample + '_1.fq'
			splitfq2 = Clean_dir + '/Split_' + sample + '_2.fq'
			script.write("%(seqtk)s sample -s100 %(rawfq1)s %(SplitData_n)s > %(splitfq1)s && gzip %(splitfq1)s && %(seqtk)s sample -s100 %(rawfq2)s %(SplitData_n)s > %(splitfq2)s && gzip %(splitfq2)s && "\
				%{'seqtk':seqtk,'rawfq1':rawfq1,'rawfq2':rawfq2,'splitfq1':splitfq1,'splitfq2':splitfq2,'SplitData_n':SplitData_n})
			script.write("%(SOAPnuke)s filter -l 10 -q 0.2 -n 0.05 -Q 2 -G -T 1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA  -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -1 %(splitfq1)s.gz -2 %(splitfq2)s.gz -C %(cleanfq1)s -D %(cleanfq2)s -o %(Clean_dir)s \n"\
				%{'SOAPnuke':SOAPnuke,'splitfq1':splitfq1,'splitfq2':splitfq2,'cleanfq1':cleanfq1,'cleanfq2':cleanfq2,'Clean_dir':Clean_dir})
		else:
			script.write("%(SOAPnuke)s filter -l 10 -q 0.2 -n 0.05 -Q 2 -G -T 1 -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA  -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG -1 %(rawfq1)s -2 %(rawfq2)s -C %(cleanfq1)s -D %(cleanfq2)s -o %(Clean_dir)s \n"\
			%{'SOAPnuke':SOAPnuke,'rawfq1':rawfq1,'rawfq2':rawfq2,'cleanfq1':cleanfq1,'cleanfq2':cleanfq2,'Clean_dir':Clean_dir})
	else:
		fqtype_error()
	return

def bwaaln(script,barcode,fqtype,read_len):
	barcode_dir = result_dir + '/' + barcode
	Clean_dir = barcode_dir + '/01.Clean'
	Align_dir = barcode_dir + '/02.Align'
	create_dirs(Align_dir)
	seed_len = str(int(int(read_len)*0.9))
	if fqtype == 'SE':
		cleanfq = Clean_dir + '/Clean_' + barcode + '.fq.gz'
		script.write("%(bwa)s aln -l %(seed_len)s -t 3 -f %(Align_dir)s/%(barcode)s.sai %(database)s/nCoV.fa %(cleanfq)s && %(bwa)s samse -r \"@RG\\tID:PE100\\tPL:MGISEQ\\tPU:PE100\\tLB:mutPCR\\tSM:%(barcode)s\\tCN:BGI\" %(database)s/nCoV.fa %(Align_dir)s/%(barcode)s.sai %(cleanfq)s | %(samtools)s view -b - | %(samtools)s sort -T %(Align_dir)s/%(barcode)s.sort -o %(Align_dir)s/%(barcode)s.sort.bam - \n" \
			%{'bwa':bwa,'samtools':samtools,'database':database,'cleanfq':cleanfq,'Align_dir':Align_dir,'barcode':barcode,'seed_len':seed_len})
		script.write("%(samtools)s index %(Align_dir)s/%(barcode)s.sort.bam \n" %{'samtools':samtools,'Align_dir':Align_dir,'barcode':barcode})
	elif fqtype == 'PE':
		cleanfq1 = Clean_dir + '/Clean_' + barcode + '_1.fq.gz'
		cleanfq2 = Clean_dir + '/Clean_' + barcode + '_2.fq.gz'
		script.write("%(bwa)s aln -l %(seed_len)s -t 3 -f %(Align_dir)s/%(barcode)s_1.sai %(database)s/nCoV.fa %(cleanfq1)s && %(bwa)s aln -l %(seed_len)s -t 3 -f %(Align_dir)s/%(barcode)s_2.sai %(database)s/nCoV.fa %(cleanfq2)s && %(bwa)s sampe -a 1000 -r \"@RG\\tID:PE100\\tPL:MGISEQ\\tPU:PE100\\tLB:mutPCR\\tSM:%(barcode)s\\tCN:BGI\" %(database)s/nCoV.fa %(Align_dir)s/%(barcode)s_1.sai %(Align_dir)s/%(barcode)s_2.sai %(cleanfq1)s %(cleanfq2)s | %(samtools)s view -b - | %(samtools)s sort -T %(Align_dir)s/%(barcode)s.sort -o %(Align_dir)s/%(barcode)s.sort.bam - \n" \
			%{'bwa':bwa,'samtools':samtools,'database':database,'cleanfq1':cleanfq1,'cleanfq2':cleanfq2,'Align_dir':Align_dir,'barcode':barcode,'seed_len':seed_len})
		script.write("%(samtools)s index %(Align_dir)s/%(barcode)s.sort.bam \n" %{'samtools':samtools,'Align_dir':Align_dir,'barcode':barcode,'seed_len':seed_len})
	return

def CovDep(script,barcode):
	barcode_dir = result_dir + '/' + barcode
	Align_dir = barcode_dir + '/02.Align'
	covdep_dir = barcode_dir + '/03.covdep'
	lambda_covdep_dir = covdep_dir + '/lambda_cov'
	create_dirs(covdep_dir,lambda_covdep_dir)
	## bamqc with all map bam
	script.write("%(bamdst)s --cutoffdepth 1 --maxdepth 1000000 -q 1 -p %(virusbed)s -o %(covdep_dir)s %(Align_dir)s/%(barcode)s.sort.bam \n"\
		%{'bamdst':bamdst,'virusbed':virusbed,'covdep_dir':covdep_dir,'Align_dir':Align_dir,'barcode':barcode})
	script.write("%(bamdst)s --cutoffdepth 1 --maxdepth 1000000 -q 1 -p %(lambdabed)s -o %(covdep_dir)s/lambda_cov %(Align_dir)s/%(barcode)s.sort.bam \n"\
		%{'bamdst':bamdst,'lambdabed':lambdabed,'covdep_dir':covdep_dir,'Align_dir':Align_dir,'barcode':barcode})
	return

def Statistics(script,fqtype):
	if fqtype == 'PE':
		script.write("export PYTHONPATH=%(python3_lib)s:$PYTHONPATH && %(python3)s %(bin)s/CoV_stat.py -t PE -d %(result_dir)s -l %(barcode_file)s\n"\
			%{'bin':bin,'result_dir':result_dir,'barcode_file':barcode_file,'python3':python3,'python3_lib':python3_lib})
	elif fqtype == 'SE':
		script.write("export PYTHONPATH=%(python3_lib)s:$PYTHONPATH && %(python3)s %(bin)s/CoV_stat.py -t SE -d %(result_dir)s -l %(barcode_file)s\n"\
			%{'bin':bin,'result_dir':result_dir,'barcode_file':barcode_file,'python3':python3,'python3_lib':python3_lib})
	return

def CutPrimer(script,fqtype_p,sample):
	sample_dir = result_dir + '/' + sample
	Align_dir = sample_dir + '/02.Align'
	CutPrimer_dir = sample_dir + '/04.CutPrimer'
	create_dirs(CutPrimer_dir)
	script.write("export PYTHONPATH=%(python3_lib)s:$PYTHONPATH && %(python3)s %(bin)s/Cut_Multi_Primer.py -p %(primer_list)s -b %(Align_dir)s/%(sample)s.sort.bam -s %(sample)s -o %(CutPrimer_dir)s -t %(fqtype_p)s \n"\
		%{'bin':bin,'primer_list':primer_list,'Align_dir':Align_dir,'sample':sample,'CutPrimer_dir':CutPrimer_dir,'lib':lib,'python3':python3,'fqtype_p':fqtype_p,'python3_lib':python3_lib})
	return

def AlignVariant(script,fqtype,cutprimer_list,consensus_depth):
	with open(cutprimer_list,'r') as fi:
		for line in fi:
			i = line.split()
			sample = i[0]
			if fqtype == 'PE':
				fq1 = i[1]
				fq2 = i[2]
			elif fqtype == 'SE':
				fq = i[1]
			sample_dir = result_dir + '/' + sample
			Align_dir = sample_dir + '/02.Align'
			CutPrimer_dir = sample_dir + '/04.CutPrimer'
			Stat_dir = sample_dir + '/05.Stat'
			create_dirs(Stat_dir)
			if fqtype == 'PE':
				script.write("%(bwa)s mem -Y -M -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\" -t 3 %(ref)s %(fq1)s %(fq2)s | %(samtools)s view -b - | %(samtools)s sort -T %(Stat_dir)s/%(sample)s -o %(Stat_dir)s/%(sample)s.bam -\n"\
					%{'bwa':bwa,'samtools':'samtools','sample':sample,'ref':ref,'fq1':fq1,'fq2':fq2,'Stat_dir':Stat_dir})
			elif fqtype == 'SE':
				script.write("%(bwa)s mem -Y -M -R \"@RG\\tID:%(sample)s\\tSM:%(sample)s\" -t 3 %(ref)s %(fq)s | %(samtools)s view -b - | %(samtools)s sort -T %(Stat_dir)s/%(sample)s -o %(Stat_dir)s/%(sample)s.bam -\n"\
					%{'bwa':bwa,'samtools':'samtools','sample':sample,'ref':ref,'fq':fq,'Stat_dir':Stat_dir})
			script.write("%(samtools)s index %(Stat_dir)s/%(sample)s.bam\n"%{'samtools':samtools,'Stat_dir':Stat_dir,'sample':sample})
			script.write("%(mosdepth)s -n --fast-mode --by 100 %(Stat_dir)s/depth %(Stat_dir)s/%(sample)s.bam\n"%{'mosdepth':mosdepth,'Stat_dir':Stat_dir,'sample':sample})
			script.write("less %(Stat_dir)s/depth.regions.bed.gz|awk  '{print NR\"\\t\"log($4+0.1)}' > %(Stat_dir)s/%(sample)s.draw.depth\n"%{'Stat_dir':Stat_dir,'sample':sample})
			script.write("%(samtools)s depth -d 100000000 -a -b %(virusbed)s %(Align_dir)s/%(sample)s.sort.bam > %(Stat_dir)s/%(sample)s.depth\n"%{'samtools':samtools,'virusbed':virusbed,'Align_dir':Align_dir,'sample':sample,'Stat_dir':Stat_dir})
			script.write("export LD_LIBRARY_PATH=/ldfssz1/MGI_BIT/RUO/meizhiying/project/Multi_PCR_2019-nCoV/SARS-CoV-2_pipeline/bin/AlignVariant/../../lib/lib64:$LD_LIBRARY_PATH && export R_LIBS=%(R_lib)s:$R_LIBS && %(Rscript)s %(bin)s/line.depth.R %(Stat_dir)s/%(sample)s.draw.depth %(Stat_dir)s/Windows.Depth.svg\n"%{'Rscript':Rscript,'bin':bin,'Stat_dir':Stat_dir,'R_lib':R_lib,'sample':sample})
			script.write("%(bin)s/Consensus.pl %(Stat_dir)s/%(sample)s.depth %(ref)s %(consensus_depth)s %(Stat_dir)s/%(sample)s.reference1.fa\n"%{'bin':bin,'Stat_dir':Stat_dir,'sample':sample,'ref':ref,'consensus_depth':consensus_depth})
			script.write("export LD_LIBRARY_PATH=/ldfssz1/MGI_BIT/RUO/meizhiying/project/Multi_PCR_2019-nCoV/SARS-CoV-2_pipeline/bin/AlignVariant/../../lib/anaconda2/lib:$LD_LIBRARY_PATH && %(freebayes)s -t %(variantbed)s %(freebayes_param)s -f %(ref)s %(Stat_dir)s/%(sample)s.bam > %(Stat_dir)s/%(sample)s.raw.vcf && %(bcftools)s view --include 'FMT/GT=\"1\" && QUAL>=100 && FMT/DP>=100' %(Stat_dir)s/%(sample)s.raw.vcf > %(Stat_dir)s/%(sample)s.vcf\n%(bgzip)s -f %(Stat_dir)s/%(sample)s.vcf\n%(tabix)s %(Stat_dir)s/%(sample)s.vcf.gz\n"%{'freebayes':freebayes,'variantbed':variantbed,'Stat_dir':Stat_dir,'sample':sample,'ref':ref,'bgzip':bgzip,'tabix':tabix,'bcftools':bcftools,'freebayes_param':freebayes_param})
			script.write("%(bcftools)s consensus -f %(Stat_dir)s/%(sample)s.reference1.fa -o %(Stat_dir)s/%(sample)s.Consensus.fa %(Stat_dir)s/%(sample)s.vcf.gz\n"%{'bcftools':bcftools,'Stat_dir':Stat_dir,'sample':sample})
			script.write("sed -i \"s/MN908947.3 Wuhan seafood market pneumonia virus isolate Wuhan-Hu-1, complete genome/%(sample)s/g\" %(Stat_dir)s/%(sample)s.Consensus.fa\n"%{'Stat_dir':Stat_dir,'sample':sample})
			script.write("less %(Stat_dir)s/%(sample)s.vcf.gz | grep -v '^#'|awk '{print $1\"\\t\"$2-1\"\\t\"$2\"\\t\"$4\"\\t\"$5}'| %(bedtools)s intersect -a - -b %(bed2)s -loj |cut -f 1,3-5,9 > %(Stat_dir)s/%(sample)s.vcf.anno\n"%{'Stat_dir':Stat_dir,'sample':sample,'bedtools':bedtools,'bed2':bed2})
			script.write("rm %(Stat_dir)s/%(sample)s.draw.depth %(Stat_dir)s/%(sample)s.reference1.fa\n"%{'Stat_dir':Stat_dir,'sample':sample})
	return

def GetReport(script,sample):
	sample_dir = result_dir + '/' + sample
	Stat_dir = sample_dir + '/05.Stat'
	#script.write("perl %(bin)s/txt2excel %(Stat_dir)s/QC.txt %(Stat_dir)s/QC.xlsx\nperl %(bin)s/txt2excel %(Stat_dir)s/Identification.txt %(Stat_dir)s/Identification.xlsx\n%(python3)s %(bin)s/etiology/generate_rem_report.py %(Stat_dir)s %(sample)s \n"\
	script.write("%(python3)s %(bin)s/etiology/generate_rem_report.py %(Stat_dir)s %(sample)s \n"\
		%{'bin':bin,'Stat_dir':Stat_dir,'sample':sample,'python3':python3})
	return

def MainShell(script_file,step1shell,step2shell,step3shell,step4shell,step5shell,step6shell,step7shell):
	script = open(script_file,'w')
	script.write('''echo "start step1  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 1G --lines 1 --maxjob 300 %(stepshell)s && echo "finish step1 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step1shell})
	script.write('''echo "start step2  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 6G --lines 2 --maxjob 300 %(stepshell)s && echo "finish step2 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step2shell})
	script.write('''echo "start step3  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 1G --lines 1 --maxjob 300 %(stepshell)s && echo "finish step3 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step3shell})
	script.write('''echo "start step4  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 1G --lines 1 --maxjob 300 %(stepshell)s && echo "finish step4 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step4shell})
	script.write('''echo "start step5  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 1G --lines 1 --maxjob 300 %(stepshell)s && echo "finish step5 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step5shell})
	script.write('''echo "start step6  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 6G --lines 14 --maxjob 300 %(stepshell)s && echo "finish step6 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step6shell})
	script.write('''echo "start step7  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(watchdog)s --mem 1G --lines 1 --maxjob 300 %(stepshell)s && echo "finish step7 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'watchdog':watchdog,'stepshell':step7shell})
	script.close()
	return

def MainShell_qsubsge(script_file,step1shell,step2shell,step3shell,step4shell,step5shell,step6shell,step7shell):
	script = open(script_file,'w')
	script.write('''echo "start step1  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=1"  --jobprefix step1 --lines 1 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step1 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step1shell})
	script.write('''echo "start step2  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=3"  --jobprefix step2 --lines 2 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step2 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step2shell})
	script.write('''echo "start step3  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=1"  --jobprefix step3 --lines 1 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step3 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step3shell})
	script.write('''echo "start step4  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=1"  --jobprefix step4 --lines 1 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step4 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step4shell})
	script.write('''echo "start step5  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=1"  --jobprefix step5 --lines 1 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step5 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step5shell})
	script.write('''echo "start step6  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=6G -P %(subproject)s -l num_proc=1"  --jobprefix step6 --lines 14 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step6 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step6shell})
	script.write('''echo "start step7  at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`" && perl %(qsubsge)s --queue %(queue)s --resource="vf=1G -P %(subproject)s -l num_proc=1"  --jobprefix step7 --lines 1 --reqsub --interval 5 --convert no -maxjob 500 %(stepshell)s && echo "finish step7 at `date +'%%Y-%%m-%%d %%H:%%M:%%S %%z'`"\n'''\
		%{'qsubsge':qsubsge,'queue':queue,'subproject':subproject,'stepshell':step7shell})
	script.close()
	return

if __name__ == '__main__':
	if len(sys.argv) < 2:
		usage()
		sys.exit(1)

	try:
		opts,args = getopt.getopt(sys.argv[1:],"hi:s:")
	except getopt.GetoptError:
		print("ERROR: Param error")
		sys.exit(1)

	subprj = 'local'

	for key, value in opts:
		if key == '-h':
			usage()
			sys.exit()
		if key == '-i':
			jsonfile = value
		if key == '-s':
			subprj = value
	file_json = open(jsonfile,'r')
	jsonobj = json.load(file_json)
	rootpath = os.path.dirname(sys.path[0])
	bin = rootpath + '/bin'
	lib = rootpath + '/lib'
	database = rootpath + '/database'
	tools = rootpath + '/tools'

	# tools
	try:
		freebayes_param = jsonobj["freebayes_param"]
	except:
		freebayes_param = '-p 1 -q 20 -m 60 --min-coverage 100'
	python3 = jsonobj["python3"]
	python3_lib = jsonobj["python3_lib"]
	Rscript = jsonobj["Rscript"]
	R_lib = jsonobj["R_lib"]
	seqtk = jsonobj["seqtk"]
	bwa = jsonobj["bwa"]
	samtools = jsonobj["samtools"]
	freebayes = jsonobj["freebayes"]
	bcftools = jsonobj["bcftools"]
	bgzip = jsonobj["bgzip"]
	tabix = jsonobj["tabix"]
	bedtools = jsonobj["bedtools"]
	mosdepth = jsonobj["mosdepth"]
	bamdst = jsonobj["bamdst"]
	SOAPnuke = jsonobj["SOAPnuke"]
	#seqtk = '%s/seqtk'%(tools)
	#bwa = '%s/bwa'%(tools)
	#samtools = '%s/samtools'%(tools)
	#freebayes = '%s/freebayes'%(tools)
	#bcftools = '%s/bcftools'%(tools)
	#bgzip = '%s/bgzip'%(tools)
	#tabix = '%s/tabix'%(tools)
	#bedtools = '%s/bedtools'%(tools)
	#mosdepth = '%s/mosdepth'%(tools)
	#bamdst = '%s/bamdst'%(tools)
	#SOAPnuke = '%s/SOAPnuke'%(tools)

	fqtype_p = jsonobj["FqType"]
	fqtype = fqtype_p[0:2]
	read_len = fqtype_p[2:]
	barcode_file = jsonobj["sample_list"]
	#snpbed = jsonobj["targetregion"]
	virusbed = database + '/nCoV.virus.bed'
	lambdabed = database + '/nCoV.lambda.bed'
	variantbed = database + '/nCoV.variant.bed'
	bed2 = database + '/wuhanRef.bed'
	ref = database + '/nCov.fasta'
	try:
		consensus_depth = jsonobj["consensus_depth"]
	except:
		consensus_depth = 100
	watchdog = bin+'/localsubmit/bin/watchDog_v1.0.pl'
	qsubsge = rootpath+'/bin/qsub-sge.pl'
	#queue = 'mgi.q'
	queue = jsonobj["queue"]
	#subproject = 'P18Z18000N0394'
	subproject = jsonobj["project"]
	work_dir = jsonobj["workdir"]
	primer_list = database + '/nCoV.primer.xls'
	try:
		SplitData = jsonobj["SplitData"]
	except:
		SplitData = ''
	result_dir = os.path.abspath(work_dir)+"/result"
	shell_dir = os.path.abspath(work_dir)+"/shell"
	create_dirs(work_dir,result_dir,shell_dir)

	step1shell = open(shell_dir+ '/step1.filter.sh','w')
	step2shell = open(shell_dir+ '/step2.bwa.sh','w')
	step3shell = open(shell_dir+ '/step3.bamdst.sh','w')
	step4shell = open(shell_dir+ '/step4.CutPrimer.sh','w')
	step5shell = open(shell_dir+ '/step5.statistic.sh','w')
	step6shell = open(shell_dir+ '/step6.AlignVariant.sh','w')
	step7shell = open(shell_dir+ '/step7.GetReport.sh','w')

	file_cut_primer_list = open('%s/CutPrimer.list'%(work_dir),'w')

	with open(barcode_file,'r') as fb:
		for line in fb:
			sample, barcode, raw_data_path = line.split()
			CleanData(step1shell,sample,barcode,fqtype,SplitData)
			bwaaln(step2shell,sample,fqtype,read_len)
			CovDep(step3shell,sample)
			CutPrimer(step4shell,fqtype_p,sample)
			GetReport(step7shell,sample)
			if fqtype == 'PE':
				file_cut_primer_list.write('%s\t%s/%s/04.CutPrimer/%s_1.cutprimer.fq.gz\t%s/%s/04.CutPrimer/%s_2.cutprimer.fq.gz\n'%(sample,result_dir,sample,sample,result_dir,sample,sample))
			elif fqtype == 'SE':
				file_cut_primer_list.write('%s\t%s/%s/04.CutPrimer/%s.cutprimer.fq.gz\n'%(sample,result_dir,sample,sample))
	file_cut_primer_list.close()
	Statistics(step5shell,fqtype)
	#check_call('perl %s/AlignVariant/Align.Variant.pl --list %s/CutPrimer.list --type %s --outdir %s'%(bin,work_dir,fqtype,result_dir),shell=True)
	AlignVariant(step6shell,fqtype,'%s/CutPrimer.list'%(work_dir),consensus_depth)
	step1shell.close()
	step2shell.close()
	step3shell.close()
	step4shell.close()
	step5shell.close()
	step6shell.close()
	step7shell.close()

	finalshell = work_dir + "/main.sh"
	if subprj == 'local':
		MainShell(finalshell,shell_dir+ '/step1.filter.sh',shell_dir+ '/step2.bwa.sh',shell_dir+ '/step3.bamdst.sh',shell_dir+ '/step4.CutPrimer.sh',shell_dir+ '/step5.statistic.sh',shell_dir+ '/step6.AlignVariant.sh',shell_dir+ '/step7.GetReport.sh')
	elif subprj == 'qsubsge':
		MainShell_qsubsge(finalshell,shell_dir+ '/step1.filter.sh',shell_dir+ '/step2.bwa.sh',shell_dir+ '/step3.bamdst.sh',shell_dir+ '/step4.CutPrimer.sh',shell_dir+ '/step5.statistic.sh',shell_dir+ '/step6.AlignVariant.sh',shell_dir+ '/step7.GetReport.sh')
	else:
		print('ERROR: invalid -s param,use local/qsubsge')
		sys.exit(1)
