# SARS-CoV-2_Multi-PCR_v1.0
SARS-CoV-2 analysis pipeline for multiplex-PCR MPS(Massive Parrallel Sequencing) data.

## Introduction
This pipeline could accurately and efficiently identify SARS-CoV-2 reads from multiplex PCR sequencing data, and report the infection status of sequencing samples with positive/negative/uncertain label. The pipeline could also get the variant information such as SNP/INDEL and generate the consensus sequence.

![Image](https://github.com/MGI-tech-bioinformatics/SARS-CoV-2_Multi-PCR_v1.0/blob/master/Pipeline.png)

## Updates:
May 11, 2020
1. Adjust the min depth threshold of freebayes from 100 to 30  
2. Updated Cut_Multi_Primer.py to save more memory  
3. Fixed some errors in the HTML report
4. Fixed a bug in consensus fasta

May 26, 2020
1. Added 'SOAPnuke_param' in the json file,users users can now customize the parameters of SOAPnuke.  
2. Fixed a bug in consensus fasta,which cause an error during the generation of consensus sequence when there is an INDEL in vcf file.  
3. Fixed a bug in generate_rem_report.py,which caused Identification.txt to display abnormally in the HTML report.  

Jun 2, 2020
1. Fixed a bug in Windows.Depth.svg

Jun 8, 2020
1. Prepared install.sh, users can install required software by running this script.  
2. Optimized Cut_Multi_Primer.py, this script only keeps virus reads now, which makes the primer cut step to run more efficiently.  

Jun 24, 2020
1. Fixed a bug in step6, we use zcat to read .gz files now.  

Jul 14, 2020
1. Fixed a bug in Cut_Multi_Primer.py
2. One sample corresponds to multiple barcodes, pipeline will merge the fastq files.  
3. Add some statistical result in Identification.txt

Aug 21, 2020
1. Apply variant annotation by snpEff  

Feb 24, 2021
1. Use bwa-mem instead of bwa-aln in alignment 
2. Update SARS-CoV-2 positive criteria: SARS-CoV-2 reads pct >= 0.1% AND (>= 1X Coverage ) >= 1% 
3. Update Freebayes version: v1.3.4

May 7, 2021
1. Use variant annotation excel instead of VCF file in HTML report
2. Optimized depth distribution SVG in HTML report.
3. Mark the primer base quality as 0 instead of removing primer sequence
4. Update primer sequence information
5. Reduce software running time
6. Upload a docker version of this software

## Requirements:
Before running this pipeline, you need to make sure that several pieces of software and/or modules are installed on the system:  

Perl: >=v5.22.0  
Python: >=v3.4.3  
R: >=v3.3.2

Library for Python3 and R:  
* Python3: pysam,pandas,openpyxl  
* R: Cairo

Softwares for data quality control:  
* seqtk v1.2 (https://github.com/lh3/seqtk)
* SOAPnuke v1.5.6 (https://github.com/BGI-flexlab/SOAPnuke)  

Software for alignment and bam file statistics:
* BWA v0.7.16 (https://github.com/lh3/bwa)
* Samtools v1.3 (https://github.com/samtools/samtools)
* bamdst v1.0.6 (https://github.com/shiquan/bamdst)  

Software for variant calling:  
* freebayes v1.3.0 (https://github.com/ekg/freebayes)  

Other required softwares:  
* bedtools v2.26.0 (https://bedtools.readthedocs.io/en/latest/)
* bcftools v1.6 (https://github.com/samtools/bcftools/)
* tabix v1.9 (https://github.com/samtools/tabix/)
* bgzip v1.9 (https://github.com/samtools/tabix/)
* mosdepth v0.2.9 (https://github.com/brentp/mosdepth)

## Installation

To clone the repository:

    git clone https://github.com/MGI-tech-bioinformatics/SARS-CoV-2_Multi-PCR_v1.0.git

To install the required software:

    cd SARS-CoV-2_Multi-PCR_v1.0; sh install.sh


Notes: 
* The install.sh will install the required software to SARS-CoV-2_Multi-PCR_v1.0/tools, if the software are not working,the above dependent software needs to be installed separately according to their instructions. The Perl,Python,R and their library need to be installed by users. After installing, users should edit the input.json file and change the software path to your own path. 
* This software is adapted to the MGI product *ATOPlex RNA Library Prep Set* which includes different versions depending on the primers, Current versions are as follows:
  ATOPlex RNA Library Prep. Set 1000027431, V2.0
  ATOPlex RNA Library Prep. Set 1000023556, V1.1
* According to the different versions of the kit, you need to specify the corresponding version information in the json file through the *primer_version* field.[1.1/2.0]


## Usage
### 1.Prepare the input.json file
The details for input.json file are as follows:
* FqType, sequencing type(PE100/SE50...). 
* sample_list, sample list file(sample_name/barcode_information/data_path).  
* workdir, analysis result directory.  
* SplitData, downsampling size of each sample(1G/1M/1K). 
* SOAPnuke_param, param of SOAPnuke.
* freebayes_param, param of freebayes.In particular,the parameter '-H -p 1' is necessary. 
* consensus_depth, threshold of point depth for consensus sequence.[1~20.Default:10] 
* python3, path to python3. 
* python3_lib, path to python3 library. 
* Rscript, path to Rscript. 
* R_lib, path to R library. 
* primer_version, version informations of data.[1.1/2.0.Default:2.0]
* tools(bwa,samtools....), path to this tool. 

### 2.Run the pipeline.
```
python3 Main_SARS-CoV-2.py -i input.json 
cd path/to/workdir
nohup sh main.sh &
```
### 3.Analysis result.
1.Quality control result
```
path/to/workdir/result/*/05.Stat/QC.xlsx
```
2.Identification result
```
path/to/workdir/result/*/05.Stat/Identification.xlsx
```
3.Variant calling result
```
path/to/workdir/result/*/05.Stat/*.vcf.gz
path/to/workdir/result/*/05.Stat/*.snpEff.anno.xlsx
```
4.HTML report
```
path/to/workdir/result/*/05.Stat/*.html
```
## With Docker

To pull a docker repository:  

    docker pull meizhiying/mgi-sars-cov-2:v1.2

Running  

	docker run -d \
	--name $WORKNAME \
	-v $workdir:$workdir \
	-v $datadir:$datadir \
	meizhiying/mgi-sars-cov-2:v1.2 \
	/SARS-CoV-2_pipeline/bin/python3/bin/python3 /SARS-CoV-2_pipeline/bin/Main_SARS-CoV-2.mgi_use.py -i $json

Notes
1. All requirements and software are install in docker image, there is no need to configure the software path in json file.  

	Json Demo: 
```json
    {
      "FqType": "PE", 
      "sample_list": "/data/MGI-SARS-CoV-2_V1.2/sample.list", 
      "workdir": "/data/MGI-SARS-CoV-2_V1.2/analysis", 
      "SplitData": "1M", 
      "freebayes_param": "-H -p 1 -q 20 -m 60 --min-coverage 20 -F 0.6", 
      "consensus_depth": "10" 
    }
```
2. $workdir is defined in json with "workdir", $datadir is defined in sample.list
