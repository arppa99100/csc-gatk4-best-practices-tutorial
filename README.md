# csc-gatk4-best-practices-tutorial.  
(Data Pre-processing + Germline SNPs+Indels workflows, optional: seq-format-conversion workflow).  

Introduction and workflows for Germline short variant discovery (SNPs + Indels) with GATK4 on CSC cPouta VM. Start from Fastq or ubam files to VCF file. 

This document offering the step-by-step introduction about calling germline short variants with GATK4.0.4.0 (and GATK4.0.6.0) Docker images, Samtools and Python 2.7 Docker images and Cromwell version 33 on CSC cPouta VM. It mainly contains four wdl scripts and matched input json files: paired-fastq-to-unmapped-bam.wdl (if needed), processing-for-variant-discovery-gatk4.wdl, haplotypecaller-gvcf-gatk4.wdl and joint-discovery-gatk4-local.wdl.  

The main steps´description can be checked from [GATK Best Practice](https://software.broadinstitute.org/gatk/best-practices/workflow?id=11165) and its matched [Github Page](https://github.com/gatk-workflows).  

In this [CSC Github Page](https://github.com/github-fish/csc-gatk4-best-practices-tutorial), we offers workflow scripts which suitable for running on cPouta VM. Tested with sample NA12878 on cPouta VM, 24 cores, 117.2GB RAM. More benchmarking scrips and results (**Sample NA12878, NA12891 and NA12892**) for different pipelines are seperately displayed in repo: [csc-seq-format-conversion](https://github.com/github-fish/csc-seq-format-conversion), [csc-gatk4-data-pre-processing](https://github.com/github-fish/csc-gatk4-data-pre-processing) and [csc-gatk4-germline-SNPs-Indels](https://github.com/github-fish/csc-gatk4-germline-SNPs-Indels).  
***
## Step by step tutorial to run this workflow (Example: NA12878 30x wgs) on cPouta VM. 
Containing 5 pipelines:  
1. Linux command line to divide FASTQ file according flowcell:lane (RG tag).  
2. paired-fastq-to-unmapped-bam.wdl + paired-fastq-to-unmapped-bam.inputs.NA12878.json.  
3. processing-for-variant-discovery-gatk4.wdl + processing-for-variant-discovery-gatk4.hg38.wgs.inputs.NA12878.json.  
4. haplotypecaller-gvcf-gatk4.wdl + haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.NA12878.json.  
5. joint-discovery-gatk4-local.wdl + joint-discovery-gatk4-local.hg38.wgs.inputs.json. 
***

### Tool Prerequisites    
Virtul machine (VM) on CSC cPouta    
Docker-ce: docker-ce-17.06.0.ce-1.el7.centos.x86_64.rpm     
Git  

### Used testing VM in cPOUTA   
* Flavor Name: hpc-gen2.24core / CentOS 7
*  RAM 117.2GB 
* VCPUs 24
* Swap 99GB (necessary),  can be adjust to smaller according your analysis scale   
 **Important Tips**:
1. To ensure docker can mount file with hard link which can save tons of space for you, store inputs data and cromwell execution directory in the Same Filesystem! 
Run the pipeline with sudo right. Unless, you will find that the intermediate files generate by docker with root right, then hard link will be invalid.
(Another possible solve method: https://gatkforums.broadinstitute.org/wdl/discussion/9477/localization-via-hard-link-has-failed)
2. Swap is necessary. I used 99G for test whole workflow. Its size can be adjust to smaller according your analysis scale.
### Running Instructions  

#### 1. Linux command line (If needed)
It separates FASTQ file to several sub-FASTQ files according different `flowcell:line` (ReadGroup tag).
```bash
mkdir /media/volume/data/fastq # Make a directory to run workflows then change into that directory.
cd /media/volume/data/fastq
mkdir /media/volume/data/fastq/NA12878 # Make a directory to store input files of sample NA12878.
```
Download fastq files of NA12878 from 1000 genomes project.
```bash
cd /media/volume/data/fastq/NA12878
wget ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR194/­ERR194147/­ERR194147_1.­fastq.­gz
wget ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­fastq/­ERR194/­ERR194147/­ERR194147_2.­fastq.­gz
```
Change file permission.
```bash
sudo chmod +w /media/volume/data/fastq/NA12878/*
sudo chmod +x /media/volume/data/fastq/NA12878/*
# Or, sudo chmod ugo+rwx /media/volume/data/fastq/NA12878/*
```
Grep then write all flowcell:lane tags into flowcell_lanes_NA12878.txt.
```bash
zcat ERR194147_1.fastq.gz |grep C0D8DACXX |cut -d ":" -f 3-4 | sort |uniq > flowcell_lanes_NA12878.txt
```
Separate fastq files according its flowcell_lane tags. Needed for pipeline #2. This flowcell_lane tag will be ReadGroup Tag in later analysis-ready BAM file.
```bash
zcat ERR194147_1.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:1 > ERR194147_C0D8DACXX.1_1.fastq
zcat ERR194147_2.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:1 > ERR194147_C0D8DACXX.1_2.fastq
zcat ERR194147_1.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:2 > ERR194147_C0D8DACXX.2_1.fastq
zcat ERR194147_2.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:2 > ERR194147_C0D8DACXX.2_2.fastq
zcat ERR194147_1.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:3 > ERR194147_C0D8DACXX.3_1.fastq
zcat ERR194147_2.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:3 > ERR194147_C0D8DACXX.3_2.fastq
zcat ERR194147_1.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:4 > ERR194147_C0D8DACXX.4_1.fastq
zcat ERR194147_2.fastq.gz |grep -A 3 --no-group-separator  C0D8DACXX:4 > ERR194147_C0D8DACXX.4_2.fastq
```
Count reads number to ensure this process not loss data. (Optional)
```bash
# (print number of lines == number of reads * 4)
$ wc -l ERR194147_C0D8DACXX.1_1.fastq
790978820 ERR194147_C0D8DACXX.1_1.fastq
$ wc -l ERR194147_C0D8DACXX.1_2.fastq
790978820 ERR194147_C0D8DACXX.1_2.fastq

$ wc -l ERR194147_C0D8DACXX.2_1.fastq
781019996 ERR194147_C0D8DACXX.2_1.fastq
$ wc -l ERR194147_C0D8DACXX.2_2.fastq
781019996 ERR194147_C0D8DACXX.2_2.fastq

$ wc -l ERR194147_C0D8DACXX.3_1.fastq
790165312 ERR194147_C0D8DACXX.3_1.fastq
$ wc -l ERR194147_C0D8DACXX.3_2.fastq
790165312 ERR194147_C0D8DACXX.3_2.fastq

$ wc -l ERR194147_C0D8DACXX.4_1.fastq
786896308 ERR194147_C0D8DACXX.4_1.fastq
$ wc -l ERR194147_C0D8DACXX.4_2.fastq
786896308 ERR194147_C0D8DACXX.4_2.fastq

$ zcat ERR194147_1.fastq.gz | echo $((`wc -l`/4))
787265109
# Number of reads == (790978820+781019996+790165312+786896308)/4 == 787,265,109. 
```
After checking, original ERR194147_1.fastq.gz and ERR194147_2.fastq.gz files can be deleted to save storage space. 
```bash
rm -rf ERR194147_1.fastq.gz ERR194147_2.fastq.gz 
```
Compress previous made separated fastq files to save storage space. (Optional)
```bash
gzip ERR194147_C0D8DACXX.1_1.fastq
gzip ERR194147_C0D8DACXX.1_2.fastq
gzip ERR194147_C0D8DACXX.2_1.fastq
gzip ERR194147_C0D8DACXX.2_2.fastq
gzip ERR194147_C0D8DACXX.3_1.fastq
gzip ERR194147_C0D8DACXX.3_2.fastq
gzip ERR194147_C0D8DACXX.4_1.fastq
gzip ERR194147_C0D8DACXX.4_2.fastq
```

#### 2. paired-fastq-to-unmapped-bam
This WDL convert paired FASTQ to uBAM and add RG (read group) tag. uBAM is the input format of pipeline #3, GATK Data Pre-processing workflow of Best Practices.

**Inputs**:  
Paired-end sequencing data in FASTQ format, already separated via sequence flowcell:lane and orientation.

**Outputs**:  
Set of unmapped BAMs, one per read group  
A unmapped_bam.list file containing a list of the generated unmapped BAMs

Run pipeline #2 to generate a set of ubam files and a NA12878_unmapped_bam.list file.
Download Cromwell, the java excutable that will run the WDL.
```bash
cd /media/volume/
wget https://github.com/broadinstitute/cromwell/releases/download/33.1/cromwell-33.1.jar
```
Setup your working directory and make a directory to store your input files.
```bash
mkdir /media/volume/gatk4-best-practices
cd /media/volume/gatk4-best-practices
```
Git clone needed github repo.
```bash
git clone git://github.com/github-fish/csc-gatk4-best-practices-tutorial
```  
Edit the paired-fastq-to-unmapped-bam.inputs.NA12878.json file so that all file paths are replaced by your VM local file paths.   

Run the pipeline with sudo right.   
```bash
# Just an example, here. My tesing work directory was different. You can get detailed info via the TestingResult folder of different repo mentioned ahead. [csc-seq-format-conversion](https://github.com/github-fish/csc-seq-format-conversion)
$ cd /media/volume/data/fastq/NA12878
$ sudo java -jar /media/volume/cromwell-33.1.jar run /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/paired-fastq-to-unmapped-bam.wdl -i /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/paired-fastq-to-unmapped-bam.inputs.NA12878.json
# if you want to monitor the whole process time costing, add "time" parameter after "sudo".
```
Cromwell will generate two folders in current working directory (/media/volume/data/fastq/NA12878) named: cromwell-executions and cromwell-workflow-logs which separately store pipeline outputs, intermediate files and cromwell-worklog file (if this pipeline is not successfully finished). Otherwise, the cromwell-workflow-logs folder will be empty.

#### 3. processing-for-variant-discovery-gatk4
The processing-for-variant-discovery-gatk4 WDL pipeline implements data pre-processing according to the Board Institute GATK Best Practices (June 2016).

**Inputs**:  
Pair-end sequencing data in unmapped BAM (uBAM) format  
One or more read groups, one per uBAM file, all belonging to a single sample (SM)  
Input uBAM files must additionally comply with the following requirements:  
1. filenames all have the same suffix (we use “.unmapped.bam”). 
2. files must pass validation by ValidateSamFile. 
3. reads are provided in query-sorted order. 
4. all reads must have RG tag

**Outputs**:  
Set of unmapped BAMs, one per read group.   
A unmapped_bam.list file containing a list of the generated unmapped BAMs.  

Make a directory to run workflows and a directory to store required input files.
```bash
mkdir /media/volume/data-pro
mkdir /media/volume/data-pro/inputs
mkdir /media/volume/data-pro/inputs/NA12878_uBAM
```
Download/prepare all needed input files declared in processing-for-variant-discovery-gatk4.hg38.wgs.inputs.NA12878.json file, store them in /media/volume/data-pro/inputs and /media/volume/data-pro/inputs/NA12878_uBAM folder. Edit the json file so that all file paths are replaced by your VM local file paths. [Google Bucket](https://console.cloud.google.com/storage/browser/broad-references/hg38/v0).      

Run the pipeline with sudo right.  
```bash
$ mkdir /media/volume/data-pro/NA12878
$ cd /media/volume/data-pro/NA12878
$ sudo java -jar /media/volume/cromwell-33.1.jar run /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/processing-for-variant-discovery-gatk4.wdl -i /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/processing-for-variant-discovery-gatk4.hg38.wgs.inputs.NA12878.json
```

#### 4. haplotypecaller-gvcf-gatk4
The haplotypecaller-gvcf-gatk4 workflow runs HaplotypeCaller from GATK4 in GVCF mode on a single sample according to the Board Institute GATK Best Practices (June 2016), scattered across intervals.

**Inputs**:  
One analysis-ready BAM file for a single sample (as identified in RG:SM). 
Set of variant calling intervals lists for the scatter, provided in a file. 

**Outputs**:  
One GVCF file and its index  

Run pipeline #4 to generate one g.vcf file for each analysis-ready BAM file or CRAM file.   
Make a directory to run workflows and a directory to store required input files.   
```bash
mkdir /media/volume/g-snps-indels/
mkdir /media/volume/g-snps-indels/inputs
```
Download/prepare all needed input files declared in json file, store them in /media/volume/g-snps-indels/inputs folder. Edit the haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.json file so that all file paths are replaced by your VM local file paths.    
Run the pipeline with sudo right. [Google Bucket1](https://console.cloud.google.com/storage/browser/broad-references/hg38/v0) and [Google Bucket2](https://console.cloud.google.com/storage/browser/gatk-test-data/intervals)
```bash
$ cd /media/volume/g-snps-indels
$ sudo java -jar /media/volume/cromwell-33.1.jar run /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/haplotypecaller-gvcf-gatk4.wdl -i /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/haplotypecaller-gvcf-gatk4.hg38.wgs.inputs.json
```
#### 5. joint-discovery-gatk4-local
The joint-discovery-gatk4-local.wdl implements the joint discovery and VQSR filtering portion of the GATK Best Practices (June 2016) for germline SNP and Indel discovery in human whole-genome sequencing (WGS) and exome sequencing data.

**Inputs**:  
One or more GVCFs produced by HaplotypeCaller in GVCF mode.  
Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.   
When deteriming disk size in the json, use the guideline below.  
    small_disk = (num_gvcfs / 10) + 10.  
    medium_disk = (num_gvcfs * 15) + 10.  
    huge_disk = num_gvcfs + 10.  

**Outputs**:  
A VCF file and its index, filtered using variant quality score recalibration (VQSR) with genotypes for all samples present in the input VCF. All sites that are present in the input VCF are retained; filtered sites are annotated as such in the FILTER field.   

run pipeline #5 to generate one VCF file and its index file from one or several g.vcf files. (Here use one NA12878 g.vcf file as example).  

Download/prepare all needed input files declared in json file, store them in /media/volume/g-snps-indels/inputs folder. Edit the joint-discovery-gatk4-local.hg38.wgs.inputs.json file so that all file paths are replaced by your VM local file paths. [Google Bucket1](https://console.cloud.google.com/storage/browser/broad-references/hg38/v0) and [Google Bucket2](https://console.cloud.google.com/storage/browser/gatk-test-data/intervals).     
       
Run the pipeline with sudo right.  
```bash
$ mkdir /media/volume/g-snps-indels/joint_calling
$ cd /media/volume/g-snps-indels/joint_calling
$ sudo java -Dsystem.input-read-limits.lines=500000 -jar /media/volume/cromwell-33.1.jar run /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/joint-discovery-gatk4-local.wdl -i /media/volume/gatk4-best-practices/csc-gatk4-best-practices-tutorial/joint-discovery-gatk4-local.hg38.wgs.inputs.json
```
**Important Tips**:
For this pipeline, because intervals-hg38.even.handcurated.20k.intervals is larger than 128000 Bytes, the default Dsystem.input-read-limits.lines is not enough.
Add `-Dsystem.input-read-limits.lines=500000` in your command line or add it in your config file. 
Here is the Cromwell Reference Config File which contains all the default _override_ settings for cromwell. (Latest update 14 Aug 2018)
 https://github.com/broadinstitute/cromwell/blob/0ad9cb61f90276e3576f6ba135cab65c02f09a9e/core/src/main/resources/reference.conf#L157