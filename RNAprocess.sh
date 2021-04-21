#!/bin/bash

### RNA Seq Preprocessing Script for Paired-End Data ###

#--------------------
# Setting Working Directory and Inputs:
#--------------------

function usage () {
	echo
	echo 'RNA Seq PreProcessing Pipeline'
	echo
	echo 'to run:'
        echo '$bash RNAprocess.sh -c <path to config file>'
        echo
	echo 'Required Arguments:'		# TO DO: Its not really required at this time, you can do conf or arguments so either make it mandatory or fix --help instructions
        echo '-c : path to the config file'
	echo
	echo 'Optional Arguments - to be used in place of config file'
	echo '-a : accesion file name, contains the SRR values from NCBI that will be used for fastqdump'
        echo '-o : name of the input file that contains required input, will also be used for output'
	echo '-d : true to download fastq files from SRA accession file, otherwise false'
	echo '-p : true if data is paired end otherwise set to false'
	echo '-r : set to true if a reference genome index is required, otherwise false; designed for hg38 human reference genome'
	echo
}

# Arguments: Can be used in place of the .conf file
if [ "$#" == 0 ]; then
	echo
	echo "Please provide the path to a valid config file (preprocess.conf) or adjust arguments. For additional help use -h / --h" 
	echo
	exit
else
	while [ "$#" -gt 0 ]; do
	case $1 in
		-h|--h)
			usage; exit
			;;
		-c|--c)
			confFile=$2
			. $confFile
			shift
			;;
		-o|--o)
                	if [ -d "$2" ]; then
				inoutdir="$2"; echo "Input Directory : $inoutdir"
			else
				echo "Please provide a valid directory for the input folder (-o / --o)"; exit
			fi
			shift
			;;
		-a|--a) #TO DO: Add error if -a input is not a csv format
			accFileName=$2; echo "Accession File : $accFileName"
			shift
			;;
		-d|--d)
			download=$2
			shift
			;;
		-p|--p)
			pairedEnd=$2
			shift
			;;
		-r|--r)
			refGenomeIndex=$2
			shift
			;;

		-i|--i)
			refIndex=$2
			shift
			;;
		*)
			echo; echo "$1 $2 is not an appropriate argument, please see usage instructions : "; usage; exit
	esac
	shift
	done
fi

#Log File:
echo >> stdout.txt
echo $(date) >> stdout.txt
echo >> stdout.txt

#--------------------
# Downloading Fastq Files:
#--------------------

# TO DO: Test this section once again, made updates (20.04.2021)
echo
echo " 1 - Data Acquisition "

if $download ; then
	echo "Fastq-Dump "$(date) >> stdout.txt
	if [ ! -f $accFileName ]; then
	        echo ; echo "Please provide the name of  a valid SRR file in .csv format"; echo; exit
	fi

	if $pairedEnd ; then
		for i in $(cut -d "," -f1 $accFileName); do
			fast1-dump -I --split-files -O $inoutdir "$i" >> stdout.txt
		done
	#else	#TO DO: Add an else statement for single end data download
	fi

	echo "Fastq Download Complete"; echo "Log: fastq-dump download complete: "$(date) >> stdout.txt
fi

#--------------------
# Data PreProcessing:
#--------------------

echo
echo " 2 - Data PreProcessing: Initial Quality Metrics "

# Fastqc on select samples

#--------------------
# Alignment
#--------------------

echo
echo " 3 - Alignment "

# TO DO: Remove gzip from .fa/.gtf file downloads, reference index wont work otherwise even with --readFilesCommand gunzip -c

# Generating a Genome Reference Index:
if $refGenomeIndex ; then
	if [ ! -d $inoutdir/hg38 ]; then
		mkdir $inoutdir/hg38
                echo "Generated Directory "$inoutdir"/hg38"
	else
		echo $inoutdir"/hg38 exists, continuing."
	fi

	echo ; echo "Generating Genome Index"; echo "Generating Reference Genome Index "$(date) >> stdout.txt
	# Downloading Reference Genome
	if [ ! -f $inoutdir/hg38/hg38.fa.gz ]; then
		wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P $inoutdir/hg38/ >> stdout.txt
	fi
	# Downloading Annotation File
	if [ ! -f $inoutdir/hg38/hg38.ensGene.gtf.gz ]; then
		wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz -P $inoutdir/hg38/ >> stdout.txt
	fi
	# Creating reference index - note that the number of threads can be adjusted accordingly here
	if [ -f $inoutdir/hg38/hg38.fa ] && [ -f $inoutdir/hg38/hg38.ensGene.gtf ]; then
		STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $inoutdir/hg38/ --genomeFastaFiles $inoutdir/hg38/hg38.fa --sjdbGTFfile $inoutdir/hg38/hg38.ensGene.gtf >> stdout.txt # 2>> stderr.txt
	fi

	refIndex=$inoutdir/hg38

	echo ; echo "Complete"; echo "Genome Index Created "$(date) >> stdout.txt

fi

# Alignment with STAR:
echo "Generating list of files for mapping of multiple in one run"$(date) >> stdout.txt

files1=$inoutdir/*1.fastq; set -- $files1
files2=$inoutdir/*2.fastq; set -- $files2
starFiles1=""
starFiles2=""
for i in $files1; do starFiles1="${starFiles1},$i"; done
for i in $files2; do starFiles2="${starFiles2},$i"; done
starFiles1="${starFiles1:1}"; starFiles2="${starFiles2:1}"

echo "Read 1 list: "$starFiles1"\nRead 2 list: "$starFiles2 >> stdout.txt

if [ ! -d $inoutdir/alignments }; then mkdir $inoutdir/alignments; fi

STAR --runThreadN 4 --readFilesCommand gunzip -c --genomeDir $refIndex --readFilesIn $starFiles1 $starFiles2 --outFilterMultimapNmax 1 \
	--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix $inoutdir/alignments 2>> stderr.txt

# Note the fastq files are not in gzip format post accession step - need to zip them prior to this alignment step
# Not finding the enome file /home/knorwood/data/GSE157103/GSE157103/hg38//genomeParameters.txt
# TO DO: Need a --outSAMattrRGline for corresponding read groups? (in the above STAR alignment command)
# TO DO: Dont think its necessary to add --clip5pNbases but need to double check, will need to make this so that it could be user input?
