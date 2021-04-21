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
			. /$confFile
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
		*)
			echo; echo "$1 $2 is not an appropriate argument, please see usage instructions : "; usage; exit
	esac
	shift
	done
fi

#Log File:
echo >> stdout.txt
echo $(date) >> stdout.txt
echo >> std.out

#--------------------
# Downloading Fastq Files:
#--------------------

echo
echo " 1 - Data Acquisition "

if [ $download==true ]; then
	echo "Fastq-Dump "$(date) >> stdout.txt
	if [ ! -f $accFileName ]; then
	        echo ; echo "Please provide the name of  a valid SRR file in .csv format"; echo; exit
	fi

	if [ $pairedEnd==true ]; then
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

if [ $refGenomeIndex==true ]; then # TO DO: Make it so that the reference gene index does not just get saved in the same file as everything else
	echo "Generating Reference Genome Index "$(date) >> stdout.txt
	# Downloading Reference Genome
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P $inoutdir; gzip -d $inoutdir/hg38.fa.gz
	# Downloading Annotation File
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz -P $inoutdir; gzip -d $inoutdir/hg38.ensGene.gtf.gz
	# Creating reference index - note that the number of threads can be adjusted accordingly here
	STAR --runMode genomeGenerate --runThreadN 1 --genomeDir $inoutdir --readFilesCommand gunzip -c --genomeFastaFiles $inoutdir/hg38.fa.gz --sjdbGTFfile $inoutdir/hg38.ensGene.gtf.gz

	echo "Genome Index Created "$(date) >> stdout.txt
fi

# Alignment with STAR:
#for [ i".fa" in $inoutdir ] # or something like this to only take fast files, may be easier to make a new fast file directory earlier

