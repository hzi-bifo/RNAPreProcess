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
	echo '-t : trim number for alignment'
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

		-t|--t)
			trim=$2
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
echo "Initializing Pipeline... "$(date) >> stdout.txt
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
			fastq-dump -I --split-files -O $inoutdir "$i" >> stdout.txt
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

#Setting up list of data files:
#--------------------

# TO DO: Clean this up - want to create list of files here, will later be used with STAR but need to do differently for both paired and single end data
if $pairedEnd ;then
	files1=$inoutdir/*1.fastq; set -- $files1; files2=$inoutdir/*2.fastq; set -- $files2
else
	files=$inoutdir/*.fastq; set -- $files
fi

# FastQC for Paired and Single End Data
#--------------------

if $fastqcAnalysis ;then
	if [ ! -d $inoutdir/fastqcAnalysis ]; then mkdir $inoutdir/fastqcAnalysis; fi
	if $pairedEnd ;then
		echo "Running fastqc on paired end data: "$(date) >> stdout.txt
		if [ -z "$files2" ]; then
			echo "No _2.fastq files detected running fastqc on only one file - double check input if unexpected" >> stdout.txt
			fastqcFile=(${files1[0]})
			fastqc -o $inoutdir/fastqcAnalysis/ $fastqcFile >> stdout.txt
		else
			fastqcFile1=(${files1[0]}); fastqcFile2=(${files2[0]})
			fastqc -o $inoutdir/fastqcAnalysis/ $fastqcFile1 $fastqcFile2 >> stdout.txt
		fi
	else
		echo "Fastqc only available for paired end data right now."
	fi
echo "\nFastQC Analysis Complete - please check output, update config file and rerun with appropriate trim"
exit
fi

# Adapter Trimming Check
if [ -z "$trim" ] && [[ $fastqcAnalysis == "false" ]]; then echo "Please add a trim quantity to the config file!"; exit; fi

#--------------------
# Alignment
#--------------------

echo
echo " 3 - Alignment "

# Generating a Genome Reference Index: So have precreated and then add option to path
#--------------------
if $refGenomeIndex ; then
	if [ ! -d $inoutdir/hg38 ]; then mkdir $inoutdir/hg38; echo "Generated Directory "$inoutdir"/hg38"; fi

	echo ; echo "Generating Genome Index"; echo "Generating Reference Genome Index "$(date) >> stdout.txt
	# Downloading Reference Genome
	if [ ! -f $inoutdir/hg38/hg38.fa ]; then
		wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P $inoutdir/hg38/
	       	gzip -d hg38.fa.gz >> stdout.txt 
	fi
	# Downloading Annotation File
	if [ ! -f $inoutdir/hg38/hg38.ensGene.gtf ]; then
		wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz -P $inoutdir/hg38/
	       	gzip -d hg38.ensGene.gtf.gz >> stdout.txt
	fi
	# Creating reference index - note that the number of threads can be adjusted accordingly here
	if [ -f $inoutdir/hg38/hg38.fa ] && [ -f $inoutdir/hg38/hg38.ensGene.gtf ]; then
		STAR --runMode genomeGenerate --runThreadN 4 --genomeDir $inoutdir/hg38/ --genomeFastaFiles $inoutdir/hg38/hg38.fa --sjdbGTFfile $inoutdir/hg38/hg38.ensGene.gtf #>> stdout.txt # 1>> stderr.txt
	fi

	refIndex=$inoutdir/hg38

	echo ; echo "Complete"; echo "Genome Index Created "$(date) >> stdout.txt

fi

# Alignment with STAR:
#--------------------

#TO DO: add a naming system before here as it should work with both paired/single end data, make sure it aligns with proteomic data
# TO DO: 30.04.2021 - need to merge file1 with file2 in tab deliminated manner for alignment

if $pairedEnd ;then

	# Creating File List
	echo "Generating list of files for mapping of multiple in one run "$(date) >> stdout.txt

	> $inoutdir/files1.txt; > $inoutdir/files2.txt
	for i in $files1; do echo "$i" >> $inoutdir/files1.txt; done; for i in $files2; do echo "$i" >> $inoutdir/files2.txt; done
	paste -d " " $inoutdir/files1.txt $inoutdir/files2.txt > $inoutdir/starfile.txt
	echo "Complete" >> stdout.txt

	# Alignment
	echo "STAR Alignment" >> stdout.txt
	echo $refIndex
	cat $inoutdir/starfile.txt | while read line; do
		base=$(echo "${line%% *}"); starname=$(basename -s .fastq $base | cut -f1 -d"_")
		echo $starname; echo $starname >> stdout.txt
		STAR --runThreadN 4 --genomeDir $refIndex --readFilesIn $line --outFilterMultimapNmax 1 \
			--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile $refIndex/hg38.ensGene.gtf \
			--clip5pNbases 13 --outFileNamePrefix $inoutdir/alignments/$starname #--genomeLoad LoadAndRemove 2>> stderr.txt
		echo -e "$starname\talignment complete"; echo -e "$line\talignment complete"$(date) >> stdout.txt
	done
	echo "Alignment Complete"; echo "Completed Paired Alignment "$(date) >> stdout.txt
fi
# Note the fastq files are not in gzip format post accession step - need to zip them prior to this alignment step (use --readFilesCommand gunzip -c )
# TO DO: Need a --outSAMattrRGline for corresponding read groups? (in the above STAR alignment command)

# Read Counts
#--------------------

# TO DO: Add read counts, print to stdout.txt as well. 

#--------------------
# Normalization
#--------------------

echo
echo " 4 - Inter and Intra Normalization "

# Pull in R script here
