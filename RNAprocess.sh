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
	echo '-i : path to the reference genome index for alignment'
	echo '-n : true if running alignment, otherwise set to false'
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
		-n|--n)
			align=$2
			shift
			;;
		*)
			echo; echo "$1 $2 is not an appropriate argument, please see usage instructions : "; usage; exit
	esac
	shift
	done
fi

#Log File:
echo > stdout.txt
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

if $pairedEnd ;then
	files1=$inoutdir/*1.fastq; set -- $files1; files2=$inoutdir/*2.fastq; set -- $files2
else
	files=$inoutdir/*.fastq; set -- $files
fi

# FastQC for Paired and Single End Data
#--------------------

# TO DO: Change so that it downloads gzip files

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

	echo ; echo "Generating Genome Index"; echo -e "\nGenerating Reference Genome Index "$(date) >> stdout.txt
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

	echo ; echo "Complete"; echo -e "\nGenome Index Created "$(date) >> stdout.txt

fi

# Alignment with STAR:
#--------------------

# TO DO: add gzip function to STAR, to work with zipped files
#TO DO: add a naming system before here as it should work with both paired/single end data, make sure it aligns with proteomic data
# TO DO: remove the align option - just temporary for testing

if $align ;then
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
		if [ ! -d $inoutdir/alignments ]; then mkdir $inoutdir/alignments; echo "Generated Directory "$inoutdir"/alignments" >> stdout.txt; fi
		
		cat $inoutdir/starfile.txt | while read line; do
			base=$(echo "${line%% *}"); starname=$(basename -s .fastq $base | cut -f1 -d"_")
			echo $starname; echo $starname >> stdout.txt
			STAR --runThreadN 4 --genomeDir $refIndex --readFilesIn $line --outFilterMultimapNmax 1 \
				--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile $refIndex/hg38.ensGene.gtf \
				--clip5pNbases $trim --outFileNamePrefix $inoutdir/alignments/$starname"_" #--genomeLoad LoadAndRemove 2>> stderr.txt
			echo -e "$starname\talignment complete"; echo -e "$line\talignment complete"$(date) >> stdout.txt
		done
		echo "Alignment Complete"; echo "Completed Paired Alignment "$(date) >> stdout.txt
		rm -rf $inoutdir/files1.txt $inoutdir/files2.txt
	fi
fi

# Renaming files by Group 
#--------------------

# TO DO: Rename bam and tab files based on groups, no need to rename fastq files? 
# TO DO: Uncomment the Read Counts and Indexing Sections
# TO DO: How to make the renaming structure more univerisal - have user input what column in SRA metadata table is used for group (3rd arg)

python $SOFTWAREPATH/rename_group.py $inoutdir/alignments/ $metaFile "7" # No need for 3rd arg at this time

# Note the fastq files are not in gzip format post accession step - need to zip them prior to this alignment step (use --readFilesCommand gunzip -c )
# TO DO: Need a --outSAMattrRGline for corresponding read groups? (in the above STAR alignment command)

# Read Counts
#--------------------

# TO DO: Add the gzip function (gzip -d ...) to fastq read counts, this will need to be added to accession as well
# TO DO: Replace " " with "\t" in read_counts.txt file - gsub 

echo " 4 - Calculating Read Statistics"; echo "Read Statistics" >> stdout.txt 

bamFiles=$inoutdir/alignments/*.bam; set -- $bamFiles
tabFiles=$inoutdir/alignments/*ReadsPerGene.out.tab; set -- $tabFiles

if [ ! -d $inoutdir/read_counts ]; then mkdir $inoutdir/read_counts; fi

# Read Counts:
if $pairedEnd ;then
	> $inoutdir/read_counts/fastq1_reads.txt; for line in $files1; do rcount1=$(cat $line | wc -l); echo $(basename $line) $(( rcount1 / 4 )) >> $inoutdir/read_counts/fastq1_reads.txt; done
	> $inoutdir/read_counts/fastq2_reads.txt; for line in $files2; do rcount2=$(cat $line | wc -l); echo $(basename $line) $(( rcount2 / 4 )) >> $inoutdir/read_counts/fastq2_reads.txt; done
fi

# Aligned Read Counts:
> $inoutdir/read_counts/bam_reads.txt; for line in $bamFiles; do bamcount=$(samtools view $line | wc -l); echo $(basename $line) $bamcount >> $inoutdir/read_counts/bam_reads.txt; done

# Exonic Reads Count:
> $inoutdir/read_counts/exonic_reads.txt; for line in $tabFiles; do tabcount=$(awk '{if ($1 ~ "ENSG") {sum += $4}} END {print sum}' $line); echo $(basename $line) $tabcount >> $inoutdir/read_counts/exonic_reads.txt; done
grep -v -e '^$'  $inoutdir/read_counts/exonic_reads.txt

# Read Stats:
paste $inoutdir/read_counts/fastq1_reads.txt $inoutdir/read_counts/fastq2_reads.txt $inoutdir/read_counts/bam_reads.txt $inoutdir/read_counts/exonic_reads.txt > $inoutdir/read_counts/read_counts.txt

#rm -rf $inoutdir/read_counts/fastq1_reads.txt $inoutdir/read_counts/fastq2_reads.txt $inoutdir/read_counts/bam_reads.txt $inoutdir/read_counts/exonic_reads.txt

cat $inoutdir/read_counts/exonic_reads.txt >> stdout.txt; echo -e "\nRead Count Stats - Complete"$(date) >> stdout.txt

# Indexing
#--------------------

echo -e "\nIndexing Bam Files\n" >> stdout.txt
for line in $bamFiles; do samtools index $line; done
echo $inoutdir/alignments/*.bam.bai >> stdout.txt

# Bigwig file creation
#echo "Creation of Bigwig files" >> stdout.txt 
#for line in $bamFiles; do bamCoverage -b $line -o "$line.bw" --normalizeUsingRPKM; done # bamCoverage command not found
echo -e "\nComplete" >> stdout.txt

#--------------------
# Normalization
#--------------------

echo -e "\nNormalization"$(date) >> stdout.txt

# Target Files:
#--------------------

# TO DO: Delete or keep the tabfiles.txt file as it could be used for target_files.txt rather than making a new one

if [ ! -d $inoutdir/normalization ]; then mkdir $inoutdir/normalization; fi; echo -e "\nNormalization Directory Created" >> stdout.txt

> $inoutdir/normalization/target_file.txt
echo -e "files\tgroup\tshort_name" >> $inoutdir/normalization/target_file.txt
for i in $tabFiles; do group=$(basename -s .out $i | cut -f2 -d"_"); srr=$(basename -s .out $i | cut -f1 -d "_"); short_name=$srr"_"$group;
	echo -e $i"\t"$group"\t"$short_name >> $inoutdir/normalization/target_file.txt; done
# echo -e $tabFiles"\t"$(basename -s .out $tabfiles | cut -f2 -d"_")"\t"$(basename -s .fastq $base | cut -f1 -d"_")"_"$(basename -s .fastq $base | cut -f2 -d"_")  >> $inoutdir/normalization/target_files.txt

echo
echo " 5 - Inter and Intra Normalization "

# Normalization
#--------------------

# TO DO: Add a gene annotation file as the fourth argument

baseName=$(basename $inoutdir)

echo -e "\n" >> stdout.txt
RScript $SOFTWAREPATH/normalization.r $inoutdir/alignments $inoutdir/normalization/ $geneLength $baseName >> stdout.txt

