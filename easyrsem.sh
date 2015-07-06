#!/bin/bash

#This script runs RSEM automatically and will merge the resulting files into a table for use in DESeq, it was also generate a TPM file for abundance comparisons.

#set arguments
#samples must be .fastq and prefixed by samplename then _1 or _2 before .fastq to indicate the read pairing

#sort out .fastq and .fq differences

echo -e "\nWelcome to the automated RSEM pipeline\n"
echo -e "\nMake sure you have the following in your current working directory\n\n1. Reads with prefix and _1 or _2 for pairing\n2. If you have it .gtf annotation and .fa genome\n3. If you have it .fa transcriptome\n\nEnter 'no' if you do not have these\n"
echo -e "\nPlease enter the name your bowtie index is to be built by e.g. maize\n"
read name
echo -e "\ndo you have a reference annotation, 'y' or 'n', enter any other letter to skip prepare reference\n"
read anno
echo -e "\nPlease enter whether you have 'paired' or 'single' reads\n"
read pairs
echo -e "\nPlease enter your bowtie path in the format ~/foldername etc, on chris computer, this is, ~/seqsoftware/bowtie-1.0.0/ \n"
read path
echo -e "\nPlease enter the number of replicates\n"
read reps

#RSEM (alignment and calculate expression)

echo -e "\nPlease enter prefix 1\n"
read id1
echo -e "Please enter prefix 2\n"
read id2
echo -e "Please enter prefix 3\n"
read id3
echo -e "Please enter prefix 4\n"
read id4

if [ "$reps" == "3" ]; then
	echo -e "\nPlease enter prefix 5\n"
	read id5
	echo -e "Please enter prefix 6\n"
	read id6
fi

#clean reads

logfile=./log.txt

cp ~/seqsoftware/Trimmomatic/Trimmomatic-0.27/trimmomatic-0.27.jar ./

echo -e "\nenter the phred quality coding, either '64' or '33', (64 for setaria, 33 for chang, usually 33)\n"
read qual
echo -e "\nare you trimming adaptors? enter 'y' if yes or 'n' if no\n" 
read adapt
mkdir trimmed
echo -e "\ndo you need a fastqc report of uncleaned reads? if yes 'y' if no 'n'\n"
read fastq

if [ "$fastqc" == "y" ]; then
	for UNCLEANFILE in ./*.fq; do
		echo -e "\nproducing fastqc reports"
		$HOME/seqsoftware/FastQC/fastqc $UNCLEANFILE
	done
	echo -e  "\ndone with fastqc\n"
fi
#done

echo -e "\nquality trimming - writing log into file logfile.txt - can collect cleaning % from here\n"

exec > $logfile 2>&1

mkdir theunclean

function trim {
	for FILE in ./*_1.fq; do
		java -jar trimmomatic-0.27.jar PE -threads 8 -phred"$qual" $FILE "${FILE/_1.fq/_2.fq}" ./trimmed/"${FILE/_1.fq/_1_paired.fq}" ./trimmed/"${FILE/_1.fq/_1_unpaired.fq}" ./trimmed/"${FILE/_1.fq/_2_paired.fq}" ./trimmed/"${FILE/_1.fq/_2_unpaired.fq}" LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:60
		mv ./trimmed/"${FILE/_1.fq/_1_paired.fq}" ./
		mv ./trimmed/"${FILE/_1.fq/_2_paired.fq}" ./
		mv FILE ./theunclean
		mv "${FILE/_1.fq/_2.fq}" ./theunclean
	done
}

function trimwithadap {
	for FILE in ./*_1.fq; do
		java -jar trimmomatic-0.27.jar PE -threads 8 -phred"$qual" $FILE "${FILE/_1.fq/_2.fq}" ./trimmed/"${FILE/_1.fq/_1_paired.fq}" ./trimmed/"${FILE/_1.fq/_1_unpaired.fq}" ./trimmed/"${FILE/_1.fq/_2_paired.fq}" ./trimmed/"${FILE/_1.fq/_2_unpaired.fq}" ILLUMINACLIP:adaptors.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:10 MINLEN:60
		mv ./trimmed/"${FILE/_1.fq/_1_paired.fq}" ./
		mv ./trimmed/"${FILE/_1.fq/_2_paired.fq}" ./
		mv FILE ./theunclean
		mv "${FILE/_1.fq/_2.fq}" ./theunclean
	done
}

if [ "$adapt" == "y" ]; then
	trimwithadap
fi

if [ "$adapt" == "n" ]; then
	trim
fi

echo -e "\ndone trimming\n"

#include fastqc report of all cleaned reads

for CLEANFILE in ./*.fq; do
	echo -e "\nproducing fastqc reports"
	$HOME/seqsoftware/FastQC/fastqc $CLEANFILE
done

echo -e "\nwe are done \n"

#done

mkdir aligned%

######

if [ "$anno" == "y" ]; then
	rsem-prepare-reference --bowtie-path ~/seqsoftware/bowtie-1.0.0/ --gtf *.gtf *.fa $name
fi

if [ "$anno" == "n" ]; then
	rsem-prepare-reference --bowtie-path ~/seqsoftware/bowtie-1.0.0/ *.fa $name
fi

#run rsem

if [ "$pairs" == "paired" ]; then
	echo -e "\nRunning RSEM\n"
	for FILE in ./*_1.fastq; do
		rsem-calculate-expression -p 8 --bowtie-path $path --paired-end $FILE "${FILE/_1.fastq/_2.fastq}" $name "${FILE/_1.fastq/}"
	done
fi

if [ "$pairs" == "single" ]; then
	echo -e "\nRunning RSEM\n"
	for FILE in ./*.fq; do
		rsem-calculate-expression -p 8 --bowtie-path $path $FILE $name "${FILE/.fq/}"
	done
fi

echo -e "\nwe are done\n"

#merging data outputs for DESeq (remove top line, remove no counts and add one pseudocount - to prevent dividing by zero)

echo -e "\nremoving line1\n"

if [ "$reps" == "2" ]; then
	awk 'FNR>1' "$id1".genes.results > "$id1".genes.results2.txt
	awk 'FNR>1' "$id2".genes.results > "$id2".genes.results2.txt
	awk 'FNR>1' "$id3".genes.results > "$id3".genes.results2.txt
	awk 'FNR>1' "$id4".genes.results > "$id4".genes.results2.txt
fi

if [ "$reps" == "3" ]; then
	awk 'FNR>1' "$id1".genes.results > "$id1".genes.results2.txt
	awk 'FNR>1' "$id2".genes.results > "$id2".genes.results2.txt
	awk 'FNR>1' "$id3".genes.results > "$id3".genes.results2.txt
	awk 'FNR>1' "$id4".genes.results > "$id4".genes.results2.txt
	awk 'FNR>1' "$id5".genes.results > "$id5".genes.results2.txt
	awk 'FNR>1' "$id6".genes.results > "$id6".genes.results2.txt
fi

if [ "$reps" == "2" ]; then
	echo -e "\nmerging counts \n"
	paste ./"$id1".genes.results2.txt ./"$id2".genes.results2.txt ./"$id3".genes.results2.txt ./"$id4".genes.results2.txt | awk -F "\t" '{print $1"\t"$5"\t"$12"\t"$19"\t"$26 }' > RSEMcounts.txt
	echo -e "\ndeleting null features and adding pseudocount=1 \n"
	awk -F "\t" '{if ($2>0 || $3>0 || $4>0 || $5>0) print $1"\t"$2+1"\t"$3+1"\t"$4+1"\t"$5+1}' RSEMcounts.txt > RSEMcountsnozeros.txt
fi

if [ "$reps" == "3" ]; then
	echo -e "\nmerging counts \n"	
	paste ./"$id1".genes.results2.txt ./"$id2".genes.results2.txt ./"$id3".genes.results2.txt ./"$id4".genes.results2.txt ./"$id5".genes.results2.txt ./"$id6".genes.results2.txt | awk -F "\t" '{print $1"\t"$5"\t"$12"\t"$19"\t"$26"\t"$33"\t"$40 }' > RSEMcounts.txt					 										  
	echo -e "\ndeleting null features and adding pseudocount=1 \n"
	awk -F "\t" '{if ($2>0 || $3>0 || $4>0 || $5>0 || $6>0 || $7>0) print $1"\t"$2+1"\t"$3+1"\t"$4+1"\t"$5+1"\t"$6+1"\t"$7+1}' RSEMcounts.txt > RSEMcountsnozeros.txt
fi

#TPMs

if [ "$reps" == "3" ]; then
	echo -e "\nmerging TPM valyes \n"
	paste ./"$id1".genes.results2.txt ./"$id2".genes.results2.txt ./"$id3".genes.results2.txt ./"$id4".genes.results2.txt ./"$id5".genes.results2.txt ./"$id6".genes.results2.txt | awk -F "\t" '{print $1"\t"$6"\t"$13"\t"$20"\t"$27"\t"$34"\t"$41 }' > TPMvalues.tsv
fi

if [ "$reps" == "2" ]; then
	echo -e "\nmerging TPM valyes \n"	
	paste ./"$id1".genes.results2.txt ./"$id2".genes.results2.txt ./"$id3".genes.results2.txt ./"$id4".genes.results2.txt | awk -F "\t" '{print $1"\t"$6"\t"$13"\t"$20"\t"$27 }' > TPMvalues.tsv
fi

#header

echo -e "\nadding a header \n"

if [ "$reps" == "2" ]; then
	echo -e "gene""\t$id1""\t$id2""\t$id3""\t$id4" >> header.txt
fi

if [ "$reps" == "3" ]; then
	echo -e "gene""\t$id1""\t$id2""\t$id3""\t$id4""\t$id5""\t$id6" >> header.txt
fi

cat header.txt RSEMcountsnozeros.txt > RSEMcountsb4rounding.txt
cat header.txt TPMvalues.tsv > TPMs.tsv

#code to round to whole numbers

echo -e "\nrounding to whole numbers\n"

if [ "$reps" == "2" ]; then
	cat RSEMcountsb4rounding.txt | awk -F "\t" '{printf $1"\t" "%.f\t%.f\t%.f\t%.f\n", $2, $3, $4, $5}' > RSEMcountsprefinal.txt
fi

if [ "$reps" == "3" ]; then
	cat RSEMcountsb4rounding.txt | awk -F "\t" '{printf $1"\t" "%.f\t%.f\t%.f\t%.f\t%.f\t%.f\n", $2, $3, $4, $5, $6, $7}' > RSEMcountsprefinal.txt
fi

awk 'FNR>1' RSEMcountsprefinal.txt > RSEMcountsprefinal2.txt
cat header.txt RSEMcountsprefinal2.txt > RSEMcountsfinal.txt
mkdir final
mv TPMs.tsv ./final
cp RSEMcountsfinal.txt ./final
cp ./final/RSEMcountsfinal.txt ./final/data.tsv

#counting aligned

mkdir aligned%

for file in ./*.transcript.bam; do
	aligned=$(samtools view -F 0x4 $file | cut -f 1 | sort | uniq | wc -l)
	lines=$(cat "${file/.transcript.bam/_1.fastq}" | wc -l)
	reads=$(($lines/4))
	echo "file: " $file >> aligned%/%.txt
	echo "Reads in:	" $reads >> aligned%/%.txt	# divide by 4 because there are 4 lines per read
	echo "Reads aligned: " $aligned >> aligned%/%.txt
	echo "Percentage reads aligned: " $((($aligned*100/$reads)))"%" >> aligned%/%.txt
done
		
#####

echo -e "\ndone\n"

#done








