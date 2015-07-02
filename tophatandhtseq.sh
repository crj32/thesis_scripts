#!/bin/bash

#run in directory with .fastq, .gtf, and index (bowtieindex must be under prefix genome)

echo -e "\nwelcome to the tophat2 and htseq alignment and counting pipeline\n"

#Enter option for single end or paired reads

#For Tophat alignment

echo -e "\nis this program for single or paired end files, enter 'single' or 'paired'\n"
read pairs
echo -e "\nis your file a .gtf or .gff, enter 'gff' or 'gtf'\n"
read ext
echo -e "\ndefault ID = 50 and SD = 20, for setaria use ID = 200 and SD = 50\n"
echo -e "\nplease enter inner distance\n"
read inner
echo -e "\nplease enter standard deviation\n"
read sd
echo -e "\nplease enter prefix of genome build\n"
read build

#For HT-SEQ counting

echo -e "\nare you working with a .gtf or .gff file? Please enter 'gtf' or 'gff'\n"
read ext
echo -e "\ndo you need the -i Parent extension for htseq count, please enter '-i Parent' if yes or just hit return\n"
read parent
echo "enter number of replicates: "
read replicates
echo -e "enter ID of sample1 (i.e. prefix to bam file)"
read id1
echo "enter ID of sample2"
read id2
echo "enter ID of sample3"
read id3
echo "enter ID of sample4"
read id4

if [ "$replicates" == "3" ]; then 
	echo "enter ID of sample5"
	read id5
	echo "enter ID of sample6"
	read id6
fi

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

#For singles 

if [ "$pairs" == "single" ]; then
	for FILE in ./*.fastq; do
		tophat -p 8 -G ./*."$ext" -o ./"${FILE/.fastq/}" --no-novel-juncs $build $FILE
		mv ./"${FILE/.fastq/}"/accepted_hits.bam ./"${FILE/.fastq/}"/"${FILE/.fastq/.bam}" 
		#move all the .bam files back into the current directory
		cp ./"${FILE/.fastq/}"/"${FILE/.fastq/.bam}" ./
		#continue
		readsin=$(cat $FILE | wc -l)
       	readsaligned=$(samtools view ./"${FILE/.fastq/}"/"${FILE/.fastq/.bam}" | cut -f 1 | sort | uniq | wc -l)
		mkdir ./"${FILE/.fastq/}"_aligned%/
		echo -e "${FILE/.fastq/}" >> ./"${FILE/.fastq/}"_aligned%/%.txt
        echo "Total number of reads: " $(($readsin/4)) >> ./"${FILE/.fastq/}"_aligned%/%.txt
        echo "Reads aligned: " $readsaligned >> ./"${FILE/.fastq/}"_aligned%/%.txt
        echo "Percent reads aligned: " $(($readsaligned*400/$readsin)) >> ./"${FILE/.fastq/}"_aligned%/%.txt
        done
fi

#For pairs

if [ "$pairs" == "paired" ]; then
	for FILE2 in ./*_1.fastq; do
		tophat -p 8 -G ./*."$ext" -o ./"${FILE2/_1.fastq/}" --no-novel-juncs -r $inner --mate-std-dev $sd $build $FILE2 "${FILE2/_1.fastq/_2.fastq}"
		mv ./"${FILE2/_1.fastq/}"/accepted_hits.bam ./"${FILE2/_1.fastq/}"/"${FILE2/_1.fastq/.bam}"
		#move all the .bam files back into the current directory
		cp ./"${FILE2/_1.fastq/}"/"${FILE2/_1.fastq/.bam}" ./
		#continue
		echo -e "\ncounting reads in fastq file\n"
		readsin=$(wc -l < "$FILE2")
		readsinover4=$(($readsin/4))
		echo -e "\ncounting aligned reads\n"
		readsaligned=$(samtools view ./"${FILE2/_1.fastq/}"/"${FILE2/_1.fastq/.bam}" | cut -f 1 | sort | uniq | wc -l) # lifted from website
		echo -e "\ncalculating percentage aligned reads and outputing to a text file\n"
		echo -e "${FILE2/_1.fastq/}" >> aligned%/%.txt
		echo "Reads in:	" $readsinover4	>> aligned%/%.txt	# divide by 4 because there are 4 lines per read
		echo "Reads aligned: " $readsaligned >> aligned%/%.txt
		echo "Percentage reads aligned: " $(($readsaligned*100/$readsinover4))"%" >> aligned%/%.txt
		echo -e "\ndone\n"		
	done
fi

echo -e "\nwe are done with tophat2\n"

#sort the bam files based on read name

echo -e "\nsorting \n"

for file in ./*.bam ; do
	echo $file
	samtools sort -n $file $file
done

#rename the files to remove extra .bam

rename s/.bam.bam/_sorted.bam/ *.bam.bam

#convert multiple bam files to sam files

echo -e "\nconverting\n"

for file2 in ./*_sorted.bam ; do
	echo $file2
	samtools view -h $file2 > ${file2/.bam/.sam}
done

#run htseq on the sorted sam file

echo -e "\nrunning htseq on the sam file\n"

for file3 in ./*.sam ; do
	echo $file3
	htseq-count -s no $parent $file3 *."$ext" >> ${file3/.sam/.counts.txt}
done

echo -e "\n$(date) -- finished --\n"

#merge count tables into one file - seperates based on a tab - MODIFY THIS FOR THE SAMPLE NUMBER AND NAME EVERYTIME

if [ "$replicates" == "2" ]; then
	echo -e "\nmerging counts \n"
	paste ./"$id1"_sorted.counts.txt ./"$id2"_sorted.counts.txt ./"$id3"_sorted.counts.txt ./"$id4"_sorted.counts.txt | awk -F "\t" '{print $1"\t"$2"\t"$4"\t"$6"\t"$8 }' > All_htseqCounts.txt
	echo -e "\ndeleting null features and adding pseudocount=1 \n"
	awk -F "\t" '{if ($2>0 || $3>0 || $4>0 || $5>0) print $1"\t"$2+1"\t"$3+1"\t"$4+1"\t"$5+1}' All_htseqCounts.txt | grep -v no_feature | grep -v ambiguous | grep -v too_low_aQual | grep -v 		not_aligned | grep -v alignment_not_unique > All_Counts_nozero_1pseudocount.txt
fi

if [ "$replicates" == "3" ]; then
	echo -e "\nmerging counts \n"	
	paste ./"$id1"_sorted.counts.txt ./"$id2"_sorted.counts.txt ./"$id3"_sorted.counts.txt ./"$id4"_sorted.counts.txt ./"$id5"_sorted.counts.txt ./"$id6"_sorted.counts.txt | awk -F "\t" '{print $1"\t"$2"\t"$4"\t"$6"\t"$8"\t"$10"\t"$12 }' > All_htseqCounts.txt						 										  
	echo -e "\ndeleting null features and adding pseudocount=1 \n"
	awk -F "\t" '{if ($2>0 || $3>0 || $4>0 || $5>0 || $6>0 || $7>0) print $1"\t"$2+1"\t"$3+1"\t"$4+1"\t"$5+1"\t"$6+1"\t"$7+1}' All_htseqCounts.txt | grep -v no_feature | grep -v ambiguous | grep -v too_low_aQual | grep -v not_aligned | grep -v alignment_not_unique > All_Counts_nozero_1pseudocount.txt
fi

#add a header to the count file

echo -e "\nadding a header \n"

if [ "$replicates" == "2" ]; then
	echo -e "gene""\t$id1""\t$id2""\t$id3""\t$id4" >> header.txt
fi

if [ "$replicates" == "3" ]; then
	echo -e "gene""\t$id1""\t$id2""\t$id3""\t$id4""\t$id5""\t$id6" >> header.txt
fi

cat header.txt All_Counts_nozero_1pseudocount.txt > All_Counts_no0_1pc_Header.txt

echo -e "\nputting in new directory\n"

mkdir final
cp All_Counts_no0_1pc_Header.txt ./final
mv ./final/All_Counts_no0_1pc_Header.txt ./final/data.tsv

echo -e "\ncleaning up\n"

########## enter code

rm *_sorted.sam
rm *_sorted.bam

echo -e "\nfile is now ready for DE analysis \n"

#Done










