#!/bin/bash
## condition_specific.sh

my=/home/chenzhuod/galaxy/database/files/expression/condition_specificity/

g=""
cat /home/chenzhuod/galaxy/database/files/expression/condition_specificity/RNA-Seq_info | while read line
do
	IFS=";"
	arr=($line)
	c1=${arr[1]}
	c2=${arr[2]}
	IFS=","
	con=($c1)
	
	if [[ "${c2}" == "SE" ]]
	then
		for s in ${con[@]}
		do
			g="${g} ${my}${s}_exp.isoforms.results"
			echo -n ${g} > ${my}sample
			
			java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred$8 ${my}${s}.fastq ${my}${s}_trim.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
			/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 2 --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s}_trim.fastq $2 ${my}${s}_exp	
			
		done	
	fi
	if [[ "${c2}" == "PE" ]]
	then
		for s in ${con[@]}
		do
			g="${g} ${my}${s}_exp.isoforms.results"
			echo -n ${g} > ${my}sample

			s1=${s}_1.fastq
			s2=${s}_2.fastq
			s1_trim=${s}_1_trim.fastq
			s2_trim=${s}_2_trim.fastq
			java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred$8 ${my}${s1} ${my}${s2} ${my}${s1_trim} ${my}${s}_1_unpaired.fastq ${my}${s2_trim} ${my}${s}_2_unpaired.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
			/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 2 --paired-end --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s1_trim} ${my}${s2_trim} $2 ${my}${s}_exp	  	
		done	
	fi	
done

/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/obtain_transcripts_expression_matrix `cat ${my}sample` > ${my}results
cat ${my}results >> $1
