#!/bin/bash
##tissue_specific.sh

my=/home/chenzhuod/galaxy/database/files/Expression/expression_evidence/
g=""

cat /home/chenzhuod/galaxy/database/files/HISAT/mapping-operation_info | while read line
do
	IFS=";"
	arr=($line)
	c0=${arr[0]}
	s1_trim=Trimmed_${c0}_1.fastq
	s2_trim=Trimmed_${c0}_2.fastq
	/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p $3 --paired-end --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 /home/chenzhuod/galaxy/database/files/HISAT/${s1_trim} /home/chenzhuod/galaxy/database/files/HISAT/${s2_trim} $2 ${my}${c0}_exp		
	/home/chenzhuod/galaxy/tools/CAFU/bedtools2/bin/bamToBed -i ${my}${c0}_exp.transcript.bam > ${my}${c0}_exp.transcript.bam.bed

	g="${g} ${my}${c0}_exp.isoforms.results"
	echo -n ${g} > ${my}sample
	
done
	
/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/obtain_transcripts_expression_matrix `cat ${my}sample` > ${my}results
cat ${my}results > $1
