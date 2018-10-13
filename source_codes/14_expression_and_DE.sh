#!/bin/bash

## $1:$output
## $2:$reference_file.files_path/$reference_name
## $3:$window_size $4:$required_quality $5:$minlen $6:$leading $7:$trailing $8:$phred
## $9:/home/malab10/software/galaxy/database/files/$ngvec_name.ngvec
my=/home/chenzhuod/galaxy/database/files/Expression/DE_analysis/

##########################
echo -n "" > ${my}results
##########################

IFS=' '

cat /home/chenzhuod/galaxy/database/files/Expression/DE_analysis/RNA-Seq_info | while read line
do
	IFS=";"
	arr=($line)
	c0=${arr[0]}
	c1=${arr[1]}
	c2=${arr[2]}
	c3=${arr[3]}
	
	IFS=","
	lefts=($c1)
	rights=($c2)
	
	g="${g} ${my}${s}_exp.isoforms.results"
		
	arrLeft=""
	for str in ${lefts[@]}
	do
	
	##########################
		echo -n " ${my}${str}_exp.isoforms.results" >> ${my}sample
	##########################
	
		arrLeft="${arrLeft} ${my}${str}_exp.isoforms.results"
	done
	arrRight=""
	for str in ${rights[@]}
	do
	
	##########################
		echo -n " ${my}${str}_exp.isoforms.results" >> ${my}sample
	##########################
	
		arrRight="${arrRight} ${my}${str}_exp.isoforms.results"
	done	
	
	
	IFS=" "
	lefts_len=${#lefts[@]}
	rights_len=${#rights[@]}
	
	
	if [[ "${c3}" == "SE" ]]
	then
		for s in ${lefts[@]}
		do
				java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred$8 ${my}${s}.fastq ${my}${s}_trim.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
				/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 10 --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s}_trim.fastq $2 ${my}${s}_exp	
			
		done
		for s in ${rights[@]}
		do
				java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar SE -phred$8 ${my}${s}.fastq ${my}${s}_trim.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
				/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 10 --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s}_trim.fastq $2 ${my}${s}_exp	

		done		
	fi
	
	if [[ "${c3}" == "PE" ]]
	then
		for s in ${lefts[@]}
		do
			s1=${s}_1.fastq
			s2=${s}_2.fastq
			s1_trim=${s}_1_trim.fastq
			s2_trim=${s}_2_trim.fastq
			java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred$8 ${my}${s1} ${my}${s2} ${my}${s1_trim} ${my}${s}_1_unpaired.fastq ${my}${s2_trim} ${my}${s}_2_unpaired.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
			/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 10--paired-end --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s1_trim} ${my}${s2_trim} $2 ${my}${s}_exp	  	
		done
		for s in ${rights[@]}
		do
			s1=${s}_1.fastq
			s2=${s}_2.fastq
			s1_trim=${s}_1_trim.fastq
			s2_trim=${s}_2_trim.fastq
			java -jar /home/chenzhuod/galaxy/tools/CAFU/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred$8 ${my}${s1} ${my}${s2} ${my}${s1_trim} ${my}${s}_1_unpaired.fastq ${my}${s2_trim} ${my}${s}_2_unpaired.fastq SLIDINGWINDOW:$3:$4 MINLEN:$5 LEADING:$6 TRAILING:$7
			/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-calculate-expression -p 10 --paired-end --bowtie2 --bowtie2-path /home/chenzhuod/galaxy/tools/CAFU/bowtie2-2.3.4.3-linux-x86_64 ${my}${s1_trim} ${my}${s2_trim} $2 ${my}${s}_exp	  	
		done	
		
	fi
	arrLeftRight="${arrLeft} ${arrRight}"
	/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-generate-data-matrix ${arrLeftRight} > ${my}${c0}Mat
	/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-run-ebseq --ngvector $9 ${my}${c0}Mat $lefts_len,$rights_len ${my}${c0}_Results
	/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/rsem-control-fdr ${my}${c0}_Results 0.05 ${my}${c0}_DEResults
	cat ${my}${c0}_DEResults > ${my}${c0}_DEResults_tmp
	
	sed -i 1d ${my}${c0}_DEResults_tmp
	sed -i 's/"//g' ${my}${c0}_DEResults_tmp
	cat ${my}${c0}_DEResults_tmp | while read line
	do
		IFS="	"
		arr=($line)
		r1=${arr[0]}
		r4=${arr[3]}
		[ $(echo "${r4}>=2" | bc) = 1 ] && {
				echo "${r1}	${r4}	${c0}" >> $1
		}
		[ $(echo "${r4}<=0.5" | bc) = 1 ] && {
				echo "${r1}	${r4}	${c0}" >> $1
		}
		
	done

done

awk '{ print $1 }' $1 | sort | uniq > ${10}

########################################
/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/obtain_transcripts_expression_matrix `cat ${my}sample` > ${my}results 
cat ${my}results >> ${11}

####################################

