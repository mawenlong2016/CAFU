
#!/bin/bash

## $1:$output
## $2:$reference_file.files_path/$reference_name
## $3:/home/malab10/software/galaxy/database/files/$ngvec_name.ngvec
my=/home/chenzhuod/galaxy/database/files/expression/DE_analysis
path1=/home/chenzhuod/galaxy/database/files/expression/expression_evidence/

##########################
echo -n "" > ${my}results
##########################

IFS=' '

cat /home/chenzhuod/galaxy/database/files/expression/DE_analysis/RNA-Seq_info | while read line
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
			
	arrLeft=""
	for str in ${lefts[@]}
	do
	
	##########################
		echo -n " ${path1}${str}_exp.isoforms.results" >> ${my}sample
	##########################
	
		arrLeft="${arrLeft} ${path1}${str}_exp.isoforms.results"
	done
	arrRight=""
	for str in ${rights[@]}
	do
	
	##########################
		echo -n " ${path1}${str}_exp.isoforms.results" >> ${my}sample
	##########################
	
		arrRight="${arrRight} ${path1}${str}_exp.isoforms.results"
	done	
	
	
	IFS=" "
	lefts_len=${#lefts[@]}
	rights_len=${#rights[@]}
	
	
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

awk '{ print $1 }' $1 | sort | uniq > ${4}

########################################
/home/chenzhuod/galaxy/tools/CAFU/RSEM-1.3.1/obtain_transcripts_expression_matrix `cat ${my}sample` > ${my}results 
cat ${my}results > ${5}	
####################################

