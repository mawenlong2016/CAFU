#!/bin/bash

# 2018-06-19
# ZHAO Xu-yang
# YangLing

# init start ...
#
# ${my}:		the root path,create by XML
# $1:			output_R1 to XML
# $2:			output_R2 to XML
# $3:			minIntron
# $4:			maxIntron
# ${resultR1}:	to $1
# ${resultR2}:	to $2
# ...
# get prarm start ...

minIntron=$3
maxIntron=$4
threads=$5

echo "minIntron=${3}"
echo "maxIntron=${4}"

# get prarm end ...
my=/home/chenzhuod/galaxy/database/files/HISAT/
mkdir -p ${my}hisat_idx/
resultR1=${my}resultR1
resultR2=${my}resultR2
echo -n "" > ${resultR1}
echo -n "" > ${resultR2}
# init end ...



# create index start ...
# XML upload collection to "/home/chenzhuod/galaxy/database/files/HISAT/collection/"
# 
# file="thisfile.txt"
# filename=${file%.*}
# extension=${file##*.}
# 
#...
cpath=${my}genome_collection/
for file in `ls ${cpath}`
do
	filename=${file%.*}
	mkdir -p /home/chenzhuod/galaxy/database/files/HISAT/hisat_idx/${filename}/
	/home/chenzhuod/galaxy/tools/CAFU/hisat2-2.1.0/hisat2-build -p $threads  ${cpath}${file} /home/chenzhuod/galaxy/database/files/HISAT/hisat_idx/${filename}/${filename}_idx
done
# create index end ...



# hisat start ...
#
# ${my}table:	XML upload
# ${c0}:		column 1 in every row
# SRRxxxxxx: 	XML upload to "/home/chenzhuod/galaxy/database/files/HISAT/"
# ${c1}:		column 2 in every row
# ${con}:		all indexID in c1
# ${faR1}:		temp of R1 for every element in every row
# ${faR2}:		temp of R1 for every element in every row
# ${R1}:		temp of R1 in every row
# ${R2}:		temp of R2 in every row
# ...
cat ${my}mapping-operation_info | while read line
do
	IFS=";"
	arr=($line)
	c0=${arr[0]}
	c1=${arr[1]}
	
	# R1,R2 save as temp start ...
	cat ${my}Trimmed_${c0}_1.fastq > ${my}${c0}_R1_temp  
	cat ${my}Trimmed_${c0}_2.fastq > ${my}${c0}_R2_temp
	R1_temp=${my}${c0}_R1_temp
	R2_temp=${my}${c0}_R2_temp
	# R1,R2 save as temp end ...
	
	faR1=${my}R1
	faR2=${my}R2

	IFS=","
	con=($c1)
	for s in ${con[@]}
	do
	
		faID=${my}hisat_idx/${s}/${s}_idx
		temp_output_sam=${my}${c0}_output.sam
		temp_output_bam=${my}${c0}_output.bam
		temp_unmap_bam=${my}${c0}_unmap.bam
		unamp_bam=${my}${c0}_unamp_sort.bam

	
		/home/chenzhuod/galaxy/tools/CAFU/hisat2-2.1.0/hisat2 --min-intron $minIntron --max-intron $maxIntron -p $threads -x ${faID} -1 ${R1_temp} -2 ${R2_temp} -S ${temp_output_sam}
		
		/home/chenzhuod/galaxy/tools/CAFU/samtools-1.8/samtools view -@ $threads -bS ${temp_output_sam} > ${temp_output_bam}

		/home/chenzhuod/galaxy/tools/CAFU/samtools-1.8/samtools view -@ $threads -b -f 4 ${temp_output_bam} > ${temp_unmap_bam}
		
		/home/chenzhuod/galaxy/tools/CAFU/samtools-1.8/samtools sort -@ $threads ${temp_unmap_bam} > ${unamp_bam}

		/home/chenzhuod/galaxy/tools/CAFU/bedtools2/bin/bedtools bamtofastq -i ${unamp_bam} -fq ${faR1} -fq2 ${faR2}

		cat ${faR1} > ${R1_temp}
		cat ${faR2} > ${R2_temp}
	
	done
	
	cat ${faR1} >> ${resultR1}
	cat ${faR2} >> ${resultR2}
done
# hisat start ...

cat ${resultR1} > $1
cat ${resultR2} > $2
