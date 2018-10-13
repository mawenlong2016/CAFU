ls /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/ | grep '_exp.transcript.bam.bed' | sed 's/_exp.transcript.bam.bed//g' > /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/RNA-Seq_ID
mkdir -p /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/bam_bed
mv /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/*_exp.transcript.bam.bed /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/bam_bed

ls /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/transBed | grep 'bed' | sed -r 's/\.bed//g' > /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/transcript_ID

for i in $(cat /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/transcript_ID)
do
mkdir -p /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/coverage/${i}
for j in $(cat /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/RNA-Seq_ID)
do
/home/chenzhuod/galaxy/tools/CAFU/bedtools2/bin/bedtools coverage -a /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/transBed/${i}.bed -b /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/bam_bed/${j}_exp.transcript.bam.bed -d > /home/chenzhuod/galaxy/database/files/Expression/expression_evidence/coverage/${i}/${j}_coverage
done
done
