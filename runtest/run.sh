#!/bin/bash

gatk=gatk-package-distribution-3.6.jar
#gatk=huaweiyun_gatk.jar
lib_path=JNILib

db=/home/hust/storage/data/known_database
reference_file=/mnt/sdc/reference_sequence/hs37d5.fasta
dbsnp_del100=$db/dbSNP/dbsnp_138.b37.vcf.gz
mills_1kg=$db/1000G_gold_standard/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
cosmic=$db/COSMIC/b37_cosmic_v73_061615.vcf.gz
targetBed=./target.bed
#targetBed=./chr22.bed
#targetBed=./target_chr1.bed
#targetBed=./test_for_chenzhuo.bed
#targetBed=./test.bed
out_dir=/mnt/sdc/gene_plus_realign_output
tmp_dir=./tmp


####################################################################################################
##############################################   Mutect2  ##########################################
####################################################################################################
#echo  "#######################################   Mutect2  ##########################################"


time java \
-d64 -server -XX:+UseG1GC -XX:ParallelGCThreads=2 -Xms3g -Xmx30g -Djava.io.tmpdir=$tmp_dir -Djava.library.path=$lib_path -jar $gatk -T FastMuTect2 -L $targetBed -R $reference_file -I:tumor $out_dir/case_sort_markdup_realign.bam -I:normal $out_dir/normal_sort_markdup_realign.bam --dbsnp $dbsnp_del100 --cosmic $cosmic -contamination 0  --max_alt_alleles_in_normal_count 3 -cTable $out_dir/case.grp -nTable $out_dir/normal.grp --max_alt_alleles_in_normal_qscore_sum 40 --max_alt_allele_in_normal_fraction 0.02 -dt NONE -o fpga_speedup.vcf -nct 6 -dtt 2 -ntLib 36 --pair_hmm_implementation FPGA
#time java -d64 -server -XX:+UseG1GC -Xms3g -Xmx20g -Djava.io.tmpdir=../tmp -Djava.library.path=$lib_path -jar $gatk -T MuTect3 -R $reference_file -I:tumor $out_dir/case_sort_markdup_realign.bam -I:normal $out_dir/normal_sort_markdup_realign.bam -L /home/lq/software/script/target.bed --dbsnp $dbsnp_del100 --cosmic $cosmic -contamination 0 --max_alt_alleles_in_normal_count 3 --max_alt_alleles_in_normal_qscore_sum 40 --max_alt_allele_in_normal_fraction 0.02 -dt NONE -dtt 1 -nct 10 -o gatk_geneplus_pr_mutect_xilinx.vcf --pair_hmm_implementation FPGA 
####################################################################################################
#cho "##################################      END!!!!!!!!!!!!!!!      ##############################"
