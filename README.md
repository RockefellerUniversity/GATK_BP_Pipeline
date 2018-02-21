# GATK_BP_Pipeline

A variants call pipeline based on GATK best practise guidance:

bwaAlign
 
 bwa mem -R "@RG\tID:$postprocID\tSM:$sampleID\tLB:$sampleID\tPL:illumina" -M -t 4 $reference $R1.fastq.gz $R2.fastq.gz |java -jar -Djava.io.tmpdir=$TMPDIR -Xmx4G $PICARDTOOLS SortSam INPUT=/dev/stdin OUTPUT=picard.sorted.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate TMP_DIR=$TMPDIR 
 
 
 picardMarkDup
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx32G $PICARDTOOLS MarkDuplicates $inputfiles REMOVE_DUPLICATES=false CREATE_INDEX=true ASSUME_SORTED=true OUTPUT=picard.sort.merged.rmdup.bam METRICS_FILE=picard.sort.merged.rmdup.bam.metric_file 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx32G $PICARDTOOLS BuildBamIndex INPUT=picard.sort.merged.rmdup.bam
 
 picardCalculateHsMetrics  ( A QC step )
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx16G $PICARDTOOLS CollectHsMetrics VALIDATION_STRINGENCY=SILENT INPUT=picard.sort.merged.rmdup.bam OUTPUT=picard.hs.metrics.txt BAIT_INTERVALS=$intervalFile TARGET_INTERVALS=$intervalFile R=$reference
head -n 8 picard.hs.metrics.txt  | tsp > picard.hs.metrics.tsp.txt
 
 gatkQscoreRecalibration
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T BaseRecalibrator -nct $threads -I picard.sort.merged.rmdup.bam -o recal_data.table -R $reference -l INFO -knownSites $dbSNP -knownSites $g1k_indel_vcf -knownSites $g1k_indel_phase1_vcf 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T PrintReads -nct $threads -I picard.sort.merged.rmdup.bam -o realigned-recalibrated.bam -R $reference -l INFO -BQSR recal_data.table 
 
 gatkRawVariantsCall (run by each chrom)

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T HaplotypeCaller -R $reference -I realigned-recalibrated.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -rf BadCigar --dontUseSoftClippedBases --min_base_quality_score 20 --emitRefConfidence GVCF -L ${chr}' -o raw_variants.chr${chr}.g.vcf --dbsnp $dbSNP 

or for high coverage samples ( > 250X )

java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T HaplotypeCaller -R $reference -I realigned-recalibrated.bam --genotyping_mode DISCOVERY -dfrac 0.99 -stand_emit_conf 10 -stand_call_conf 30 -rf BadCigar --min_base_quality_score 20 --emitRefConfidence GVCF -nct 4 -L ${chr} -o raw_variants.${chr}.g.vcf --dbsnp $dbSNP  

 gatkRawVariants (combination of all chrom)
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T CombineGVCFs all_chr_files("--variant ", "raw_variants.${chr}.g.vcf") -R $reference -o allchr.g.vcf 
 
 gatkJointGenotyping (snp and indels only output)
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $GATK -T GenotypeGVCFs -R $reference all_chr_files("--variant ", "raw_variants.${chr}.g.vcf") -o gatk.snp.indel.vcf" 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx16G $GATK -T VariantEval -L $captureKitFile -o raw.eval.txt --eval gatk.snp.indel.vcf --dbsnp $dbSNP -R $reference 
bedtools intersect -a gatk.snp.indel.vcf -b $captureKitFile -wa > gatk.snp.indel.target.vcf 
 
 annovar
 
convert2annovar.pl -format vcf4 gatk.snp.indel.vcf > gatk.snp.indel.annovar 
table_annovar.pl --protocol dbnsfp30a,segdup,cg46,avsnp147,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,clinvar_20160302,refgene,cosmic70,ensGene,exac03 --operation f,r,f,f,f,f,f,f,f,f,f,f,f,f,g,f,g,f --buildver hg19 --remove --nastring gatk.snp.indel.annovar $annovarHDB
 
 snpEff
 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $snpEff eff -v -i vcf -o vcf -c $snpEffConfig -spliceSiteSize 7 hg19 gatk.snp.indel.vcf > var.annotated.refseq.vcf 
java -jar -Djava.io.tmpdir=$TMPDIR -Xmx24G $snpEff eff -v -motif -nextprot -i vcf -o vcf -c $snpEffConfig GRCh37.75 gatk.snp.indel.vcf > var.annotated.ens.vcf 
