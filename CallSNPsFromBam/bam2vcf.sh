#!/bin/bash

BAM=$1
OUTDIR=$2
if [ -z $3 ]
then
	SAMPLENAME=$1
else
	SAMPLENAME=$3
fi

mkdir $OUTDIR

# user defined variables
DIR=/path/to/Bam2Ethnicity
SAMTOOLS=/software/npg/current/bin/samtools
BGZIP=/software/team170/miniconda3/bin/bgzip
TABIX=/software/team170/miniconda3/bin/tabix


for CHR in {1..22}
do
	$SAMTOOLS view -F 0x0100 $BAM $CHR | \
        	$DIR/QCFilterBam/qcFilterBam stdin -skipMissing F -maxMismatch 3 -maxGapOpen 0 -maxBestHit 1 -minQual 10 -singleEnd | \
        	sort -k 6 | awk 'BEGIN{CBUB="X";OFS="\t"}{if(CBUB!=$6){print $1,$2,$3,$4,$5; CBUB=$6}}' | sort -k 2,2n | \
        	$DIR/QCFilterBam/countAS /nfs/team205/nk5/1000G/20181129/chr"$CHR"_SNP_header.gz | \
		awk '{print $5","$6}' | gzip > $OUTDIR/tmp.gz

	$DIR/QCFilterBam/zpaste /nfs/team205/nk5/1000G/20181129/chr"$CHR"_SNP_header.gz $OUTDIR/tmp.gz | \
		awk 'BEGIN{OTF="\t"}{for(i=1;i<=8;i++){printf $i"\t"};printf "AS\t"; print $10}' | $BGZIP > $OUTDIR/chr$CHR.AS.gz
	$TABIX -f -s 1 -b 2 -e 2 $OUTDIR/chr$CHR.AS.gz
	
	cp $DIR/CallSNPsFromBam/header $OUTDIR/header1
	echo $SAMPLENAME > $OUTDIR/id.txt
	awk '{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"$1}' $OUTDIR/id.txt >> $OUTDIR/header1
	$DIR/CallSNPsFromBam/call_genotype $OUTDIR/chr$CHR.AS.gz $CHR $DIR/CallSNPsFromBam/th.txt $OUTDIR/header1 -gt -gp -as --min-cov 4.0 | $BGZIP > $OUTDIR/chr$CHR.vcf.gz;
	$TABIX -f -s 1 -b 2 -e 2 $OUTDIR/chr$CHR.vcf.gz
done


rm $OUTDIR/header1 $OUTDIR/tmp.gz $OUTDIR/id.txt $OUTDIR/chr*.AS.gz $OUTDIR/chr*.AS.gz.tbi

