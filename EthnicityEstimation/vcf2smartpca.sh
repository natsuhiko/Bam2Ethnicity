
VCFDIR=$1
OUTDIR=$1/SmartPCA
mkdir $OUTDIR

TABIX=/software/team170/miniconda3/bin/tabix 

# converting vcf to eigen format
rm $OUTDIR/eigen.orig.gz
for CHR in {1..22}
do
	$TABIX $VCFDIR/chr$CHR.vcf.gz $CHR | awk '{print $1"\t"$2-1"\t"$2}' | gzip > $OUTDIR/tmp.bed.gz
	$TABIX -R $OUTDIR/tmp.bed.gz /nfs/team205/nk5/1000G/20181129/chr$CHR.vcf.gz | /nfs/team205/nk5/Applications/EthnicityEst/vcf2eigen CONV -N 2548 -dose -MAF 0.05 | gzip > $OUTDIR/tmp1.eigen.gz
	$TABIX -R $OUTDIR/tmp.bed.gz $VCFDIR/chr$CHR.vcf.gz | /nfs/team205/nk5/Applications/EthnicityEst/vcf2eigen CONV -N 1 -dose -MAF -1.0 | gzip > $OUTDIR/tmp2.eigen.gz
	join -1 1 -2 1 <(zcat $OUTDIR/tmp1.eigen.gz | sort -k1 -b) <(zcat $OUTDIR/tmp2.eigen.gz | sort -k1 -b) | sed 's/_/ /g' | sed 's/ /\t/g' | sort -k2 -n | gzip >> $OUTDIR/eigen.orig.gz
done
rm $OUTDIR/tmp1.eigen.gz $OUTDIR/tmp2.eigen.gz

# snp file
zcat $OUTDIR/eigen.orig.gz | awk '{print $1"_"$2,$1,0.0,$2,$3,$4}' > $OUTDIR/snp

# genotype file
zcat $OUTDIR/eigen.orig.gz | awk '{print $5$6}' > $OUTDIR/genotype

# individual file
for I in {1..2549}
do
       echo $I F Case
done > $OUTDIR/ind

# SMARTPCA
/nfs/team205/nk5/Applications/EIG/bin/smartpca.perl -i $OUTDIR/genotype -a $OUTDIR/snp -b $OUTDIR/ind -e $OUTDIR/eval -l $OUTDIR/log -o $OUTDIR/out -p $OUTDIR/plot -m 0
bash /nfs/team205/nk5/Applications/EthnicityEst/GMC.sh $OUTDIR

