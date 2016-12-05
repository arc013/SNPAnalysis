vcf=$1
for file in $1; do
  for sample in `bcftools view -h $file | grep "^#CHROM" | cut -f10-`; do
	echo SAMPLE=$sample  
  bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
  done
done
