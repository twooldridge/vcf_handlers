#!/bin/bash

invcf=$1
refgenome=$2
sampleID=$3

cat <<EOF > ${sampleID}.vcffilter.slurm
#!/bin/bash
#SBATCH -p hoekstra,commons,shared,serial_requeue
#SBATCH -t 56:00:00
#SBATCH --mem=20000
#SBATCH -n 2
#SBATCH -e ./logs/${sampleID}.var.select.e
#SBATCH -o ./logs/${sampleID}.var.select.o
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -J ${sampleID}.var.select

module load centos6/0.0.1-fasrc01
module load java
module load bcftools
#export PATH="/n/home01/lassance/Software/bcftools-1.7/:$PATH"
export BCFTOOLS_PLUGINS=/n/home01/lassance/Software/bcftools-1.7/plugins

### VARIANT ###
java -Xmx10G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/GATK/GenomeAnalysisTK-3.8.jar -T SelectVariants -R ${refgenome} -V ${invcf} --selectTypeToExclude NO_VARIATION -o ${sampleID}.var.vcf
bash /n/home11/twooldridge/scripts/zip_index_vcf.sh ${sampleID}.var.vcf
#
#bcftools view -i 'ALT !~ "*" && TYPE="snp"' ${sampleID}.var.vcf.gz -O v  -o ${sampleID}.snps.vcf 
#bcftools view -e 'ALT !~ "*" && TYPE="snp"' ${sampleID}.var.vcf.gz -O v  -o ${sampleID}.indels.vcf 
#
java -Xmx5G -jar ~/Software/picard/2.18.4/picard.jar SplitVcfs I=${sampleID}.var.vcf.gz SNP_OUTPUT=${sampleID}.snps.vcf INDEL_OUTPUT=${sampleID}.indels.vcf STRICT=false
#
java -Xmx5G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/GATK/GenomeAnalysisTK-3.8.jar -R ${refgenome}  -T VariantFiltration --variant ${sampleID}.snps.vcf --filterExpression "QD < 2.0 || MQ < 40.0 || FS > 10.0 || ReadPosRankSum < -20.0 || SOR > 3.0 || MQRankSum < -12.5" --filterName snps_filter -o sf.${sampleID}.snps.vcf 
~/scripts/zip_index_vcf.sh sf.${sampleID}.snps.vcf 
#
java -Xmx5G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/GATK/GenomeAnalysisTK-3.8.jar -R ${refgenome}  -T VariantFiltration --variant ${sampleID}.indels.vcf --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 3.0" --filterName indels_filter -o sf.${sampleID}.indels.vcf 
~/scripts/zip_index_vcf.sh sf.${sampleID}.indels.vcf 
#
bcftools concat --allow-overlaps sf.${sampleID}.snps.vcf.gz sf.${sampleID}.indels.vcf.gz -O v -o sf.${sampleID}.var.vcf 
#
~/scripts/zip_index_vcf.sh sf.${sampleID}.var.vcf 
#
#
#### INVARIANT ###
java -Xmx10G -XX:ParallelGCThreads=1 -Djava.io.tmpdir=`pwd`/tmp -jar ~/Software/GATK/GenomeAnalysisTK-3.8.jar -T SelectVariants -R ${refgenome} -V ${invcf} --selectTypeToInclude NO_VARIATION -o ${sampleID}.invar.vcf
bash /n/home11/twooldridge/scripts/zip_index_vcf.sh ${sampleID}.invar.vcf
#
bcftools filter --include '%QUAL>=20' --soft-filter poorQual -O v -o sf.${sampleID}.invar.vcf ${sampleID}.invar.vcf.gz 
~/scripts/zip_index_vcf.sh sf.${sampleID}.invar.vcf 
#
#
#### COMBINE AND DEPTH MASK ###
bcftools concat --allow-overlaps sf.${sampleID}.var.vcf.gz sf.${sampleID}.invar.vcf.gz -O v -o sf.${sampleID}.vcf 
bash ~/scripts/zip_index_vcf.sh sf.${sampleID}.vcf
#
bcftools +setGT sf.${sampleID}.vcf.gz -- -t q -e 'FMT/DP>4' -n . | bcftools view -f .,PASS > hf.${sampleID}.vcf 
bash ~/scripts/zip_index_vcf.sh hf.${sampleID}.vcf


### GET SNPS AND MONOMORMPHIC SITES### 
bcftools view hf.${sampleID}.vcf.gz -i 'ALT !~ "*"  && STRLEN(REF) == 1 && (TYPE="snp" || ALT == ".")' -M2 -O v -o ${sampleID}.snp.mono.vcf 
bash ~/scripts/zip_index_vcf.sh ${sampleID}.snp.mono.vcf 
EOF

sbatch ${sampleID}.vcffilter.slurm
