# Skim-Seq-Population-Genotyping-for-Linkage-Mapping
Genotyping a larger plant population using skim-sequencing such as Nextera sequencing has gained popularity as a low-cost, efficient and high accuracy sequencing method. However, no proper analysis method to generate genotyping file from the raw sequence files have been documented in public domain.  This protocol describe how we call variants from skim-seq population when we have deep-sequencing data of the parental lines and the skim-sequencing data of the progenies or mapping population. 

# Generating Parent and Progeny VCF files

### A. Variant calling in whole-genome-sequencing (WGS) parents
  Two monococcum wheat parents;  wild (TA4342_L95) and domesticated (TA4342_L96) and their progenies (RIL7) were used to test the pipeline. 
  
 Trimming WGS data using fastp to remove adapters 
```
fastp -i sample.R1.fq -I sample.R2.fq -o sample.fp.R1.fq.gz -O sample.fp.R2.fq.gz --thread=7 --html=sample.html --json=sample.json --detect_adapter_for_pe --qualified_quality_phred=10 --length_required=150
```


The WGS data of two monococcum parents were aligned to the TA299 wild monoroccom reference genome using Hisat2:
```
hisat2-2.1.0/hisat2 -p 12 -x TA299_validated_v1.0.monoc.ref -1 TA4342_L95.fp.R1.fq.gz -2 TA4342_L95.fp.R2.fq.gz--no-spliced-alignment --no-unal -S TA4342_L95.sam
hisat2-2.1.0/hisat2 -p 12 -x TA299_validated_v1.0.monoc.ref -1 TA4342_L96.fp.R1.fq.gz -2 TA4342_L96.fp.R2.fq.gz--no-spliced-alignment --no-unal -S TA4342_L96.sam
```

Retriving concordant unique reads:
```
cat <(samtools view -H TA4342_L95.sam) <(awk '/YT:Z:CP/ && /NH:i:1/' TA4342_L95.sam) | samtools sort -o TA4342_L95.s.bam
samtools index -c TA4342_L95.s.bam
cat <(samtools view -H TA4342_L96.sam) <(awk '/YT:Z:CP/ && /NH:i:1/' TA4342_L96.sam) | samtools sort -o TA4342_L96.s.bam
samtools index -c TA4342_L96.s.bam
```

Variant calling:
```
bcftools mpileup --annotate AD,DP,INFO/AD --skip-indels -f TA299_validated_v1.0.monoc.ref.fasta -b bamFilesList.txt -B | bcftools call -m --variants-only  --skip-variants indels --output-type v -o monococcum.parents.with.TA299.vcf --group-samples -
```
Variant filter:

The variants called on parents were filtered so as to remove any loci with het genotype call, missing call and monomorphic loci


#### B. Genotyping of variants identified between parents in a recombinant inbred line (RIL) population.

The SNP positions are listed in a file which is used in BCFtools:
```
grep -v '^#' monococcum.parents.with.TA299.vcf | awk '{print $1"\t"$2"\t"$4","$5}' | bgzip -c > parentSNP_positions.tsv.gz && tabix -s1 -b2 -e2 parentSNP_positions.tsv.gz
```

Remove adapters from RILs skim-sequencing data as we did with parents WGS:
```
fastp -i sample.R1.fq -I sample.R2.fq -o sample.fp.R1.fq.gz -O sample.fp.R2.fq.gz --thread=5 --html=sample.html --json=sample.json --detect_adapter_for_pe --qualified_quality_phred=10 --length_required=150
```

The adapter trimmed reads from the doubled haploid lines were aligned to the reference genome using Hisat2 and filtered to recover unique concordant reads as parent. The 48 sorted bam file names are listed in `bamFile_list.txt` per line and used for genotyping using BCFtools:
```
bcftools mpileup -T parentSNP_positions.tsv.gz --annotate AD,DP,INFO/AD --skip-indels -f TA299_validated_v1.0.monoc.ref.fasta -b bamFile_list.txt -B | bcftools call -m --constrain alleles -T parentSNP_positions.tsv.gz --variants-only --skip-variants indels --output-type v -o monococcum.RILs.vcf --group-samples -
```

# Finding allelic disributions in RILs

https://user-images.githubusercontent.com/49244360/195334810-6f2bc70e-96e9-4e96-b420-882fd309f5b5.png

Merge RILs and Parents VCF so that two parents (P1 and P2) are in the last two columns:
```
module load  BCFtools
bgzip -c monococcum.RILs.vcf  > monococcum.RILs.vcf.gz
bgzip -c monococcum.parents.with.TA299.vcf > monococcum.parents.with.TA299.vcf.gz
bcftools merge *vcf.gz -Oz -o Merged.RIL.aegilorefTA299.validated.vcf.gz
unzip  Merged.RIL.aegilorefTA299.validated.vcf.gz > Merged.RIL.aegilorefTA299.validated.vcf
```

Convert Merge vcf file to txt file using BCFtools:

```
module load BCFtools
grep '#CHROM' Merged.RIL.aegilorefTA299.validated.vcf > Merged.RIL.aegilorefTA299.validated.header.row.txt ``` ## get header row of the txt file separately
cut --complement  -f3,6,7,8,9. Merged.RIL.aegilorefTA299.validated.header.row.txt > Merged.RIL.aegilorefTA299.validated.header.row.trimmed.txt ## remove unwanted header columns
bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' Merged.RIL.aegilorefTA299.validated.vcf > Merged.RIL.aegilorefTA299.validated.vcf.txt
cat Merged.RIL.aegilorefTA299.validated.header.row.trimmed.txt  Merged.RIL.aegilorefTA299.validated.vcf.txt >  Merged.RIL.aegilorefTA299.validated.genotyped.header.txt
```
