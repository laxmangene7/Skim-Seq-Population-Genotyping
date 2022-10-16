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
hisat2-2.1.0/hisat2 -p 12 -x Hisat2.indexed.monoc.ref -1 TA4342_L95.fp.R1.fq.gz -2 TA4342_L95.fp.R2.fq.gz--no-spliced-alignment --no-unal -S parent1.sam
hisat2-2.1.0/hisat2 -p 12 -x Hisat2.indexed.monoc.ref -1 TA4342_L96.fp.R1.fq.gz -2 TA4342_L96.fp.R2.fq.gz--no-spliced-alignment --no-unal -S parent2.sam
```

Retriving concordant unique reads:
```
cat <(samtools view -H parent1.sam) <(awk '/YT:Z:CP/ && /NH:i:1/' parent1.sam) | samtools sort -o parent1.s.bam
samtools index -c parent1.s.bam
cat <(samtools view -H parent2.sam) <(awk '/YT:Z:CP/ && /NH:i:1/' parent2.sam) | samtools sort -o parent2.s.bam
samtools index -c parent2.s.bam
```

Variant calling:
```
bcftools mpileup --annotate AD,DP,INFO/AD --skip-indels -f monoc.ref.fasta -b bamFilesList.txt -B | bcftools call -m --variants-only  --skip-variants indels --output-type v -o monococcum.parents.RILs.vcf --group-samples -
```
Variant filter:

The variants called on parents were filtered so as to remove any loci with het genotype call, missing call and monomorphic loci


#### B. Genotyping of variants identified between parents in a recombinant inbred line (RIL) population.

The SNP positions are listed in a file which is used in BCFtools:
```
grep -v '^#' monococcum.parents.RILs.vcf | awk '{print $1"\t"$2"\t"$4","$5}' | bgzip -c > parentSNP_positions.tsv.gz && tabix -s1 -b2 -e2 parentSNP_positions.tsv.gz.  ## though we did not use -b2 flag here
```

Remove adapters from RILs skim-sequencing data as we did with parents WGS:
```
fastp -i sample.R1.fq -I sample.R2.fq -o sample.fp.R1.fq.gz -O sample.fp.R2.fq.gz --thread=5 --html=sample.html --json=sample.json --detect_adapter_for_pe --qualified_quality_phred=10 --length_required=150
```

The adapter trimmed reads from the RILs were aligned to the reference genome using Hisat2. The RIL's sorted bam file names are listed in `bamFile_list.txt` per line and used for genotyping using BCFtools:
```
bcftools mpileup -T parentSNP_positions.tsv.gz --annotate AD,DP,INFO/AD --skip-indels -f monoc.ref.fasta -b bamFile_list.txt -B | bcftools call -m --constrain alleles -T parentSNP_positions.tsv.gz --variants-only --skip-variants indels --output-type v -o monococcum.RILs.vcf --group-samples -
```

# Finding allelic disributions in RILs

https://user-images.githubusercontent.com/49244360/195334810-6f2bc70e-96e9-4e96-b420-882fd309f5b5.png

Merge RILs and Parents VCF so that two parents (P1 and P2) are in the last two columns:
```
module load  BCFtools
bgzip -c monococcum.RILs.vcf  > monococcum.RILs.vcf.gz
bgzip -c monococcum.parents.RILs.vcf > monococcum.parents.RILs.vcf.gz
bcftools merge *vcf.gz -Oz -o Merged.mono.parents.RILs.vcf.gz
unzip  Merged.mono.parents.RILs.vcf.gz > Merged.mono.parents.RILs.vcf
```

Convert Merge vcf file to txt file using BCFtools:

```
module load BCFtools
grep '#CHROM' Merged.mono.parents.RILs.vcf >   ## get header row of the txt file separately

cut --complement  -f3,6,7,8,9  Merged.mono.parents.RILs.vcf.header.row.txt > header.row.trimmed.txt  ## remove unwanted header columns

bcftools query -f '%CHROM %POS  %REF  %ALT [ %GT]\n' Merged.mono.parents.RILs.vcf > Merged.mono.parents.RILs.vcf.txt ## this file does not have header

cat header.row.trimmed.txt  Merged.mono.parents.RILs.vcf.header.row.txt >  mono.parents.RILs.genotyped.header.txt
```

Convert alleles as either P1, P2, or H: 

```
sed -e "/1$/s/0\/0/P1/g" -e "/1$/s/1\/1/P2/g" -e "/0$/s/1\/1/P1/g" -e "/0$/s/0\/0/P2/g" -e "s/0\/1/H/g"  mono.parents.RILs.genotyped.header.txt > mono.parents.RILs.alleles.identified.txt
```



