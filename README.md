# SimpleHunter
A light-weight java program to find loci of interest in a genomics mutant screen. 


This code was written to support Emily Tipper (JIC) for a mutational genomics experiment in her PhD thesis.
It is designed to find loci of interest in a mutant screen. SNPs are called in individual samples as well as a wildtype self-mapping. Then data is integrated highlighting loci in a genome where multiple samples have SNPs at different positions.

## Prerequisites:

Burrows-Wheeler Aligner [https://bio-bwa.sourceforge.net/](https://bio-bwa.sourceforge.net/)
Samtools [http://www.htslib.org/](http://www.htslib.org/)

Java Runtime Environment jre-1.8 [https://www.java.com/en/download/manual.jsp](https://www.java.com/en/download/manual.jsp)


### Step 1:

Run BWA. Then sort, index and convert to mpileup format using samtools.


#### create bwa index and fasta index for reference
```
bwa index <reference.fasta>
samtools faidx <reference.fasta>
```

#### for each samples, create mpileup file

```
bwa mem <reference.fasta> <sample_read1.fq> <sample_read2.fq> | samtools sort -o sample.bam -
samtools index <sample.bam>
samtools mpileup -BQ0 -f <reference.fasta> <sample.bam> > <sample.mpileup>
```


### Step 2:

Call SNPs in individual samples. Download QuickSNP from [https://github.com/steuernb/oat_mutseq](https://github.com/steuernb/oat_mutseq) 


```
java -jar QuickSNP.jar -i <sample.mpileup> -o <sample.snp.txt> -a <maximum_reference_allele_frequency> -c <minimum_coverage>
```

Options:

Parameter | | Description  
--- | --- | ---
 -i	| String	 | Input File in mpileup format
 -o	 | String	 | Output File (default:stdout)
 -a	 | float	 | Maximum reference allele frequency (default 0.1)
 -c	 | int	 | Minimum coverage to record (default: 5)
 -e		 |  | only output EMS specific SNPs. These are C->T and G->A


### Step 3

Integrate SNP calls

```
java -jar SimleHunter.jar -w <wt.snps.txt> -i <m1.snps.txt> [<m2.snps.txt> ...] -o <outputFile> -l <window_size> -m <min_number_mutants> -a <maximum_reference_allele_frequency> -n <maxNumberMutantsSharingPosition>
```

Options:

Parameter | | Description
--- | --- | ---
 -w	| String	 | SNP calls from self mapping of wildtyp
 -i	 | List of Strings	 | SNP calls from individual mappings of mutants
 -l | int | window size; the length of a locus to be reported
 -m | int | minimum number of mutants having a SNP in a window to report the locus
 -a	 | float	 | Maximum reference allele frequency
 -n	 | int	 | maximum number of mutants that may share a SNP position. Otherwise position will be disregarded
 -g		 | String | only regard SNPs in exons specified by this input file in GFF format

 
