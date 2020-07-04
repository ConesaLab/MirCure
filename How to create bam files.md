# How to create a bam file:

```
Bowtie2-build GRCh38.fa humanGenome
bowtie2 -L 18 -N 0 -x humanGenome  -U $i  -S $i.sam
samtools view -S -b $i > $i.bam
samtools merge tissue.bam *.bam
```
**1.** You need to remove adaptor before this step. You can use **cutadapt** to do that: https://github.com/jnktsj/DNApi
**2.** **GRCh38.fa** is the reference genome in FASTA format; you need to change to your FASTA genome file;**humanGenome** is the index name, and you can type a name you want.
**3.** You need to change to your reference genome index(not **humanGenome**)
**4.** **$I** is the RNA-seq data in FASTA or FASTQ format. You need to change to your own RNA-seq data (not *cut). **$I.sam** is the sam file name, you can choose what you want.

