# How to create a bam file:

```
bowtie2-build GENOME_FILE.fa GENOME

bowtie2 -L 18 -N 0 -x GENOME  -U $i  -S $i.sam

i=FASTQ_FILES.fastq

samtools view -S -b $i > $i.bam

samtools merge tissue.bam *.bam
```


* Adapters should be removed  before these steps. You can use **cutadapt** to do that: https://github.com/jnktsj/DNApi

