# is6110

This tool is designed to identify the presence of the IS6110 insertion sequence in Mycobacterium tuberculosis genomes. It uses a BAM file as input and outputs a VCF file with the identified regions. 

> [!IMPORTANT]
> The tool is still under development and may not be fully functional.

## Install

```bash
mamba create -n is6110 pathogen-profiler 
pip install git+https://github.com/jodyphelan/is6110.git
```

## Use

```bash
# With a single region (mmpR5)
is6110 -b <file.bam> -o out.vcf -r <reference.fasta> -g <genes.gff> -r Chromosome:778385-779715
# With a list of regions in bed format
is6110 -b <file.bam> -o out.vcf -r <reference.fasta> -g <genes.gff> -T test.bed
# Across whole genome
is6110 -b <file.bam> -o out.vcf -r <reference.fasta> -g <genes.gff>
```