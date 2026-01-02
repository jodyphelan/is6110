# IS6110

This tool is designed to identify the presence of insertion sequences (e.g. IS6110) in _Mycobacterium tuberculosis_ genomes. It uses a BAM file as input and outputs a VCF file with the identified regions. 

## Motivation
IS6110 is a mobile genetic element that plays a significant role in the evolution and adaptation of Mycobacterium tuberculosis. It has been implecated in drug resistance through disruption of genes such as _mmpR5_ (https://pubmed.ncbi.nlm.nih.gov/34287057/). Accurate identification of IS6110 insertions is crucial to perform accurate genomic drug resistance prediction. Though there are existing tools to identify IS6110 insertions, I wanted to create a lightweight and easy-to-use tool that could be easily integrated into existing pipelines. Key to this is to:

 1. Take input in BAM format to avoid realigning to the reference genome
 2. Output in VCF format to allow for easy integration into downsteam analysis pipelines.

At first the tool is designed to identify IS6110 insertions, but I have adapted it to be able to identify other insertion sequences by providing a custom FASTA file of the IS of interest. I might change the name later to reflect this.

> [!IMPORTANT]
> The tool is still under development and may not be fully functional. It does not detect deletions of IS at this time. For a more comprehansive characterisation including deletions, see [ISMapper](https://https://github.com/jhawkey/IS_mapper).

## Install

```bash
mamba create -n is6110 is6110 -c conda-forge -c bioconda
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

## Output

The output is a VCF file. If a GFF file is provided, the tool will annotate the VCF file with gene disruption information. Genes which are disrupted by an IS insertion will be annotated with a similar format to that used by SnpEff, and with the `INS:ME:<IS_NAME>` as the alternate. INS:ME is used to specify insertion of a mobile element relative to the reference (see [VCF specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)).
```
Chromosome      3074302 .       T       <INS:ME:IS6110> .       .       END=3074302;ANN=<INS:ME:IS6110>|transcript_ablation|HIGH|Rv2764c|thyA|transcript|transcript:Rv2764c|protein_coding|1/1|c.170_171_ins[IS6110]||||||    GT:AD   1:1,241
```