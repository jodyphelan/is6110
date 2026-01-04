"""
Code to perform IS6110 variant calling
"""

__version__ = "0.3.0"


import argparse
import shutil
import pysam
from collections import Counter
import subprocess as sp
from typing import List, Tuple
import tempfile
import logging
from .gff import load_gff
import re
import os
from tqdm import tqdm

Item = Tuple[str, int]  # ('R' or 'L', position)
def cluster_positions(items: List[Item], min_distance: int = 5) -> List[List[Item]]:
    if not items:
        return []

    items = sorted(items, key=lambda x: x[1])

    clusters: List[List[Item]] = []
    letters_in_cluster: List[set[str]] = []

    for side, pos in items:
        if not clusters:
            clusters.append([(side, pos)])
            letters_in_cluster.append({side})
            continue

        last_cluster = clusters[-1]
        last_letters = letters_in_cluster[-1]
        last_pos = last_cluster[-1][1]

        within_dist = (pos - last_pos) <= min_distance
        has_opposite = (('L' in last_letters) if side == 'R' else ('R' in last_letters))

        if within_dist and has_opposite:
            last_cluster.append((side, pos))
            last_letters.add(side)
        else:
            clusters.append([(side, pos)])
            letters_in_cluster.append({side})

    coordinates = []
    for cluster in clusters:
        if len(cluster) > 1:
            coordinates.append( (cluster[0][1], cluster[-1][1]) )
    return sorted(coordinates)


def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, _ = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True
    return False

def run_cmd(cmd: str, desc=None, log: str=None, exit_on_error: bool=True) -> sp.CompletedProcess:
    if desc:
        logging.info(desc)
    processed_cmd = cmd.replace("&&","XX")
    programs = set([x.strip().split()[0] for x in re.split("[|&;]",processed_cmd.strip()) if x!=""])
    missing = [p for p in programs if which(p)==False]
    if len(missing)>0:
        raise ValueError("Cant find programs: %s\n" % (", ".join(missing)))
    logging.debug(f"Running command: {cmd}")
    cmd = "/bin/bash -c set -o pipefail; " + cmd
    output = open(log,"w") if log else sp.PIPE
    result = sp.run(cmd,shell=True,stderr=output,stdout=output)
    if result.returncode != 0:
        logging.error(result.stderr.decode("utf-8"))
        if exit_on_error:
            raise ValueError("Command Failed:\n%s\nstderr:\n%s" % (cmd,result.stderr.decode()))
    return result




def get_gene_position(gene, position):
    """
    Get the gene position for a given position.
    """
    results = {}
    for transcript in gene.transcripts:
        i,exon = [(i,exon) for i,exon in enumerate(transcript.exons) if exon.start <= position <= exon.end][0]
        if exon.strand=="+":
            exon_relative_position = position - exon.start
            gene_relative_position = sum([e.end - e.start for e in transcript.exons[:i]]) + exon_relative_position
        else:
            exon_relative_position = exon.end - position
            gene_relative_position = sum([e.end - e.start for e in transcript.exons[i+1:]]) + exon_relative_position
        results[transcript.name] = gene_relative_position
    return results

def get_overlapping_genes(gff, position):
    """
    Get the genes overlapping with a given position.
    """
    overlapping_genes = [g for g in gff if g.feature_start <= position <= g.feature_end]
    return overlapping_genes


        
def get_hgvs(gene_position, me_name):
    """
    Get the HGVS notation for a given gene position.
    """
    return f'c.{gene_position}_{gene_position+1}_ins[{me_name}]'


def get_annotation(gff: list, var: pysam.VariantRecord, me_name:str):
    """
    Get the annotation for a given variant.
    """
    anns = []
    genes = get_overlapping_genes(gff, var.start)

    for gene in genes:
        gene_positions = get_gene_position(gene, var.start)
        for transcript,pos in gene_positions.items():
            hgvs = get_hgvs(pos, me_name)
            ann = f'{var.alts[0]}|transcript_ablation|HIGH|{gene.gene_id}|{gene.name}|transcript|{transcript}|protein_coding|1/1|{hgvs}||||||'
            anns.append(ann)
    
    return anns




def get_is_overlapping_reads(
        bam_file: str, 
        seq: str, 
        region: str = None, 
        cores: int = 1, 
        min_seed_len: int = 11, 
        min_aln_score: int = 20
    ) -> dict:
    """
    Get reads overlapping with IS sequences. Returns a dictionary mapping IS names to sets of read names.
    """
    read_names_by_is = {}
    with tempfile.TemporaryDirectory() as tmpdirname:
        logging.debug(f"Temporary directory created: {tmpdirname}")
        # copy the seq file to the temp directory
        # new file name = tmpdirname + IS.fasta
        new_seq_path = os.path.join(tmpdirname, "IS.fasta")
        shutil.copy(seq, new_seq_path)
        seq = new_seq_path
        # index the seq file
        run_cmd(f'bwa index {seq}')
        if region is None:
            first_part = f'extract-clipped-reads -b {bam_file} '
        else:
            first_part = f'samtools view -b {bam_file} {region}'
        run_cmd(f'{first_part} | samtools fastq | bwa mem -k {min_seed_len} -T {min_aln_score} -t {cores} {seq} - | samtools sort -o {tmpdirname}/temp.bam -')
        run_cmd(f'samtools index {tmpdirname}/temp.bam')
        bam = pysam.AlignmentFile(f'{tmpdirname}/temp.bam', "rb")
        for read in bam:
            if not read.is_unmapped:
                is_name = read.reference_name
                if is_name not in read_names_by_is:
                    read_names_by_is[is_name] = set()
                read_names_by_is[is_name].add(read.query_name)
    return read_names_by_is

def get_most_common_positions(bam_file: str, region: str, read_names: List[str], min_depth: int = 10) -> List[int]:
    bam = pysam.AlignmentFile(bam_file)
    positions = []
    
    chrom = region.split(':')[0]
    start,end = region.split(':')[1].split('-')
    start = int(start)
    end = int(end)

    for read in bam.fetch(chrom, start, end):
        if read.query_name in read_names:
            if not read.cigarstring:
                continue
            if 'S' in read.cigarstring:
                if re.search(r'^[0-9]+S', read.cigarstring):
                    positions.append(('L',read.reference_start))
                elif re.search(r'[0-9]+S$', read.cigarstring):
                    positions.append(('R',read.reference_end))
                else:
                    pass


    # Filter positions based on minimum depth
    position_counts = Counter(positions)
    logging.debug(f"Position counts: {position_counts}")
    positions = sorted(
        [pos for pos, count in position_counts.items() if count >= min_depth]
        , key=lambda x: x[1]
    )

    
    return positions

def get_AD(bam_file: str, chrom: str, start: int, end: int = None):
    """
    Get the reads overlapping with the specified region in the BAM file.
    """

    if end is None:
        end = start + 1

    bam = pysam.AlignmentFile(bam_file)
    covering_reads = []
    clipped_reads = []

    logging.debug(f"Fetching reads from {chrom}:{start}-{end}")
    for read in bam.fetch(chrom, start-1, end):
        # if read is soft clipped 
        if read.cigarstring and 'H' in read.cigarstring:
            continue
        elif read.cigarstring and 'S' in read.cigarstring:
            clipped_reads.append(read)
        else:
            # check if read spans the region
            if read.reference_end is None or read.reference_start is None:
                continue
            if read.reference_start <= start and read.reference_end >= end:
                covering_reads.append(read)
            else:
                # logging.debug(f"Read {read.query_name} does not span the region {chrom}:{start}-{end}")
                pass

    

    return (len(covering_reads), len(clipped_reads))

def get_insertion_position(position_counts: Tuple[int,int]):
    """
    Get the most common insertion position from the list of positions.
    """
    depth = sum(count for _, count in position_counts)
    return (position_counts[0][0], depth)

def write_vcf(
    positions: List[Tuple[str,int,Tuple[int,int],str]], 
    out_file: str, ref: pysam.FastaFile, 
    sample_name: str, 
    is_names: List[str],
    gff_file: str = None
):
    """
    Write the insertion positions to a VCF file using pysam
    """
    chromosome_names = ref.references
    with pysam.VariantFile(out_file, "w") as vcf_out:
        vcf_out.header.add_meta('INFO', items=[('ID', 'IS'), ('Number', '.'), ('Type', 'String'), ('Description', 'Insertion Sequence')])
        vcf_out.header.add_meta('FORMAT', items=[('ID', 'GT'), ('Number', '1'), ('Type', 'String'), ('Description', 'Genotype')])
        vcf_out.header.add_meta('FORMAT', items=[('ID', 'AD'), ('Number', 'R'), ('Type', 'Integer'), ('Description', 'Number of reads supporting the reference and insertion alleles')])

        vcf_out.header.add_meta('FORMAT', items=[('ID', 'CR'), ('Number', '1'), ('Type', 'Integer'), ('Description', 'Number of clipped reads supporting the insertion')])
        for is_name in is_names:
            vcf_out.header.add_meta('ALT', items=[('ID', f'INS:ME:{is_name}'), ('Description', f'Insertion of {is_name} element')])
        if gff_file:
            vcf_out.header.add_meta('INFO', items=[('ID', 'ANN'), ('Number', '.'), ('Type', 'String'), ('Description', "Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ")])
            gff = load_gff(gff_file)


        vcf_out.header.add_sample(sample_name)
        
        for chrom in chromosome_names:
            chromlen = ref.get_reference_length(chrom)
            vcf_out.header.contigs.add(chrom, chromlen)
        
        for chrom,pos,ad,is_name in positions:
            pos = pos - 1
            ref_seq = ref.fetch(chrom, pos, pos+1)
            record = vcf_out.new_record(
                contig=chrom,
                start=pos,
                alleles=[ref_seq,f'<INS:ME:{is_name}>'],
                id='.',
            )
            record.samples[sample_name]['GT'] = (1)
            record.samples[sample_name]['AD'] = ad
            if gff_file:
                record.info['ANN'] = get_annotation(gff, record, is_name)
            vcf_out.write(record)

def get_sample_name(bam_file: str):
    """
    Get the sample name from the BAM file.
    """
    bam = pysam.AlignmentFile(bam_file)
    sample_name = bam.header['RG'][0]['SM']
    return sample_name


def cli_get_clipped_reads():
    parser = argparse.ArgumentParser(description="Get clipped reads from a BAM file.")
    parser.add_argument(
        "-b", "--bam", type=str, required=True, help="Input BAM file."
    )

    args = parser.parse_args()

    inbam = pysam.AlignmentFile(args.bam,'rb')
    # create filehandle to stdout
    outbam = pysam.AlignmentFile("/dev/stdout", "wb", template=inbam)

    for read in inbam:
        # check if read is clipped
        if read.cigarstring and 'S' in read.cigarstring:
            outbam.write(read)
    outbam.close()
    inbam.close()



def cli():


    parser = argparse.ArgumentParser(description="Find IS elements in a BAM file.")
    parser.add_argument(
        "-b", "--bam", type=str, required=True, help="Input BAM file."
    )
    parser.add_argument(
        "-o", "--out", type=str, required=True, help="Output file for IS counts."
    )
    parser.add_argument(
        "-r", "--ref", type=str, required=True, help="Reference genome file."
    )
    parser.add_argument(
        "-g", "--gff", type=str, required=False, help="Reference genome file."
    )
    parser.add_argument(
        "--is-sequences", type=str, required=False, help="Custom FASTA file containing IS sequences. If not provided, uses the default IS.fasta."
    )

    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument(
        "-t", "--target", type=str, required=False, help="Region to search in the format 'chr:start-end'."
    )
    group.add_argument(
        "-T", "--targets-file", type=str, required=False, help="Bed file with regions to search."
    )

    parser.add_argument(
        "-d", "--min-depth", type=int, default=10, help="Minimum depth for positions to be considered."
    )
    parser.add_argument(
        "--clipping-gap", type=int, default=10, help="Maximum gap between left and right clipped reads."
    )
    parser.add_argument(
        "--min-seed", type=int, default=10, help="Minimum seed length for BWA MEM (-k)."
    )
    parser.add_argument(
        "--min-score", type=int, default=15, help="Minimum alignment score for BWA MEM (-T)."
    )
    parser.add_argument(
        "-c", "--cpus", type=int, default=1, help="Number of CPUs to use."
    )
    parser.add_argument(
        "--debug", action="store_true", help="Enable verbose logging."
    )
    parser.add_argument(
        "--version", action="version", version=f"is6110 version {__version__}"
    )

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    else:
        # Set to WARNING to suppress debug messages
        logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')



    ref = pysam.FastaFile(args.ref)
    ref_seqname = ref.references[0]

    #get directory where this file is located
    dir_path = os.path.dirname(os.path.realpath(__file__))
    
    # Use user-provided IS sequences file if specified, otherwise use default
    if args.is_sequences:
        args.seq = args.is_sequences
        logging.info(f"Using custom IS sequences file: {args.seq}")
    else:
        args.seq = f'{dir_path}/data/IS.fasta'
        logging.info(f"Using default IS sequences file: {args.seq}")
    
    insertion_sequence = pysam.FastaFile(args.seq)
    is_names = list(insertion_sequence.references)
    
    logging.info(f"Detecting insertions for {len(is_names)} IS elements: {', '.join(is_names)}")
    
    sample_name = get_sample_name(args.bam)

    if args.targets_file:
        regions = []
        for l in open(args.targets_file):
            row = l.strip().split('\t')
            regions.append(row[0]+":"+row[1]+"-"+row[2])
    elif args.target:
        regions = [args.target]
    else:
        regions = []


    insertion_positions = set()
    
    # Process each IS sequence
    if len(regions)>0:
        for region in tqdm(regions, desc="Scanning regions"):
            # Map reads once for all IS elements
            is_overlapping_reads_by_is = get_is_overlapping_reads(args.bam, args.seq, region, min_seed_len=args.min_seed, min_aln_score=args.min_score, cores=args.cpus)
            
            # Process each IS element with its specific reads
            for is_name in is_names:
                if is_name not in is_overlapping_reads_by_is:
                    logging.debug(f"No reads mapped to {is_name} in {region}")
                    continue
                    
                is_overlapping_reads = is_overlapping_reads_by_is[is_name]
                positions = get_most_common_positions(args.bam, region, is_overlapping_reads, min_depth=args.min_depth)
                logging.debug(f"Positions for {is_name} in {region}: {positions}")

                if len(positions)>0:
                    positions = sorted(positions)
                    start = positions[0]
                    end = positions[1] if len(positions)==2 else None
                    
                    AD = get_AD(args.bam, ref_seqname, start, end)
                    insertion_positions.add((ref_seqname,start,AD,is_name))
    else:
        # Map reads once for all IS elements (whole genome)
        logging.info("Mapping reads to IS sequences...")
        is_overlapping_reads_by_is = get_is_overlapping_reads(args.bam, args.seq,  cores=args.cpus, min_seed_len=args.min_seed, min_aln_score=args.min_score)
        

        
        def get_clustered_positions(positions):
            clustered = set()
            for i, (direction, pos) in enumerate(positions):
                for j, (other_direction, other_pos) in enumerate(positions):
                    if (
                        (i != j and abs(pos - other_pos) <= args.clipping_gap)
                        and (direction != other_direction)
                    ):
                        clustered.add(pos)
                        break
            clustered = sorted(list(clustered))
            coordinates = []
            while len(clustered) > 0:
                if len(clustered) == 1:
                    raise ValueError("Unpaired clipping position found.")
                start = clustered.pop(0)
                end = clustered.pop(0)
                coordinates.append((start, end))
            return sorted(coordinates)
        
        # Process each IS element with its specific reads
        for is_name in is_names:
            logging.info(f"Processing {is_name}...")
            if is_name not in is_overlapping_reads_by_is:
                logging.info(f"No reads mapped to {is_name}")
                continue
                
            is_overlapping_reads = is_overlapping_reads_by_is[is_name]
            logging.info(f"Found {len(is_overlapping_reads)} reads mapping to {is_name}")
            positions = get_most_common_positions(args.bam, f'{ref_seqname}:1-{ref.get_reference_length(ref_seqname)}', is_overlapping_reads, min_depth=args.min_depth)
            logging.debug(f"Positions for {is_name}: {positions}")


            clustered_positions = cluster_positions(positions, min_distance=args.clipping_gap)
            logging.debug(f"Clustered Positions for {is_name}: {clustered_positions}")
            if len(clustered_positions)>0:
                for start,end in clustered_positions:      
                    AD = get_AD(args.bam, ref_seqname, start, end)
                    insertion_positions.add((ref_seqname,start,AD,is_name))
    
    insertion_positions = sorted(insertion_positions, key=lambda x: x[1])
    logging.info(f"Found {len(insertion_positions)} total insertion(s)")

    write_vcf(insertion_positions, args.out, ref, sample_name, is_names, gff_file=args.gff)