#!/usr/bin/env python3
__version__ = "0.0.1"
__author__ = "Sangjin Lee"

# modules
import ccs2mut.util
import ccs2mut.caller
import ccs2mut.vcflib 
from ccs2mut.parse_args import parse_args

def main():
    parser, options = parse_args(program_version=__version__)
    ccs2mut.util.check_num_threads(options.threads)
    ccs2mut.caller.call_somatic_substitutions(
        options.bam, # input # bam_file
        options.vcf, # germline mutations
        options.phased_vcf,  # phased germline mutations
        options.region, # target contigs/scaffolds/chromosomes
        options.region_list, # target contigs/scaffolds/chromosomes fofn
        options.ploidy, # haploid, diploid or polyploid
        options.min_mapq, # int: 0-60
        options.min_trim, # float: 0.01 - 0.1
        options.min_sequence_identity, # float: blast sequence identity
        options.min_hq_base_proportion, # float: proportion of BQ=93 bases
        options.min_alignment_proportion, # float: query alignment length / query read length
        options.common_snps, # common snps
        options.panel_of_normals, # panel of normals
        options.min_bq, # minimum base quality score: int
        options.mismatch_window, # mismatch window size
        options.max_mismatch_count, # maximum number of mismatches within a window
        options.min_ref_count, # number of reads supporting the reference base
        options.min_alt_count, # number of reads supporting the alterantive base
        options.min_hap_count, # number of reads supporting h0 and h1 haplotype
        options.threads, # maxminum number of threads
        options.phase, # bool
        __version__, # str
        options.out, # output # himut vcf file
    )


if __name__ == "__main__":
    main()
    ccs2mut.util.exit()
