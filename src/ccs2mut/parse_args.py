# modules
import sys
import warnings
import argparse

def make_wide(formatter, w=120, h=36):
    """Return a wider HelpFormatter, if possible."""
    try:
        # https://stackoverflow.com/a/5464440
        # beware: "Only the name of this class is considered a public API."
        kwargs = {'width': w, 'max_help_position': h}
        formatter(None, **kwargs)
        return lambda prog: formatter(prog, **kwargs)
    except TypeError:
        warnings.warn("argparse help formatter failed, falling back.")
        return formatter

# argparse
def parse_args(program_version, arguments=sys.argv[1:]):
    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=make_wide(argparse.ArgumentDefaultsHelpFormatter),
        description="ccs2mut calls single-base substitutions, double-base substitutions, indelsfrom PacBio CCS reads",
    )
    parser.add_argument(
        "-i",
        "--bam",
        type=str,
        required=True,
        help="minimap2 (parameters: -ax map-hifi --cs=short) aligned SAM/BAM files"
    )
    parser.add_argument(
        "--vcf",
        type=str,
        required=False,
        help="deepvariant VCF file with germline mutations",
    )
    parser.add_argument(
        "--phased_vcf",
        type=str,
        required=False,
        help="phased deepvariant VCF file",
    )
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="target chromosome",
    )
    parser.add_argument(
        "--region_list",
        type=str,
        required=False,
        help="list of target chromosomes separated by new line"
    )
    parser.add_argument(
        "--common_snps",
        type=str,
        required=False,
        help="1000G common SNPs VCF file",
    )
    parser.add_argument(
        "--panel_of_normals",
        type=str,
        required=False,
        help="panel of normal VCF file",
    )
    parser.add_argument(
        "--ploidy",
        type=str,
        default="diploid",
        required=False,
        help="haploid, diploid or polyploid"
    )
    parser.add_argument(
        "--min_mapq",
        type=int,
        default=60,
        required=False,
        help="minimum mapping quality score"
    )
    parser.add_argument(
        "--min_trim",
        type=float,
        default=0.01,
        required=False,
        help="minimum proportion of bases to be trimmed from the start and end of the read"
    )
    parser.add_argument(
        "--min_sequence_identity",
        type=float,
        default=0.99,
        required=False,
        help="minimum sequence identity threshold"
    )
    parser.add_argument(
        "--min_hq_base_proportion",
        type=float,
        default=0.5,
        required=False,
        help="minimum proportion of high quality base (BQ=93)"
    )
    parser.add_argument(
        "--min_alignment_proportion",
        type=float,
        default=0.90,
        required=False,
        help="minimum proportion of aligned CCS bases"
    )
    parser.add_argument(
        "--min_bq",
        type=int,
        default=93,
        required=False,
        help="minimum base quality score threshold"
    )
    parser.add_argument(
        "--min_ref_count",
        type=int,
        default=3,
        required=False,
        help="minimum reference allele depth at single base substitution site"
    )
    parser.add_argument(
        "--min_alt_count",
        type=int,
        default=1,
        required=False,
        help="minimum alternative allele depth at single base substitution site"
    )
    parser.add_argument(
        "--min_hap_count",
        type=int,
        default=3,
        required=False,
        help="minimum haplotype count"
    )
    parser.add_argument(
        "--mismatch_window",
        type=int,
        default=20,
        required=False,
        help="mismatch window size"
    )
    parser.add_argument(
        "--max_mismatch_count",
        type=int,
        default=0,
        required=False,
        help="maximum number of mismatches within the mismatch window"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        required=False,
        help="maximum number of threads to be used"
    )
    parser.add_argument(
        "-o",
        "--out",
        type=str,
        required=True,
        help="file to write mutations"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    if len(arguments) == 0: 
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
