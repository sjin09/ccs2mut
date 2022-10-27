import time
import math
import pysam
import bisect
import natsort
import numpy as np
import ccs2mut.util
import ccs2mut.cslib
import ccs2mut.haplib
import ccs2mut.bamlib
import ccs2mut.vcflib
import multiprocessing as mp
from dataclasses import dataclass
from typing import Dict, List, Tuple
from collections import defaultdict, Counter


def update_allelecounts(
    read,
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2bq_lst: Dict[int, Dict[int, List[int]]],
    tpos2qbase2read_lst: Dict[int, Dict[int, List[str]]]
) -> None:

    tpos2qbase = {}
    tpos = read.tstart
    qpos = read.qstart
    if read.is_primary: ## TODO
        for cstuple in read.cstuple_lst:
            state, ref, alt, ref_len, alt_len, = cstuple
            if state == 1:  # match
                for i, alt_base in enumerate(alt):
                    tpos2allelecounts[tpos + i + 1][ccs2mut.util.base2idx[alt_base]] += 1
                    tpos2qbase[tpos + i + 1] = (alt_base, read.bq_int_lst[qpos + i])
                    tpos2qbase2read_lst[tpos + i + 1][ccs2mut.util.base2idx[alt_base]].append(read.qname)
                    tpos2qbase2bq_lst[tpos + i + 1][ccs2mut.util.base2idx[alt_base]].append(read.bq_int_lst[qpos + i])
            elif state == 2:  # sub
                tpos2allelecounts[tpos + 1][ccs2mut.util.base2idx[alt]] += 1
                tpos2qbase[tpos + 1] = (alt, read.bq_int_lst[qpos])
                tpos2qbase2read_lst[tpos + 1][ccs2mut.util.base2idx[alt]].append(read.qname)
                tpos2qbase2bq_lst[tpos + 1][ccs2mut.util.base2idx[alt]].append(read.bq_int_lst[qpos])
            elif state == 3:  # insertion
                tpos2allelecounts[tpos + 1][4] += 1
                tpos2qbase2read_lst[tpos + 1][4].append(read.qname)
                tpos2qbase2bq_lst[tpos + 1][4].append(read.bq_int_lst[qpos])
            elif state == 4:  # deletion
                for j in range(len(ref[1:])):
                    tpos2allelecounts[tpos + j + 1][5] += 1
                    tpos2qbase[tpos + j + 1] = ("-", 0)
                    tpos2qbase2read_lst[tpos + j + 1][5].append(read.qname)
                    tpos2qbase2bq_lst[tpos + j + 1][5].append(0)
            tpos += ref_len
            qpos += alt_len
    else:
        tpos2qbase = ccs2mut.cslib.cs2tpos2qbase(read)
    return tpos2qbase


def get_sbs_allelecounts(
    tpos: int,
    ref: str,
    alt: str,
    tpos2allelecounts: Dict[int, np.ndarray],
    tpos2qbase2bq_lst: Dict[int, Dict[int, List[int]]]
) -> Tuple[int, int, str, float, int, int, int]:
    
    ins_count = tpos2allelecounts[tpos][4]
    del_count = tpos2allelecounts[tpos][5]
    total_count = sum(tpos2allelecounts[tpos]) - ins_count
    ref_count = tpos2allelecounts[tpos][ccs2mut.util.base2idx[ref]]
    alt_count = tpos2allelecounts[tpos][ccs2mut.util.base2idx[alt]] 
    bq = "{:.1f}".format(sum(tpos2qbase2bq_lst[tpos][ccs2mut.util.base2idx[alt]]) / float(alt_count))
    vaf = alt_count / float(total_count)
    return bq, vaf, ref_count, alt_count, ins_count, del_count, total_count


def get_somatic_substitutions(
    chrom: str,
    chrom_len: int,
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    common_snps: str,
    panel_of_normals: str,
    loci_lst: List[Tuple[str, int, int]],
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    qlen_mean: int,
    qlen_lower_limit: int,
    qlen_upper_limit: int,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    md_threshold: int,
    min_hap_count: int,
    phase: bool,
    chrom2tsbs_lst: Dict[str, List[List[Tuple[str, int, str, str, int, int, int, int, str]]]],
) -> List[Tuple[str, int, str, str, int, int, int, float, float]]:

    somatic_tsbs_lst = []
    if vcf_file.endswith(".vcf"):
        sample_snp_set = ccs2mut.vcflib.load_snp(chrom, vcf_file)
        
    if phase and ploidy == "diploid":
        (
            hpos_lst,
            hblock_lst,
            hetsnp_lst,
            hidx2hetsnp,
            hidx2hstate,
            hetsnp2bidx,
            hetsnp2hidx,
        ) = ccs2mut.vcflib.get_phased_hetsnps(phased_vcf_file, chrom, chrom_len)
 
    alignments = pysam.AlignmentFile(bam_file, "rb")
    for loci in loci_lst:
        chunkloci_lst = ccs2mut.util.chunkloci(loci)
        for chunkloci in chunkloci_lst:                
            somatic_tsbs_candidate_lst = []
            read2tpos2qbase = defaultdict(dict)
            chunk_start, chunk_end = chunkloci[1:]
            tpos2allelecounts = defaultdict(lambda: np.zeros(6)) 
            tpos2qbase2bq_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
            tpos2qbase2read_lst = defaultdict(lambda: {0: [], 1:[], 2:[], 3:[], 4:[], 5:[]})
            if vcf_file.endswith(".bgz"):
                sample_snp_set = ccs2mut.vcflib.load_bgz_snp((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), vcf_file)
            
            if common_snps.endswith(".bgz"):
                common_snp_set = ccs2mut.vcflib.load_bgz_common_snp((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), common_snps)
                
            if panel_of_normals.endswith(".bgz"):
                pon_sbs_set, _ = ccs2mut.vcflib.load_bgz_pon((chrom, chunk_start - qlen_mean, chunk_end + qlen_mean), panel_of_normals) 
                
            for line in alignments.fetch(*chunkloci):
                read = ccs2mut.bamlib.BAM(line)                
                ccs2mut.cslib.cs2tuple(read)
                ccs_somatic_tsbs_candidate_lst = [] 
                ccs_somatic_qsbs_candidate_lst = [] 
                ccs_somatic_qsbs_candidate_bq_lst = []
                ccs2mut.cslib.cs2mut(read) 
                for idx, tsbs in enumerate(read.tsbs_lst):
                    if tsbs in sample_snp_set:
                        continue
                    ccs_somatic_tsbs_candidate_lst.append(tsbs)
                    ccs_somatic_qsbs_candidate_lst.append(read.qsbs_lst[idx])
                    ccs_somatic_qsbs_candidate_bq_lst.append(read.qsbs_bq_lst[idx])
                read2tpos2qbase[read.qname] = update_allelecounts(read, tpos2allelecounts, tpos2qbase2bq_lst, tpos2qbase2read_lst)

                if not read.is_primary:
                    continue

                if read.mapq < min_mapq:
                    continue

                if read.qlen < qlen_lower_limit or read.qlen > qlen_upper_limit:
                    continue

                if read.query_alignment_proportion < min_alignment_proportion:
                    continue 

                if ccs2mut.bamlib.get_hq_base_proportion(read) < min_hq_base_proportion:
                    continue

                if ccs2mut.util.get_blast_sequence_identity(read) < min_sequence_identity:
                    continue
                
                if sum([1 for tsbs in ccs_somatic_tsbs_candidate_lst if tsbs in common_snp_set]) > 0:
                    continue
                
                if len(ccs_somatic_tsbs_candidate_lst) == 0:
                    continue
                
                for tsbs in ccs_somatic_tsbs_candidate_lst: ## TODO
                    somatic_tsbs_candidate_lst.append(tsbs)
               
                trimmed_qstart = math.floor(min_trim * read.qlen)
                trimmed_qend = math.ceil((1 - min_trim) * read.qlen)
                mismatch_lst = natsort.natsorted(ccs_somatic_tsbs_candidate_lst + read.mismatch_lst)
                mismatch_set = set(mismatch_lst)
                mpos_lst = [mismatch[1] for mismatch in mismatch_lst]
                for i, (tsbs, qsbs) in enumerate(zip(ccs_somatic_tsbs_candidate_lst, ccs_somatic_qsbs_candidate_lst)):
                    bq = ccs_somatic_qsbs_candidate_bq_lst[i]
                    if bq < min_bq:
                        continue
                    
                    if tsbs in pon_sbs_set:
                        continue

                    tpos = tsbs[1]
                    qpos = qsbs[1]
                    if qpos < trimmed_qstart:
                        continue
                    elif qpos > trimmed_qend:
                        continue

                    mismatch_start, mismatch_end = ccs2mut.util.get_mismatch_range(tpos, qpos, read.qlen, mismatch_window)
                    idx = bisect.bisect_left(mpos_lst, mismatch_start)
                    jdx = bisect.bisect_right(mpos_lst, mismatch_end)
                    if tsbs in mismatch_set:
                        mismatch_count = (jdx - idx) - 1 
                    else:
                        mismatch_count = jdx - idx 

                    if mismatch_count > max_mismatch_count:
                        continue
                    somatic_tsbs_candidate_lst.append(tsbs) 

            for tsbs in set(somatic_tsbs_candidate_lst):
                _, tpos, ref, alt = tsbs 
                if tpos <= chunk_start or tpos >= chunk_end:
                    continue

                bq, vaf, ref_count, alt_count, ins_count, del_count, total_count = get_sbs_allelecounts(tpos, ref, alt, tpos2allelecounts, tpos2qbase2bq_lst)
                if del_count != 0 or ins_count != 0:
                    continue

                if total_count > md_threshold:
                    continue

                if ref_count >= min_ref_count and alt_count >= min_alt_count:
                    if phase: # unphased single base substitutions
                        bidx2count = defaultdict(lambda: 0)
                        alt_hap2count = defaultdict(lambda: 0)
                        ref_read_lst = tpos2qbase2read_lst[tpos][ccs2mut.util.base2idx[ref]]
                        alt_read_lst = tpos2qbase2read_lst[tpos][ccs2mut.util.base2idx[alt]]
                        for alt_read in alt_read_lst:
                            bidx, alt_hap = ccs2mut.haplib.get_read_haplotype(
                                read2tpos2qbase[alt_read],
                                hpos_lst,
                                hetsnp_lst,
                                hidx2hstate,
                                hetsnp2bidx,
                                hetsnp2hidx,
                            )
                            bidx2count[bidx] += 1
                            alt_hap2count[alt_hap] += 1

                        if len(bidx2count.keys()) == 1 and len(alt_hap2count.keys()) == 1:
                            bidx = list(bidx2count.keys())[0]
                            if bidx == ".":
                                continue
                            
                            alt_hap = list(alt_hap2count.keys())[0]
                            if alt_hap == ".":
                                continue
                        else:
                            continue
                        
                        ref_hap = "1" if alt_hap == "0" else "0"
                        hap2count = ccs2mut.haplib.get_region_hap2count( 
                            ref_read_lst, read2tpos2qbase, hblock_lst[bidx], hidx2hetsnp,
                        )
                        ref_hap_count = hap2count[ref_hap]
                        alt_hap_count = hap2count[alt_hap] 
                        if ref_hap_count >= min_hap_count and alt_hap_count >= min_hap_count:
                            phase_set = hidx2hetsnp[hblock_lst[bidx][0][0]][1]
                            somatic_tsbs_lst.append((chrom, tpos, ref, alt, "PASS", bq, total_count, ref_count, alt_count, vaf, phase_set))
                    else:                         
                        somatic_tsbs_lst.append(
                            (chrom, tpos, ref, alt, "PASS", bq, total_count, ref_count, alt_count, vaf, ".")
                        )
                    
    chrom2tsbs_lst[chrom] = natsort.natsorted(somatic_tsbs_lst)
    alignments.close()


def call_somatic_substitutions(
    bam_file: str,
    vcf_file: str,
    phased_vcf_file: str,
    region: str,
    region_list: str,
    ploidy: str,
    min_mapq: int,
    min_trim: float,
    min_sequence_identity: float,
    min_hq_base_proportion: float,
    min_alignment_proportion: float,
    common_snps: str,
    panel_of_normals: str,
    min_bq: int,
    mismatch_window: int,
    max_mismatch_count: int,
    min_ref_count: int,
    min_alt_count: int,
    min_hap_count: int,
    threads: int,
    phase: bool,
    version: str,
    out_file: str,
) -> None:

    cpu_start = time.time() / 60
    ccs2mut.util.check_caller_input_exists(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        ploidy,
        common_snps,
        panel_of_normals,
        phase,
        out_file,
    )

    _, tname2tsize = ccs2mut.bamlib.get_tname2tsize(bam_file)
    chrom_lst, chrom2loci_lst = ccs2mut.util.load_loci(region, region_list, tname2tsize)
    qlen_mean, qlen_lower_limit, qlen_upper_limit, md_threshold = ccs2mut.bamlib.get_thresholds(
        bam_file, chrom_lst, tname2tsize
    )
    vcf_header = ccs2mut.vcflib.get_vcf_header(
        bam_file,
        vcf_file,
        phased_vcf_file,
        region,
        region_list,
        ploidy,         
        tname2tsize,
        common_snps,
        panel_of_normals,
        min_mapq,
        min_trim,
        min_sequence_identity,
        min_hq_base_proportion,
        min_alignment_proportion,
        min_bq,
        mismatch_window,
        max_mismatch_count,
        min_ref_count,
        min_alt_count,
        md_threshold,
        min_hap_count, 
        threads,
        phase,
        version,
        out_file,
    )
    p = mp.Pool(threads)
    manager = mp.Manager()
    chrom2tsbs_lst = manager.dict()
    get_somatic_substitutions_arg_lst = [
        (
            chrom,
            tname2tsize[chrom],
            bam_file,
            vcf_file,
            phased_vcf_file,
            common_snps,
            panel_of_normals,
            chrom2loci_lst[chrom],
            ploidy,
            min_mapq,
            min_trim,
            qlen_mean,
            qlen_lower_limit,
            qlen_upper_limit,
            min_sequence_identity,
            min_hq_base_proportion,
            min_alignment_proportion,
            min_bq,
            mismatch_window,
            max_mismatch_count,
            min_ref_count,
            min_alt_count,
            md_threshold,
            min_hap_count,
            phase,
            chrom2tsbs_lst,
        )
        for chrom in chrom_lst
    ]
    p.starmap(
        get_somatic_substitutions, get_somatic_substitutions_arg_lst,
    )
    p.close()
    p.join()
    ccs2mut.vcflib.dump_ccs2mut_sbs(chrom_lst, chrom2tsbs_lst, phase, vcf_header, out_file)  
    print("finished calling and returning substitutions with {} threads".format(threads))
    cpu_end = time.time() / 60
    duration = cpu_end - cpu_start
    print("single molecule somatic mutation detection took {} minutes".format(duration))
