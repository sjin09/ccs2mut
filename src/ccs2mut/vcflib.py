import tabix
import cyvcf2
import natsort
import numpy as np 
import ccs2mut.util
import ccs2mut.bamlib
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple


class VCF:
    def __init__(self, line):
        arr = line.strip().split()
        self.chrom = arr[0]
        self.pos = int(arr[1])
        self.id = arr[2]
        self.ref = arr[3]
        self.alt_lst = arr[4].split(",")
        self.qual = arr[5]
        self.qual = float(self.qual) if self.qual != "." else self.qual
        self.is_pass = True if arr[6] == "PASS" else False
        self.info = arr[7]
        self.format_lst = arr[8].split(":")
        self.sample_format_lst = arr[9].split(":")
        hsh = {i:j for i,j in zip(self.format_lst, self.sample_format_lst)}
        if "GT" in hsh: self.sample_gt = hsh["GT"] 
        if "PS" in hsh: self.sample_phase_set = hsh["PS"] 
        if "AD" in hsh: 
            arr = hsh["AD"].split(",")
            self.ref_count = arr[0]
            self.alt_count_arr = arr[1:]

        self.is_snp = False
        self.is_dbs = False
        self.is_indel = False
        if len(self.alt_lst) == 1:  # bi-allelic
            self.is_biallelic = True
            self.alt = self.alt_lst[0]
            if len(self.ref) == 1 and len(self.alt) == 1:  # snp
                self.is_snp = True
            elif len(self.ref) == len(self.alt) == 2:
                self.is_dbs = True
            elif len(self.ref) > len(self.alt): # del
                self.is_indel = True
            elif len(self.ref) < len(self.alt): # ins
                self.is_indel = True
        else:
            self.is_biallelic = False


def get_sample(vcf_file: str) -> str:
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#CHROM"):
                sample = line.strip().split()[-1]
                return sample
    elif vcf_file.endswith(".bgz"):
        v = cyvcf2.VCF(vcf_file)
        sample = v.samples[0]
        return sample


def load_snp(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    snp_set = set()
    for line in open(vcf_file):
        if line.startswith("#"):
            continue
        v = VCF(line)
        if v.chrom != chrom:
            continue
        if v.is_pass:
            if v.is_biallelic:
                if v.is_snp:
                    snp_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.chrom, v.pos, v.ref, alt))
    return snp_set 
      
                     
def load_bgz_snp(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:

    snp_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        if v.is_pass:
            if v.is_biallelic:
                if v.is_snp:
                    snp_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                for alt in v.alt_lst:
                    if len(v.ref) == 1 and len(alt) == 1:
                        snp_set.add((v.chrom, v.pos, v.ref, alt))
    return snp_set


def load_pon(
    chrom: str,
    vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:

    state = 0
    sbs_set = set() 
    dbs_set = set()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom: 
                if state == 0:
                    continue
                else:
                    break
            else:
                state = 1
                if v.is_pass and v.is_biallelic:
                    if v.is_snp:
                        sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
                    else:
                        if len(v.ref) == len(v.alt) == 2:
                            dbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set, dbs_set


def load_bgz_pon(
    loci: Tuple[str, int, int],
    vcf_file: str
) -> Tuple[Set[Tuple[str, int, str, str]], Set[Tuple[str, int, str, str]]]:
   
    dbs_set = set() 
    sbs_set = set() 
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for record in records:
        v = VCF("\t".join(record))            
        if v.is_pass and v.is_biallelic:
            if v.is_snp:
                sbs_set.add((v.chrom, v.pos, v.ref, v.alt))
            else:
                if len(v.ref) == len(v.alt) == 2:
                    dbs_set.add((v.chrom, v.pos, v.ref, v.alt))
    return sbs_set, dbs_set


def load_common_snp(
    chrom: str,
    vcf_file: str
) -> Set[Tuple[str, int, str, str]]:
    
    snp_set = set()
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file).readlines():
            if line.startswith("#"):
                continue
            arr = line.strip().split()
            if chrom != arr[0]:
                continue
            if arr[6] == "PASS":
                alt_lst = arr[4].split(",") 
                if len(alt_lst) == 1:            
                    ref = arr[3]
                    alt = alt_lst[0]
                    if len(ref) == 1 and len(alt) == 1:
                        snp_set.add((chrom, int(arr[1]), ref, alt))
    return snp_set


def load_bgz_common_snp(
    loci: Tuple[str, int, int],
    vcf_file: str,
) -> Set[Tuple[str, int, str, str]]:
  
    snp_set = set()
    tb = tabix.open(vcf_file)
    records = tb.query(*loci)
    for arr in records:
        alt_lst = arr[4].split(",") 
        if arr[6] == "PASS":
            if len(alt_lst) == 1:
                ref = arr[3]
                alt = alt_lst[0]
                if len(ref) == 1 and len(alt) == 1:
                    snp_set.add((arr[0], int(arr[1]), ref, alt))
    return snp_set


def load_hetsnps(
    vcf_file: str,
    chrom: str,
    chrom_len: int,
):

    hetsnp_lst = []
    hidx2hetsnp = {} 
    hetsnp2hidx = {} 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for record in records:
            v = VCF("\t".join(record))            
            if v.is_snp and v.is_pass and v.is_biallelic and (v.sample_gt == "0/1" or v.sample_gt == "1/0"):
                hetsnp_lst.append((v.chrom, v.pos, v.ref, v.alt))
    for hidx, hetsnp in enumerate(hetsnp_lst):
        hidx2hetsnp[hidx] = hetsnp
        hetsnp2hidx[hetsnp] = hidx
    return hetsnp_lst, hidx2hetsnp, hetsnp2hidx


def load_phased_hetsnp_index(
    vcf_file: str, 
    chrom: str,
    chrom_len: int,
) -> Dict[str, List[List[Tuple[int, str]]]]:

    hidx = 0  
    state = 0
    hidx2hetsnp = {}
    hetsnp2hidx = {}
    hblock_lst = defaultdict(list)
    phase_set2hblock = defaultdict(list) 
    if vcf_file.endswith(".vcf"):
        for line in open(vcf_file):
            if line.startswith("#"):
                continue
            v = VCF(line)
            if chrom != v.chrom:
                if state == 0:
                    continue
                elif state == 1:
                    break
            elif chrom == v.chrom: 
                if state == 0:
                    state = 1
                if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                    hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                    hidx2hetsnp[hidx] = hetsnp
                    hetsnp2hidx[hetsnp] = hidx
                    hstate = v.sample_gt.split("|")[0]                    
                    phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
                hidx += 1 
    elif vcf_file.endswith(".bgz"):
        tb = tabix.open(vcf_file)
        records = tb.query(chrom, 0, chrom_len)
        for hidx, record in enumerate(records):
            v = VCF("\t".join(record))           
            if v.sample_gt == "0|1" or v.sample_gt == "1|0":
                hetsnp = (v.chrom, v.pos, v.ref, v.alt)
                hidx2hetsnp[hidx] = hetsnp
                hetsnp2hidx[hetsnp] = hidx                
                hstate = v.sample_gt.split("|")[0]
                phase_set2hblock[v.sample_phase_set].append((hidx, hstate))
    hblock_lst = [hblock for hblock in phase_set2hblock.values()]
    return hblock_lst, hidx2hetsnp, hetsnp2hidx


def get_phased_hetsnps(
    vcf_file: str, 
    chrom: str,
    chrom_len: int, 
):
    hetsnp2bidx = {}
    hidx2hstate = {} 
    filtered_hblock_lst = []
    hblock_lst, hidx2hetsnp, hetsnp2hidx = load_phased_hetsnp_index(vcf_file, chrom, chrom_len)
    for hblock in hblock_lst:
        ipos = 0
        filtered_hblock = []
        for (hidx, hstate) in hblock:
            jpos = hidx2hetsnp[hidx][1]
            if jpos - ipos > 1:
                filtered_hblock.append((hidx, hstate))
            ipos = jpos
        if len(filtered_hblock) > 1:
            filtered_hblock_lst.append(filtered_hblock)
   
    phased_hetsnp_lst = []
    for bidx, hblock in enumerate(filtered_hblock_lst):
        for (hidx, hstate) in hblock:
            hidx2hstate[hidx] = hstate 
            hetsnp = hidx2hetsnp[hidx]
            hetsnp2bidx[hetsnp] = bidx 
            phased_hetsnp_lst.append(hetsnp)
    phased_hetsnp_lst = natsort.natsorted(phased_hetsnp_lst) 
    hpos_lst = [hetsnp[1] for hetsnp in phased_hetsnp_lst]
    return hpos_lst, filtered_hblock_lst, phased_hetsnp_lst, hidx2hetsnp, hidx2hstate, hetsnp2bidx, hetsnp2hidx


def dump_ccs2mut_sbs(
    chrom_lst: List[str],
    chrom2tsbs_lst: Dict[str, List[Tuple[str, int, str, str, str, int, int, int, float]]],
    phase: bool,
    header: str,
    out_file: str
) -> None:

    if not out_file.endswith(".vcf"):
        print("VCF file must have .vcf suffix")
        ccs2mut.util.exit()
        
    o = open(out_file, "w")
    o.write("{}\n".format(header))
    for chrom in chrom_lst:
        for (chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf, phase_set) in chrom2tsbs_lst[chrom]:
            if phase:
                if alt_count == 1:
                    o.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF:PS\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}:{}\n".format(
                            chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf, phase_set
                        )
                    )
            else:
                if alt_count == 1:
                    o.write(
                        "{}\t{}\t.\t{}\t{}\t.\t{}\t.\tGT:BQ:DP:AD:VAF\t./.:{}:{:0.0f}:{:0.0f},{:0.0f}:{:.2f}\n".format(
                            chrom, pos, ref, alt, annot, bq, total_count, ref_count, alt_count, vaf
                        )
                    )
    o.close()

