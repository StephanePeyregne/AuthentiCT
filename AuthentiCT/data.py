__all__ = ['InformativeFragment']

import re
from functools import lru_cache

from . import CIGAR_REGEX

MD_REGEX = re.compile(r"(\d+|[AaCcGgTtNn]|\^[AaCcTtGgNn]+)", flags=re.ASCII)
# Value for MD field must match r'\d+(([ACGT]|\^[ACGT]+)\d+)*',
# but `re` can't capture repeated groups (the future `regex` can)


@lru_cache(maxsize=512)
def ascii_to_bq(character):
    return ord(character) - 33
#    return pow(10,-(ord(character)-33) // 10)

def create_ref(md_field, cigar, read):
    md_list = MD_REGEX.findall(md_field)

    md_counter = 0
    ref_seq = ''
    aligned_read = ''
    alignment = {'ref': '', 'mol': ''}
    for e in md_list:
        try:
            e = int(e)
            ref_seq += read[md_counter:md_counter + e]
            alignment['ref'] += '.' * e
            alignment['mol'] += '.' * e
            aligned_read += read[md_counter:md_counter + e]
            md_counter += e
        except ValueError:
            if e.startswith('^'):
                e = e.lstrip('^')
                ref_seq += e
                alignment['ref'] += e.upper()
                alignment['mol'] += '-' * len(e)
                aligned_read += '-' * len(e)
                continue
            else:
                ref_seq += e
                alignment['ref'] += e.upper()
                alignment['mol'] += '.'
                aligned_read += read[md_counter]
                md_counter += 1

    if 'I' in cigar or 'S' in cigar:
        # find insertions and clips in cigar
        # update the ref and the alignment using indel and clip info
        ref_seq = ''
        aligned_read = ''
        readcounter = 0
        aligncounter = 0
        cigarcomponents = CIGAR_REGEX.findall(cigar)
        for val, op in cigarcomponents:
            val = int(val)
            if op == 'D':
                ref_seq += alignment['ref'][aligncounter:aligncounter + val]
                aligned_read += '-' * val
                aligncounter += val
            elif op in 'IS':
                ref_seq += '-' * val
                aligned_read += read[readcounter:readcounter + val]
                readcounter += val
            elif op == 'M':
                for c in range(val):
                    if alignment['ref'][aligncounter] == '.':
                        ref_seq += read[readcounter]
                    else:
                        ref_seq += alignment['ref'][aligncounter]
                    aligned_read += read[readcounter]
                    aligncounter += 1
                    readcounter += 1
            # No other Cigar operations, checked already in read_input() (cf. input.py)
            # Done earlier as we do not want to return an empty object here
    # end cigar parsing
    return aligned_read, ref_seq


class InformativeFragment:
    """Class defining a DNA fragment by:
    - name
    - chromosome
    - starting position
    - sequenced strand
    - sequence length
    - positions (in fragment and in reference) aligning to Cs or Gs on the reference
    - alleles at these positions
    - allele at the informative position

    """

    def __init__(self, qname, flag, rname, pos, cigar, seq, qual, md_value, min_bq, sites=None):

        self.name = qname
        self.rname = rname
        self.length = len(seq)
        self.sequence = ''
        self.baseq = []
        self.reference = ''
        self.ingenome = []
        self.inread = []

        # Define which strand was sequenced
        try:
            self.reverse = int(flag) & 16
        except ValueError:
            self.reverse = 'r' in flag

        aligned_read, ref_seq = create_ref(md_value, cigar, seq)
        genome_pos = int(pos) - 1
        read_pos = 0
        for a, p in zip(aligned_read, ref_seq):
            if p != '-':
                genome_pos += 1
            if a != '-':
                read_pos += 1
                if ascii_to_bq(qual[read_pos - 1]) >= min_bq:
                    try:
                        sites[(rname, str(genome_pos))][0]+=1
                    except (KeyError, NameError, TypeError) as error:
                        pass
                    if ((self.reverse and p == 'G') or (self.forward and p == 'C')):
                        self.sequence += a
                        self.baseq.append(ascii_to_bq(qual[read_pos - 1]))
                        self.reference += p
                        self.ingenome.append(genome_pos)
                        self.inread.append(read_pos)


    @property
    def forward(self):
        return not self.reverse

    @property
    def strand(self):
        return 'R' if self.reverse else 'F'


