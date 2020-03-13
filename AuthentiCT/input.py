__all__ = ['read_sites','read_parameters','read_input','count_terminal_CtoT','deamination_patterns', 'distance_btw_deam', 'count_errors']#,'get_covered_infosites']

import re
import sys

from . import CIGAR_REGEX
from .data import *

LIMITED_CIGAR = re.compile(r"[^\dMISD]", flags=re.ASCII)

def read_parameters(filename):
    parameters={}
    expression = re.compile("^(e|rss|lo|lss|lds|contam|rds|o|o2)[\t ]+(\d+\.\d+)")
    for line in filename:
        info = expression.search(line)
        if info != None:
            parameters[info.group(1)]=float(info.group(2))
    try:
        return [parameters["e"],parameters["rss"],parameters["lo"],parameters["lss"],parameters["lds"],parameters["rds"],parameters["contam"],parameters["o"],parameters["o2"]]
    except (KeyError, ValueError):
        print("\nPlease provide all parameters in the following tab-separated format:\n"
               "e           0.003170\n"
               "rss         0.602679\n"
               "lo          0.261435\n"
               "lss         0.142886\n"
               "lds         0.008414\n"
               "rds         0.035792\n"
               "contam      0.001000\n"
               "o           0.637999\n"
               "o2          0.477368")
        sys.exit(1)

def read_sites(filename):
    sites={} #Should define a list for each type of sites
    for pos in filename:
        pos = pos.rstrip().split('\t')
        try:
            if pos[1].isdigit() and pos[2].isalpha() and pos[3].isalpha():
            #store positions in a dictionary with keys (chromosome, position)
            #and values: number of overlapping fragments (0 for now) and a tuple (ancestral allele, derived allele)
                sites[(pos[0],pos[1])] = [0, (pos[2].upper(), pos[3].upper())]
            else:
                print("positions should be integers (1-based) while ancestral and derided states should be alphabetic characters.")
                sys.exit(1)
        except IndexError:
            print("Please provide the informative positions in a tab-separated file following this format: chr  pos     anc     der     type(s).")
            sys.exit(1)
    return sites


def get_coordinates(sites):
    coordinates={}
    #Make lists of informative sites, one for each chromosome
    for key in sites.keys():
        if key[0] not in coordinates:
            coordinates[key[0]]=[int(key[1])]
        else:
            coordinates[key[0]].append(int(key[1]))
    #Sort these lists
    for chromosome in coordinates:
        coordinates[chromosome].sort()
    return coordinates


def read_cigar(cigar, pos, seqlen):
    if 'I' not in cigar and 'S' not in cigar and 'D' not in cigar:
        return pos + seqlen - 1
    else:
        for n, op in CIGAR_REGEX.findall(cigar):
            if op not in 'IS':
                pos += int(n)
        return pos - 1


def read_input(min_mapq, minlength, min_bq, sites=None, filter_reads=False, fh=sys.stdin):
    if sites is not None:
        coordinates = get_coordinates(sites)
    else:
        coordinates = None
    # Assuming the aligned bam file is sorted, go through the
    # reads and get the ones overlapping the informative sites
    for line in fh:
        # Skip the header
        if line.startswith('@'):
            continue

        qname, flag, rname, pos, mapq, cigar, _, _, _, seq, qual, *opt = \
            line.rstrip().split('\t')

        # Check if the BAM file is sorted
        # Only necessary if working with informative positions

        if int(mapq) < min_mapq:
            continue
        if len(seq) < minlength:
            continue

        # Check cigar operations
        if LIMITED_CIGAR.search(cigar):
            print(qname, '-- CIGAR operation other than M I S D; excluded', file=sys.stderr)
            continue
        end_position = read_cigar(cigar, int(pos), len(seq))

        if filter_reads == True:
            if coordinates is not None:
                if rname not in coordinates or coordinates[rname] == []:
                    continue
                while coordinates[rname][0] < int(pos):
                    del coordinates[rname][0]
                    if coordinates[rname] == []:
                        break
                if coordinates == []:
                    break
                if coordinates[rname] == []:  # Skip sequences aligning to this chromosome
                    continue
                if coordinates[rname][0] > end_position:
                    continue

        opt_fields = dict(e.split(':', 1) for e in opt)

        # Recreate ref seq first from MD field
        try:
            _, md_value = opt_fields['MD'].split(':')
        except (KeyError, ValueError):
            print(qname, '-- no MD field found; excluded', file=sys.stderr)
            continue

        yield InformativeFragment(qname, flag, rname, pos, cigar, seq, qual, md_value, min_bq, sites)



def count_terminal_CtoT(reads):
    counts={"deam53" : 0 , "deam5" : 0, "deam3" : 0, "nodeam" : 0}
    for read in reads:
        deam5 = False
        deam3 = False
        if (1 not in read.inread) or (read.length not in read.inread):
            continue
        if read.forward:
            if read.reference[0]!="C" or read.reference[-1]!="C":
                continue
            if read.sequence[0]=="T": deam5 = True
            if read.sequence[-1]=="T": deam3 = True
        else:
            if read.reference[0]!="G" or read.reference[-1]!="G":
                continue
            if read.sequence[0]=="A": deam3 = True
            if read.sequence[-1]=="A": deam5 = True
        if deam5 == True and deam3 == True:
            counts["deam53"]+=1
        elif deam5 == True:
            counts["deam5"]+=1
        elif deam3 == True:
            counts["deam3"]+=1
        else:
            counts["nodeam"]+=1
    return counts       


def count_errors(reads):
    error = 0
    no_error = 0
    for read in reads:
        if read.forward:
            informative = "C"
            deam = "T"
        else:
            informative = "G"
            deam = "A"
        for i in range(len(read.reference)):
            if read.reference[i]==informative:
                if read.sequence[i]==deam:
                    continue
                if read.sequence[i]!=informative:
                    error+=1
                else:
                    no_error+=1
    return error, error + no_error


def deamination_patterns(reads):
    counts={}
    for read in reads:
        length, forward, reverse = read.length, read.forward, read.reverse
        if forward:
            seqiter = range(len(read.sequence))
        else:
            seqiter = reversed(range(len(read.sequence)))
        for i in seqiter:
            C = 0
            T = 0
            if forward:
                d5 = read.inread[i]
                d3 = -(length + 1 - read.inread[i])
                if read.reference[i]=="C" and read.sequence[i]=="T":
                    T = 1
                elif read.reference[i]=="C" and read.sequence[i]=="C":
                    C = 1
            else:
                d5 = length + 1 - read.inread[i]
                d3 = -read.inread[i]
                if read.reference[i]=="G" and read.sequence[i]=="A":
                    T = 1
                elif read.reference[i]=="G" and read.sequence[i]=="G":
                    C = 1
            for j in d5,d3:
                if j>=-50 and j<=50: 
                    try:
                        counts[j][0] += C
                        counts[j][1] += T
                    except KeyError:
                        counts[j]=[C,T]
    for d in sorted(counts.keys()):
        low, up = calcBin(counts[d][1], counts[d][0] + counts[d][1])
        yield(d,
              float(counts[d][1]) / (counts[d][0] + counts[d][1]),
              counts[d][1],
              counts[d][0] + counts[d][1],
              low,
              up)

def distance_btw_deam(reads, mask=5):
    counts = {}
    length_seq = 0
    deam = 0
    nodeam = 0
    pairs = 0
    total = 0
    for read in reads:
        length, forward, reverse = read.length, read.forward, read.reverse
        masked_seq = [(base, pos) for base, pos in zip(read.sequence, read.inread) if (pos > mask and pos < length - mask)]
        bases = [v[0] for k,v in enumerate(masked_seq)]
        if forward:
            if bases.count("T")!=2:
                continue
            else:
                deam += bases[bases.index("T")+1:].count("T")
                nodeam += bases[bases.index("T")+1:].count("C")
                positions = [v[1] for k,v in enumerate(masked_seq) if v[0]=="T"]
                length_seq += length - mask - bases.index("T")
        else:
            if bases.count("A")!=2:
                continue
            else:
                deam += bases[bases.index("A")+1:].count("A")
                nodeam += bases[bases.index("A")+1:].count("G")
                positions = [v[1] for k,v in enumerate(masked_seq) if v[0]=="A"]
                length_seq += length - mask - bases.index("A")

        for i,j in zip(positions, positions[1:]):
            pairs += 1
            try:
                counts[j-i] += 1
            except KeyError:
                counts[j-i] = 1

#        length_seq += length - 2 * mask
        total += 1
    try:
        avg_length = int(float(length_seq) / total)
    except ZeroDivisionError:
        print("No sequence passed the filtering, you might consider decreasing the number of masked positions or use more sequences") 
        sys.exit(1)

    deam_rate = float(deam) / (deam + nodeam)
    for d in sorted(counts.keys()):
        if d <= avg_length:
            yield(d,
                  float(counts[d])/pairs,
                  (1 - deam_rate) ** (d-1) * deam_rate / (1 - (1 - deam_rate) ** avg_length),
                  counts[d])


def binP(N, p, x1, x2):
    p = float(p)
    q = p/(1-p)
    k = 0.0
    v = 1.0
    s = 0.0
    tot = 0.0

    while(k<=N):
        tot += v
        if(k >= x1 and k <= x2):
            s += v
        if(tot > 10**30):
            s = s/10**30
            tot = tot/10**30
            v = v/10**30
        k += 1
        v = v*q*(N+1-k)/k
    return s/tot

def calcBin(vx, vN, vCL = 95):
    '''
    Calculate the exact confidence interval for a binomial proportion
    '''
    vx = float(vx)
    vN = float(vN)
    #Set the confidence bounds
    vTU = (100 - float(vCL))/2
    vTL = vTU

    vP = vx/vN
    if(vx==0):
        dl = 0.0
    else:
        v = vP/2
        vsL = 0
        vsH = vP
        p = vTL/100

        while((vsH-vsL) > 10**-5):
            if(binP(vN, v, vx, vN) > p):
                vsH = v
                v = (vsL+v)/2
            else:
                vsL = v
                v = (v+vsH)/2
        dl = v

    if(vx==vN):
        ul = 1.0
    else:
        v = (1+vP)/2
        vsL =vP
        vsH = 1
        p = vTU/100
        while((vsH-vsL) > 10**-5):
            if(binP(vN, v, 0, vx) < p):
                vsH = v
                v = (vsL+v)/2
            else:
                vsL = v
                v = (v+vsH)/2
        ul = v
    return dl, ul
