__all__ = ['simulate_read', 'format_simulated_reads']

from math import log, exp
from numpy import random

from .contamination_est import transition, emission

STATES = ("5ss", "ds", "ss", "3ss")
N_STATES = len(STATES)


def update_seq(sequence, base, forward):
    return sequence + base if forward else base + sequence



def update_state(previous_state, p, d, s, o2, d3):
    p5ss = transition(previous_state, "5ss", p, d, s, o2, d3, 1)
    pds = transition(previous_state, "ds", p, d, s, o2, d3, 1)
    pss = transition(previous_state, "ss", p, d, s, o2, d3, 1)
    p3ss = transition(previous_state, "3ss", p, d, s, o2, d3, 1)

    p = random.uniform(0,1)

    if p < p5ss:
        new_state = "5ss"
    elif p < p5ss + pds:
        new_state = "ds"
    elif p < p5ss + pds + pss:
        new_state = "ss"
    else:
        new_state = "3ss"

    return new_state



def sample_ref(GC, forward):
    p = random.uniform(0,1)
    if p < GC / 2:
        ref = "G"
        if not forward:
            informative = True
        else:
            informative = False
    elif p < GC:
        ref = "C"
        if forward:
            informative = True
        else:
            informative = False
    elif p < GC + (1 - GC)/2:
        ref = "T"
        informative = False
    else:
        ref = "A"
        informative = False
    return ref, informative



def sample_seq(state, ref, e, rss, rds, forward, informative):
    if informative:
        p_deam = emission(state, "deam", e, rss, rds, informative)
        p_error = emission(state, "error", e, rss, rds, informative)
        p_match = 1 - p_deam - p_error
    else:
        p_deam = 0
        p_error = emission(state, "error", e, rss, rds, informative)
        p_match = 1 - p_error

    p = random.uniform(0,1)
    if p < p_deam:
        if forward:
            seq = "T"
        else:
            seq = "A"
    elif p < p_deam + p_error:
        if forward:
            seq = random.choice([i for i in ["A","C","G","T"] if (i!=ref and ref+i!="CT")])
        else:
            seq = random.choice([i for i in ["A","C","G","T"] if (i!=ref and ref+i!="GA")])
    else:
        seq = ref
    return seq



def call_MD(reference, sequence):
    match = 0
    MD = "MD:Z:"

    for k in range(len(reference)-1):
        if reference[k] == sequence[k]:
            match += 1
        else:
            MD += str(match) + reference[k]
            match = 0

    if reference[len(reference)-1] == sequence[len(reference)-1]:
        match += 1
        MD += str(match)
    else:
        MD += str(match) + reference[len(reference)-1] + str(0)

    return MD



def simulate_read(length, e, rss, p, s, d, rds, contam, o, o2, GC, forward=True, endogenous=True):
    reference = str()
    sequence = str()
    states = "ST:Z:"

    if forward:
        seqiter = range(2,length+1)
    else:
        seqiter = reversed(range(1,length))

    state = "5ss" if (random.uniform(0,1) <= o) else "ds"
    ref_base, informative = sample_ref(GC, forward)
    seq_base = sample_seq(state, ref_base, e, rss, rds, forward, informative)

    reference = update_seq(reference, ref_base, forward)
    sequence = update_seq(sequence, seq_base, forward)
    states = update_seq(states, state[0], forward)

    for j in seqiter:
        if forward:
            d3 = length - j
        else:
            d3 = j - 1

        state = update_state(state, p, d, s, o2, d3)
        ref_base, informative = sample_ref(GC, forward)
        seq_base = sample_seq(state, ref_base, e, rss, rds, forward, informative)

        reference = update_seq(reference, ref_base, forward)
        sequence = update_seq(sequence, seq_base, forward)
        states = update_seq(states, state[0], forward)

    MD = call_MD(reference, sequence)
    
    return sequence, MD, states



def format_simulated_reads(parameters, length, minlength, GC, N):
    for i in range(1,N+1):
        forward = True
        endogenous = True
        seq_length = random.geometric(p=1./(length-minlength)) + minlength
        sequence, MD, states = simulate_read(seq_length, *parameters, GC, forward, endogenous)
        yield (
            "Read"+str(i),
            "0" if forward else "16",
            "1",#RNAME
            "1",#POS
            "30",#MAPQ
            str(seq_length)+"M",
            "*",
            "0",
            "0",
            sequence,
            "]"*seq_length,#QUAL
            MD,
            states
        )

