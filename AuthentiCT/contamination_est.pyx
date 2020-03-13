#cython: language_level=3
__all__ = ['compute_loglikelihood', 'format_loglikelihood','transition','emission']

from math import log, exp
from numpy import logaddexp
import numpy as np
from functools import lru_cache

from libcpp.vector cimport vector

STATES = ("5ss", "ds", "ss", "3ss")
N_STATES = len(STATES)

memoize = lru_cache(maxsize=None)

def transition(prev_state, state, p, d, s, o2, d3, d5):
    if prev_state == "5ss":
        if state == "5ss":
            return (1 - p) ** d5
        elif state == "ds":
            return 1 - (1 - p) ** d5
        else:
            return 0

    if prev_state == "ds":
        if state == "5ss":
            return 0
        elif state == "ds":
            return ((1 - d) ** d5) * (1 - o2 * (1 - p) ** d3)
        elif state == "ss":
            return (1- (1 - d) ** d5) * (1 - o2 * (1 - p) ** d3)
        elif state == "3ss":
            return  o2 * (1 - p) ** d3

    if prev_state == "ss":
        if state == "ds":
            return 1 - (1 - s) ** d5
        elif state == "ss":
            return (1 - s) ** d5
        else:
            return 0

    if prev_state == "3ss":
        if state == "3ss":
            return 1
        else:
            return 0


# Transition matrix:
#      5ss ds ss 3ss
# 5ss (             )
# ds  (             )
# ss  (             )
# 3ss (             )

@memoize
def tmat(p, d, s, o2, d3, d5):
#    M = np.zeros((N_STATES, N_STATES))
#    for i, s1 in enumerate(STATES):
#        for j, s2 in enumerate(STATES):
#            M[i, j] = transition(s1, s2, p, d, s, o2, d3, d5)
#    return M
    pd5 = (1 - p) ** d5
    dd5 = (1 - d) ** d5
    o2_pd3 = o2 * (1 - p) ** d3
    sd5 = (1 - s) ** d5
    M = np.array([[pd5, 1 - pd5, 0., 0.],
                  [0., dd5 * (1 - o2_pd3), (1 - dd5) * (1 - o2_pd3), o2_pd3],
                  [0., 1 - sd5, sd5, 0.],
                  [0., 0., 0., 1.]])
    return M

def emission(state, obs, e, rss, rds, informative):
    if informative:
        if state == "ds" and obs == "match":
            return ((1 - e) * (1 - rds) + (e / 3) * rds)
        elif state != "ds" and obs == "match":
            return (1 - e) * (1 - rss) + (e / 3) * rss

        elif state == "ds" and obs == "deam":
            return  (1 - e) * rds + (e / 3) * (1 - rds)

        elif state != "ds" and obs == "deam":
            return  (1 - e) * rss + (e / 3) * (1 - rss)

        elif obs == "error":
            return 2 * e / 3
        else:
            raise ValueError("unexpected observation")
    else:
        if obs == "match":
            return 1 - e
        elif obs == 'error':
            return e
        else:
            raise ValueError("unexpected observation")


def evec(obs, e, rss, rds, informative):
    return np.array([emission(state, obs, e, rss, rds, informative) for state in STATES])


cpdef read_loglikelihood(object read, double e, double rss, double p, double s, double d, double rds, double contam, double o, double o2, vector[double] match, vector[double] deam, vector[double] error, bint maximize=True):
    cdef double likelihood_endo = float()
    cdef double likelihood_cont = float()

    length, forward, reverse = read.length, read.forward, read.reverse

    if forward:
        seqiter = range(len(read.sequence))
    else:
        seqiter = reversed(range(len(read.sequence)))

    cdef int previous = 1 if forward else length

    cdef vector[double] fwd_5ss_scaled = [o]
    cdef vector[double] fwd_ds_scaled = [1-o]
    cdef vector[double] fwd_ss_scaled = [0.]
    cdef vector[double] fwd_3ss_scaled = [0.]
    cdef vector[double] fwd_scaling = [1.]
    #fwd = [np.array([o, 1-o, 0, 0])] #5ss, ds, ss, 3ss


    last = -1
    i = -1
    cdef int j
    cdef int d5
    cdef int d3
    cdef double p_d5
    cdef double p_d3
    cdef double d_d5
    cdef double s_d5
    cdef double subterm01
    cdef double subterm02
    cdef double tmp_5ss
    cdef double tmp_ds
    cdef double tmp_ss
    cdef double tmp_3ss

    for j in seqiter:
        if i== -1:
            first = j
        i += 1
        last = j
        if forward:
            d5 = read.inread[j] - previous
            d3 = length - read.inread[j]
            previous = read.inread[j]
        else:
            d5 = previous - read.inread[j]
            d3 = read.inread[j] - 1
            previous = read.inread[j]

        p_d5 = (1 - p) ** d5
        p_d3 = (1 - p) ** d3
        d_d5 = (1 - d) ** d5
        s_d5 = (1 - s) ** d5
        subterm01 = 1 - p_d3 * o2#(1 - p_d3) * o2 + (1 - o2)
        subterm02 = p_d3 * o2

        if read.sequence[j] == read.reference[j]:

            tmp_5ss = match[0] * fwd_5ss_scaled[i] * p_d5
            tmp_ds = match[1] * (fwd_5ss_scaled[i] * (1 - p_d5)
                         + fwd_ds_scaled[i] * d_d5 * subterm01
                         + fwd_ss_scaled[i] * (1 - s_d5))
            tmp_ss = match[2] * (fwd_ss_scaled[i] * s_d5
                         + fwd_ds_scaled[i] * (1 - d_d5) * subterm01)
            tmp_3ss = match[3] * (fwd_3ss_scaled[i]
                         + fwd_ds_scaled[i] * subterm02)

            likelihood_cont += log(1 - e)
        elif (forward and read.sequence[j] == "T") or (reverse and read.sequence[j] == "A"):

            tmp_5ss = deam[0] * fwd_5ss_scaled[i] * p_d5
            tmp_ds = deam[1] * (fwd_5ss_scaled[i] * (1 - p_d5)
                         + fwd_ds_scaled[i] * d_d5 * subterm01
                         + fwd_ss_scaled[i] * (1 - s_d5))
            tmp_ss = deam[2] * (fwd_ss_scaled[i] * s_d5
                         + fwd_ds_scaled[i] * (1 - d_d5) * subterm01)
            tmp_3ss = deam[3] * (fwd_3ss_scaled[i]
                         + fwd_ds_scaled[i] * subterm02)

            likelihood_cont += log(e / 3)
        else:

            tmp_5ss = error[0] * fwd_5ss_scaled[i] * p_d5
            tmp_ds = error[1] * (fwd_5ss_scaled[i] * (1 - p_d5)
                         + fwd_ds_scaled[i] * d_d5 * subterm01
                         + fwd_ss_scaled[i] * (1 - s_d5))
            tmp_ss = error[2] * (fwd_ss_scaled[i] * s_d5
                         + fwd_ds_scaled[i] * (1 - d_d5) * subterm01)
            tmp_3ss = error[3] * (fwd_3ss_scaled[i]
                         + fwd_ds_scaled[i] * subterm02)

            likelihood_cont += log(2 * e / 3)


        #transition_mat = tmat(p, d, s, o2, d3, d5)
        #tmp_prob = np.matmul(fwd[i], transition_mat) * emission_vector


        fwd_scaling.push_back(tmp_5ss + tmp_ds + tmp_ss + tmp_3ss)

        fwd_5ss_scaled.push_back(tmp_5ss / fwd_scaling[i+1])
        fwd_ds_scaled.push_back(tmp_ds / fwd_scaling[i+1])
        fwd_ss_scaled.push_back(tmp_ss / fwd_scaling[i+1])
        fwd_3ss_scaled.push_back(tmp_3ss / fwd_scaling[i+1])

        likelihood_endo += log(fwd_scaling[i + 1])

    if last == -1:
        return 0,0,0,0,0,0,0

#    if maximize == True:
#        return likelihood_endo, likelihood_cont, 0, 0, 0, 0

    d_d3 = (1 - d) ** d3
    tmp_5ss = 0
    tmp_ds = (fwd_5ss_scaled[i+1] #5ss
              + fwd_ss_scaled[i+1] #ss
              + fwd_ds_scaled[i+1] * d_d3) #ds
    tmp_ss = 0
    tmp_3ss = (fwd_ds_scaled[i+1] * (1 - d_d3) #ds
               + fwd_3ss_scaled[i+1]) #3ss

    fwd_scaling.push_back(tmp_5ss + tmp_ds + tmp_ss + tmp_3ss)
    fwd_5ss_scaled.push_back(tmp_5ss / fwd_scaling[i + 2])
    fwd_ds_scaled.push_back(tmp_ds / fwd_scaling[i + 2])
    fwd_ss_scaled.push_back(tmp_ss / fwd_scaling[i + 2])
    fwd_3ss_scaled.push_back(tmp_3ss / fwd_scaling[i + 2])

    likelihood_endo += log(fwd_scaling[i + 2])

    if maximize:
        return likelihood_endo, likelihood_cont, 0, 0, 0, 0, 0

    bwd_5ss_scaled = list(fwd_5ss_scaled)
    bwd_ds_scaled = list(fwd_ds_scaled)
    bwd_ss_scaled = list(fwd_ss_scaled)
    bwd_3ss_scaled = list(fwd_3ss_scaled)

    if forward:
        previous = read.inread[last]
        seqiter = reversed(range(first, last))
    else:
        previous = read.inread[last]
        seqiter = range(last+1,first+1)

    i += 1

    bwd_5ss_scaled[i] = 1
    bwd_ds_scaled[i] = 1
    bwd_ss_scaled[i] = 1
    bwd_3ss_scaled[i] = 1

    previous_j = last
    for j in seqiter:
        i -= 1
        if forward:
            d5 = previous - read.inread[j]
            d3 = length - previous
            previous = read.inread[j]
        else:
            d5 = read.inread[j] - previous
            d3 = previous - 1
            previous = read.inread[j]

        p_d5 = (1 - p) ** d5
        p_d3 = (1 - p) ** d3
        d_d5 = (1 - d) ** d5
        s_d5 = (1 - s) ** d5
        subterm01 = 1 - p_d3 * o2
        subterm02 = p_d3 * o2

        if read.sequence[previous_j] == read.reference[previous_j]:

            tmp_3ss = match[3] * bwd_3ss_scaled[i+1]
            tmp_ds = match[3] * (subterm02 * bwd_3ss_scaled[i+1]
                                + (1 - d_d5) * subterm01 * bwd_ss_scaled[i+1]) \
                     + match[1] * bwd_ds_scaled[i+1] * d_d5 * subterm01
            tmp_ss = match[2] * bwd_ss_scaled[i+1] * s_d5 \
                     + match[1] * bwd_ds_scaled[i+1] * (1 - s_d5)
            tmp_5ss = match[0] * bwd_5ss_scaled[i+1] * p_d5 \
                      + match[1] * bwd_ds_scaled[i+1] * (1 - p_d5)

        elif (forward and read.sequence[previous_j] == "T") or (reverse and read.sequence[previous_j] == "A"):

            tmp_3ss = deam[3] * bwd_3ss_scaled[i+1]
            tmp_ds = deam[3] * (subterm02 * bwd_3ss_scaled[i+1]
                                + (1 - d_d5) * subterm01 * bwd_ss_scaled[i+1]) \
                     + deam[1] * bwd_ds_scaled[i+1] * d_d5 * subterm01
            tmp_ss = deam[2] * bwd_ss_scaled[i+1] * s_d5 \
                     + deam[1] * bwd_ds_scaled[i+1] * (1 - s_d5)
            tmp_5ss = deam[0] * bwd_5ss_scaled[i+1] * p_d5 \
                      + deam[1] * bwd_ds_scaled[i+1] * (1 - p_d5)

        else:

            tmp_3ss = error[3] * bwd_3ss_scaled[i+1]
            tmp_ds = error[3] * (subterm02 * bwd_3ss_scaled[i+1]
                                + (1 - d_d5) * subterm01 * bwd_ss_scaled[i+1]) \
                     + error[1] * bwd_ds_scaled[i+1] * d_d5 * subterm01
            tmp_ss = error[2] * bwd_ss_scaled[i+1] * s_d5 \
                     + error[1] * bwd_ds_scaled[i+1] * (1 - s_d5)
            tmp_5ss = error[0] * bwd_5ss_scaled[i+1] * p_d5 \
                      + error[1] * bwd_ds_scaled[i+1] * (1 - p_d5)


        previous_j = j
        bwd_5ss_scaled[i] = tmp_5ss / fwd_scaling[i + 1]
        bwd_ds_scaled[i] = tmp_ds / fwd_scaling[i + 1]
        bwd_ss_scaled[i] = tmp_ss / fwd_scaling[i + 1]
        bwd_3ss_scaled[i] = tmp_3ss / fwd_scaling[i + 1]


    chain_length=len(fwd_5ss_scaled)-1
    posterior_5ss = [f * b for f, b in zip(fwd_5ss_scaled[1:chain_length], bwd_5ss_scaled[1:chain_length])]
    posterior_ds = [f * b for f, b in zip(fwd_ds_scaled[1:chain_length], bwd_ds_scaled[1:chain_length])]
    posterior_ss = [f * b for f, b in zip(fwd_ss_scaled[1:chain_length], bwd_ss_scaled[1:chain_length])]
    posterior_3ss = [f * b for f, b in zip(fwd_3ss_scaled[1:chain_length], bwd_3ss_scaled[1:chain_length])]

    if forward:
        seqiter = range(len(read.sequence))
    else:
        seqiter = list(reversed(range(len(read.sequence))))

    return likelihood_endo, likelihood_cont, posterior_5ss, posterior_ds, posterior_ss, posterior_3ss, seqiter


def compute_loglikelihood(parameters, reads, maximize=True):
    tmat.cache_clear()
    contam = parameters[6]
    fullLogLikelihood = 0

    e = parameters[0]
    rss = parameters[1]
    rds = parameters[5]
    match = evec("match", e, rss, rds, True)
    deam = evec("deam", e, rss, rds, True)
    error = evec("error", e, rss, rds, True)

    for read in reads:
        likelihood_endo, likelihood_cont, posterior_5ss, posterior_ds, posterior_ss, posterior_3ss, info_positions = read_loglikelihood(read, *parameters, match, deam, error, maximize)
#        likelihood_endo, likelihood_cont, fwd_ds_scaled, bwd_ds_scaled = read_loglikelihood(read, *parameters)
        fullLogLikelihood += logaddexp(likelihood_endo + log(1 - contam), likelihood_cont + log(contam))
    return -fullLogLikelihood


def format_loglikelihood(parameters, reads, per_position=False):
    contam = parameters[6]

    e = parameters[0]
    rss = parameters[1]
    rds = parameters[5]
    match = evec("match", e, rss, rds, True)
    deam = evec("deam", e, rss, rds, True)
    error = evec("error", e, rss, rds, True)

    for read in reads:
        likelihood_endo, likelihood_cont, posterior_5ss, posterior_ds, posterior_ss, posterior_3ss, info_positions = read_loglikelihood(read, *parameters, match, deam, error, maximize=False)
#        likelihood_endo, likelihood_cont, fwd_ds_scaled, bwd_ds_scaled = read_loglikelihood(read, *parameters)

        if per_position:
            try:
                for i,pos in enumerate(info_positions):
                    if read.forward:
                        position = read.inread[pos]
                    else:
                        position = read.inread[pos] - (read.length + 1)
                    yield (
                        position,
                        read.length,
                        read.sequence[pos],
                        posterior_5ss[i],
                        posterior_ds[i],
                        posterior_ss[i],
                        posterior_3ss[i],
                        read.name
                    )
            except TypeError:
                yield (
                    "NA",
                    read.length,
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    read.name
                )
        else:
            try:
                yield (
                    read.name,
                    ''.join([read.sequence[i] for i in info_positions]),
                    read.reverse,
                    read.length,
                    posterior_5ss,
                    posterior_ds,
                    posterior_ss,
                    posterior_3ss,
                    [read.inread[i] for i in info_positions],
                    exp(likelihood_endo) * (1 - contam) / (exp(likelihood_endo) * (1 - contam) + exp(likelihood_cont) * contam),
                    likelihood_endo,
                    likelihood_cont,
                    likelihood_endo-likelihood_cont
                )
            except TypeError:
                yield (
                    read.name,
                    read.sequence,
                    read.reverse,
                    read.length,
                    posterior_5ss,
                    posterior_ds,
                    posterior_ss,
                    posterior_3ss,
                    read.inread,
                    exp(likelihood_endo) * (1 - contam) / (exp(likelihood_endo) * (1 - contam) + exp(likelihood_cont) * contam),
                    likelihood_endo,
                    likelihood_cont,
                    likelihood_endo-likelihood_cont
                )

