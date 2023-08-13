import numpy as np
from multiprocessing.pool import Pool

def qn(phi, n) -> complex:
    phi = (np.array([phi])).astype(complex)
    qn_vector = np.exp(1j*n*phi).sum()
    return qn_vector

def eve_weight(mult, n) -> float:
    weight = 1.
    for i in range(mult - n + 1, mult + 1):
        i_float = float(i)
        weight = i_float*weight
    return weight

def eve_weight_etagap(mult1, mult2) -> float:
    return mult1*mult2

def single_ave2(qn, mult) -> float:
    qn_conj = np.conjugate(qn)
    ave2_single_event = (((qn*qn_conj).real - mult)/(eve_weight(mult, 2)))
    return ave2_single_event

def single_ave2_with_eta_gap(qn_1,qn_2,mult_1,mult_2) -> float:
    qn_2_conj = np.conjugate(qn_2)
    return float((qn_1*qn_2_conj).real)/float((mult_1*mult_2))

def single_ave4(qn,q2n,mult) -> float:
    qn_conj = np.conjugate(qn)
    q2n_conj = np.conjugate(q2n)
    abs_qn_fp = (qn*qn*qn_conj*qn_conj).real
    abs_q2n_sq = (q2n_conj*q2n).real
    single_ave4_numerator = (abs_qn_fp + abs_q2n_sq - 2*(q2n*qn_conj*qn_conj).real - 
                4*(mult - 2)*(qn*qn_conj).real + 2*mult*(mult - 3))
    single_ave4_denominator = eve_weight(mult, 4)
    ave4_single_event = single_ave4_numerator/single_ave4_denominator
    return ave4_single_event

def meanpT(pT_array):
    return np.mean(pT_array)


def runobs(phi_std,phi_sub_1,phi_sub_2,pT_array):
    mult_std = len(phi_std)
    mult_sub_1 = len(phi_sub_1)
    mult_sub_2 = len(phi_sub_2)
    single_ave2_2_eta = single_ave2_with_eta_gap(qn(phi_sub_1,2),qn(phi_sub_2,2),mult_sub_1,mult_sub_2)
    single_ave4_2 = single_ave4(qn(phi_std,2),qn(phi_std,4),mult_std)
    single_ave2_3_eta = single_ave2_with_eta_gap(qn(phi_sub_1,3),qn(phi_sub_2,3),mult_sub_1,mult_sub_2)
    single_ave4_3 = single_ave4(qn(phi_std,3),qn(phi_std,6),mult_std)
    single_ave2_4_eta = single_ave2_with_eta_gap(qn(phi_sub_1,4),qn(phi_sub_2,4),mult_sub_1,mult_sub_2)
    single_ave4_4 = single_ave4(qn(phi_std,4),qn(phi_std,8),mult_std)
    meanpt = meanpT(pT_array)
    return single_ave2_2_eta, single_ave4_2, single_ave2_3_eta, \
           single_ave4_3, single_ave2_4_eta, single_ave4_4, meanpt
