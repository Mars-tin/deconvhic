import numpy as np


def E_Tk(string1, string2):
    return np.sum(string1) * np.sum(string2) / len(string1) ** 2


def r_1k(string1, string2):
    return np.sum(string1 * string2) / len(string1) - E_Tk(string1, string2)


def r_2k(string1, string2):
    return np.sqrt(np.var(string1) * np.var(string2))


def tho_k(string1, string2):
    return r_1k(string1, string2) / r_2k(string1, string2)


def hiCRep(contact1, contact2, L=300):
    if L is None:
        L = np.shape(contact1)[0]
    r1k, r2k, thok, N = [], [], [], []
    for i in range(1, L):
        index = [j for j in range(L - i)]
        index2 = [j + i for j in range(L - i)]
        string1, string2 = contact1[index, index2], contact2[index, index2]
        if not np.isnan(string1).any() and r_2k(string1, string2) > 0:
            r1k.append(r_1k(string1, string2))
            r2k.append(r_2k(string1, string2))
            thok.append(r1k[-1] / r2k[-1])
            N.append(L - i)
    wk = [N[i] * r2k[i] for i in range(len(N))]
    wk = wk / (np.sum(wk))

    return np.sum(wk * thok)
