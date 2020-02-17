import numpy as np
import pyBigWig
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import norm
from math import log10, exp


def process_bed():
    unsorted_list = []
    fr = open("data/atac-seq/hesc_hg38.bed", 'r')
    items = [line.split() for line in fr]
    fr.close()
    for item in items:
        if item[0] == 'chr21':
            unsorted_list.append([str(item[0]), int(item[1]), int(item[2])])
    sorted_list = sorted(unsorted_list, key=lambda x: x[1])
    with open("data/atac-seq/hesc_chr21.txt", 'w') as fw:
        for item in sorted_list:
            s = item[0] + '\t' + str(item[1]) + '\t' + str(item[2]) + '\n'
            fw.write(s)


def get_open_region(fname):
    """
    bw_hESC = pyBigWig.open(bw_hESC_fname)
    bw_hESC_chrom = dict(bw_hESC.chroms())
    bw_hESC_header = dict(bw_hESC.header())
    bw_hESC.close()
    """
    fr = open(fname, 'r')
    items = [line.split() for line in fr]
    fr.close()

    open_region = []
    for item in items:
        open_region.append([int(item[1]), int(item[2])])
    return open_region


def hic_mix(matrix1, matrix2):
    matrix = np.add(matrix1, matrix2)
    return matrix


def show_HiC_plot(hic_matrix, ran=None, vmin=0, vmax=20):
    if ran:
        hic_matrix = hic_matrix[ran[0]:ran[1], ran[0]:ran[1]]
    trans = hic_matrix.transpose()
    for i in range(0, ran[1]-ran[0]):
        trans[i][i] = 0
    hic_matrix += trans
    sns.set()
    sns.heatmap(hic_matrix, vmin=vmin, vmax=vmax, cmap="Blues")
    plt.show()


def get_open_rate(open_region, resolution, region):
    size = (region[1] - region[0]) / resolution + 1

    open_rate = np.zeros(int(size))
    for interval in open_region:
        start = interval[0] - region[0]
        end = interval[1] - region[0]
        start_idx = int(start / resolution)
        end_idx = int(end / resolution)
        if start_idx != end_idx:
            for i in range(start_idx + 1, end_idx):
                open_rate[i] = 1
            open_rate[start_idx] = start_idx - start/resolution + 1
            open_rate[end_idx] = end/resolution - end_idx
        else:
            open_rate[start_idx] = (end - start) / resolution
    return open_rate


def load_hic_matrix(fname, region, resolution=2000):
    region = [int(region[0] / resolution), int(region[1] / resolution)]
    hic_matrix = np.zeros((region[1]-region[0]+1, region[1]-region[0]+1))
    with open(fname, 'r') as f:
        content = f.readlines()
    content = [x.split() for x in content]
    for data in content:
        p1 = int(data[0])/resolution - region[0]
        p2 = int(data[1])/resolution - region[1]
        v = np.float(data[2])
        hic_matrix[int(p1)][int(p2)] = v
    return hic_matrix


def get_share_function(anc11, anc12, anc21, anc22, decay):
    v1 = anc11 * anc12 / decay
    v2 = anc21 * anc22 / decay
    if v1 == 0 and v2 == 0:
        return 0.5, 0.5
    share = max(0.8, v1/(v1+v2))
    share = min(0.2, share)
    return share, 1-share


def get_share(rate1, rate2):
    sigma = 1
    spread = 3
    lense = len(rate1)
    norm_rate1 = np.zeros(lense)
    norm_rate2 = np.zeros(lense)

    for i in range(0, lense):
        value1 = rate1[i] * 100
        value2 = rate2[i] * 100
        distribution = norm(i, sigma)
        if value1 != 0:
            for j in range(i - spread, i + spread):
                if j < 0 or j > lense:
                    continue
                norm_rate1[j] += distribution.pdf(j) * value1
        if value2 != 0:
            for j in range(i - spread, i + spread):
                if j < 0 or j > lense:
                    continue
                norm_rate2[j] += distribution.pdf(j) * value2

    share_matrix1 = np.zeros((lense, lense))
    share_matrix2 = np.zeros((lense, lense))

    for i in range(0, lense):
        for j in range(i, lense):
            share1, share2 = get_share_function(
                norm_rate1[i], norm_rate1[j],
                norm_rate2[i], norm_rate2[j],
                exp(i-j)
            )
            share_matrix1[i][j] = share1
            share_matrix2[i][j] = share2
    return share_matrix1, share_matrix2


def deconvolute(share, mix):
    return np.multiply(share, mix)
