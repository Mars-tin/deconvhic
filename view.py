import numpy as np
from pyHiC.visualization import *


def merge_to_1kb_resolution(mat):
    new_mat = np.zeros((250, 250))
    for i in range(250):
        for j in range(250):
            new_mat[i, j] = np.sum(mat[i * 5:(i + 1) * 5, j * 5:(j + 1) * 5])
    return new_mat


def vis_micro_atac(cell_type, ch, pos):
    path = '/home/marstin/Documents/deconvhic/data'
    mmm = np.load(f'{path}/{cell_type}/{ch}/{cell_type}_{ch}_200bp_{pos}_{pos + 250000}.npy')

    mat = np.load(f'{path}/hesc-exps-mat.npy')
    mmm = (np.exp(mmm) - 1) * mat
    mmm = merge_to_1kb_resolution(mmm)
    mmm = np.log(mmm + 1)

    atac = np.load(f'{path}/{cell_type}/{cell_type}-atac-seq-hg38-1kb.npy')[pos // 1000: (pos + 250000) // 1000]

    visualize_HiC_epigenetics(mmm, [atac], f'{cell_type}_{ch}_{pos}.png', vmax=3, colorbar=True,
                              epi_ratio=0.15, interval_after_heatmap=0.0, epi_yaxis=False,
                              epi_colors=['black'], x_ticks=['33.35 Mb', '33.60 Mb']
                              )


def vis_micro_atac_dec(cell_type1, cell_type2, ch, pos):
    path = '/home/marstin/Documents/deconvhic/data'
    mat = np.load(f'{path}/hesc-exps-mat.npy')

    # Load cell type 1
    mmm1 = np.load(f'{path}/{cell_type1}/{ch}/{cell_type1}_{ch}_200bp_{pos}_{pos + 250000}.npy')
    mmm1 = (np.exp(mmm1) - 1) * mat
    mmm1 = merge_to_1kb_resolution(mmm1)
    mmm1 = np.log(mmm1 + 1)

    # Load cell type 2
    mmm2 = np.load(f'{path}/{cell_type2}/{ch}/{cell_type2}_{ch}_200bp_{pos}_{pos + 250000}.npy')
    mmm2 = (np.exp(mmm2) - 1) * mat
    mmm2 = merge_to_1kb_resolution(mmm2)
    mmm2 = np.log(mmm2 + 1)


    atac1 = np.load(f'{path}/{cell_type1}/{cell_type1}-atac-seq-hg38-1kb.npy')[pos // 1000: (pos + 250000) // 1000]

    visualize_HiC_epigenetics(mmm1, [atac1], f'{cell_type1}_{ch}_{pos}.png', vmax=3, colorbar=True,
                              epi_ratio=0.15, interval_after_heatmap=0.0, epi_yaxis=False,
                              epi_colors=['black'], x_ticks=['33.35 Mb', '33.60 Mb']
                              )


vis_micro_atac('hff', 'chr2', 33350000)
vis_micro_atac('hesc', 'chr2', 33350000)
