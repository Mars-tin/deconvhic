from utils import load_hic_matrix, hic_mix, show_HiC_plot, \
    get_open_region, get_open_rate, get_share, deconvolute, plot_atac_signal
from hicrep import hiCRep


def main():
    resolution = 2000
    region = [1e7, 2e7]
    ran = [2500, 3500]
    # ran = [2900, 3200]
    # ran = [0, 5000]

    # Filenames
    bw_HFF_fname = "data/atac-seq/bigwig/atac-seq-hdd-2m.bw"
    bw_hESC_fname = "data/atac-seq/bigwig/atac-seq-hesc-2m.bw"
    hic_hESC_fname = "data/hic/H1_hESC_chr21_2kb_10m.txt"
    hic_HFF_fname = "data/hic/HFF_chr21_2kb_10m.txt"
    bed_hESC_fname = "data/atac-seq/processed/hesc_chr21_pro.txt"
    bed_HFF_fname = "data/atac-seq/processed/hff_chr21_pro.txt"

    # Plot ATAC-seq signals
    plot_atac_signal(bw_hESC_fname)
    plot_atac_signal(bw_HFF_fname)

    # Process ATAC-seq
    hESC_open_region = get_open_region(bed_hESC_fname)
    HFF_open_region = get_open_region(bed_HFF_fname)

    # Load HiC data
    hic_hESC_matrix = load_hic_matrix(hic_hESC_fname, region, resolution)
    hic_HFF_matrix = load_hic_matrix(hic_HFF_fname, region, resolution)
    print("Data Loading Succeed!")

    # Simulation
    acc = hiCRep(hic_hESC_matrix, hic_HFF_matrix)
    print("The HiCRep score between 2 original HiC datasets is ", acc)

    show_HiC_plot(hic_hESC_matrix, ran)
    show_HiC_plot(hic_HFF_matrix, ran)
    print("Heatmaps generated.")

    # test = hic_hESC_matrix / (hic_hESC_matrix + hic_HFF_matrix)

    hic_mix_matrix = hic_mix(hic_hESC_matrix, hic_HFF_matrix)
    show_HiC_plot(hic_mix_matrix, ran)
    print("Mixed HiC heatmap generated (rate = 1:1).")

    hESC_open_rate = get_open_rate(hESC_open_region, resolution, region)
    HFF_open_rate = get_open_rate(HFF_open_region, resolution, region)
    hESC_share, HFF_share = get_share(hESC_open_rate, HFF_open_rate)

    decov_hESC_matrix = deconvolute(hESC_share, hic_mix_matrix)
    decov_HFF_matrix = deconvolute(HFF_share, hic_mix_matrix)
    print("Deconvolution Succeed!")

    acc = hiCRep(decov_hESC_matrix, hic_hESC_matrix)
    print("The HiCRep score of hESC HiC deconvolution is ", acc)
    acc = hiCRep(decov_hESC_matrix, hic_HFF_matrix)
    print("The HiCRep score of HFF HiC deconvolution is ", acc)
    acc = hiCRep(decov_hESC_matrix, decov_HFF_matrix)
    print("The HiCRep score between 2 deconvoluted datasets is ", acc)

    show_HiC_plot(decov_hESC_matrix, ran)
    show_HiC_plot(decov_HFF_matrix, ran)
    print("Heatmaps generated.")


if __name__ == "__main__":
    main()
