from os.path import basename,dirname
from collections import defaultdict
from operator import itemgetter
from functools import partial
import scipy.stats as ss
import numpy as np
import glob
import pysam
from scipy.stats.stats import pearsonr

try:
    import matplotlib
    matplotlib.use('Agg')
finally:
    import matplotlib.pyplot as plt
    import matplotlib.lines as lines
    import matplotlib.patches as patches
    from matplotlib.ticker import NullFormatter

from matplotlib.backends.backend_pdf import PdfPages


plot_reduction_level = 1
color_scheme = {'plot_orf': '#FF7F7F',
                'plot_inter':"#3776ab",
                'primerP1_fc':'#ff5353',
                'primerP2_fc':'#5fefc7',
                'gene_fc':'#ffe17f',
                'mismatch_normal':'#ffe17f',
                'mismatch_primer':'red',
                }

def set_plot_area(ax, x_max,max_hight=10000):

    # Setting x and y labels
    ax.set_ylabel('Depth')
    ax.set_xlabel('Genome position nt')

    ax.axhline(1, color='k', linewidth=0.5, zorder=102) # line at depth=1

    # Ticks
    #ymax = np.ceil(tbl['DEPTH'].max() / 10000) * 10000

    ax.set_xlim(0, x_max)
    ax.set_yscale('log')
    ymin = 10**-(np.log10(max_hight)/5)
    # For linear scale
    # ymin = -(max_hight)/5

    ax.set_ylim(ymin, max_hight)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())

    y_major_ticks = [el for el in ax.get_yticks() if el >=1 and el <= max_hight]
    y_minor_ticks = [el for el in ax.get_yticks(minor=True) if el >1 and el < max_hight]
    ax.set_yticks(y_major_ticks)
    ax.set_yticks(y_minor_ticks, minor=True)
    ax.set_yticklabels([str(int(el)) for el in y_major_ticks])

def plot(profile,sample,ref_len):
    fig, ax = plt.subplots(1,1, figsize=(8, 3))
    plt.subplots_adjust(right=0.85,
                        hspace=0.5,
                        bottom=0.5/3,
                        top=1-0.5/3)
    ax.set_title("coverage for sample %s"%sample)
    set_plot_area(ax,ref_len,max_hight=10000)
    x=[0]+list(range(len(profile)))[::plot_reduction_level]+[len(profile)]
    y = [0.9]+list(profile)+[0.9]
    ax.fill(x,y,color_scheme['plot_orf'],zorder=52)
    labels = [item.get_text() for item in ax.get_yticklabels()]


def main(bam_file,REF_NAME,REF_LEN,output_pdf):
    with PdfPages(output_pdf) as pdf:
        profile = np.zeros(REF_LEN)
        samfile = pysam.AlignmentFile(bam_file, "rb")
        for index,col in enumerate(samfile.pileup(REF_NAME,0,REF_LEN,stepper="nofilter")):
            profile[col.pos-1] = col.n
        plot(profile,sample,REF_LEN)
        pdf.savefig()
        plt.close()




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam_file", help="bam file resulting for mapping reads to reference sequence")
    parser.add_argument("ref_name", help="name of the reference contig")
    parser.add_argument("ref_len", help="length of reference contig")
    parser.add_argument("output", help="path/name of output pdf file")
    args = parser.parse_args()
    bam_file = args.bam_file
    REF_NAME = args.ref_name
    REF_LEN = args.ref_len
    output = args.output
    main(bam_file,REF_NAME,REF_LEN,output_pdf)
