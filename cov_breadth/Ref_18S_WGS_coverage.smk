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


# ------------------- TO CHANGE !!! --------------------

REF = "<PATH_TO_REFERENCE_SEQUENCE>/FRE_4583.fa"
DATA = "<PATH_TO_SAMPLE_FOLDER>"
# dictionary expect in this case for 2 .fq.gz to be present in individual samples folders located in DATA
# also focus on samples starting by either H,L or P. 
SAMPLES = {basename(file):glob.glob("%s/*.fq.gz"%file) for file in glob.glob("%s/[H-L-P]*"%DATA)}
ROOT = "<PATH_TO_OUTPUT_FOLDER>"
SLURM_PARTITIONS = [["ei-largemem",550000,4000000,0,64],["ei-medium,ei-long",0,256000,0,32]]


# ------------------- utils ---------------------------

def get_resource(mode, **kwargs):
    return partial(get_resource_real, SLURM_PARTITIONS=SLURM_PARTITIONS, mode=mode, **kwargs)

def get_resource_real(wildcards, input, threads, attempt, SLURM_PARTITIONS="", mode="", mult=2, min_size=10000):
    # return either partition, mem or threads, theses need to be done together
    # because if mem/cpu/threads is over partition definition, then mem/threads needs to change
    # since there can be setting for minimum mem/threads in some partition definition
    # also max threads needs 
    def return_result(mem,partition,threads,mode):
        if mode=="mem":
            return int(mem)
        if mode=="partition":
            return partition
        if mode=="threads":
            return int(threads)

    # change with each attempt
    # Where input.size//1000000 is used convert the cumulative size of input files in bytes to mb, and the tailing 2 could be any arbitrary number based on the specifics of your shell/script requirements.
    mem = max((input.size//1000000) * attempt * mult, attempt*min_size* mult) # this is mb

    # handle case where we are not on a cluster, no partition is defined
    if SLURM_PARTITIONS[0][0]=="":
        partition = ""
        return return_result(mem,partition,threads,mode)

    # choose the best partition, priority to memory
    # HIGH_MEM_PARTITIONS = [["name","min_memory","max_memory","min_threads","max_threads"]]

    #### check memory, get all partition giving mem/1000 (Gb)
    mem_selection = [part for part in SLURM_PARTITIONS if part[2]>=mem]

    # if empty, get the partition with most mem
    if not mem_selection:
        mem_selection = [max(SLURM_PARTITIONS,key=lambda x:x[2])]

    #### check threads, get all partitions giving at least threads
    threads_selection = [part for part in mem_selection if part[4]>=threads]

    # if empty, get partition with most threads
    if not threads_selection:
        threads_selection = [max(mem_selection,key=lambda x:x[4])]

    #### if more than 1 partition remain, 
    # choose the one with, in this order 
    # least min_mem, min_threads, max_mem, max_threads
    final_selection = min(threads_selection,key=lambda x:[x[1],x[3],x[2],x[4]])

    #### decide on mem/threads to use
    partition, min_memory, max_memory, min_threads ,max_threads = final_selection
    mem_final = max(mem, min_memory)
    threads_final = max(threads, min_threads)

    # output
    return return_result(mem_final,partition,threads_final,mode)


def union(intervals):
    sorted_intervals = sorted([tuple(sorted(i)) for i in intervals], key=itemgetter(0))
    if len(sorted_intervals)<1:  # no intervals to merge
        return
    if len(sorted_intervals)==1:  # no intervals to merge
        yield sorted_intervals[0][0],sorted_intervals[0][1]
        return
    # low and high represent the bounds of the current run of merges
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:  # new interval overlaps current run
            high = max(high, iv[1])  # merge with the current run
        else:  # current run is over
            yield low, high  # yield accumulated interval
            low, high = iv  # start new run
    yield low, high  # end the final run


def filter_bam(bam_file, pid_min, breadth_min):
    contig_to_map = defaultdict(list)
    contig_to_map_pid = defaultdict()
    contig_to_map_bid = defaultdict()
    samfile = pysam.AlignmentFile(bam_file, "rb")
    database = samfile.references
    # check mapping id is at 99%
    # check that length of alignment it at least 10% of contig
    fail = defaultdict(list)
    for refmap in samfile:

        fungus = database[refmap.reference_id]
        fungus_len = refmap.reference_length

        contig = refmap.qname
        contig_len = refmap.qlen

        alen = refmap.alen
        bcv = alen/fungus_len
        matchs = sum(el[1] for el in refmap.cigartuples if el[0]==0)
        pid = matchs/float(alen)
        breath_cov = alen/float(contig_len)
        if (pid>=pid_min)&(breath_cov>=breadth_min):
            contig_to_map[contig].append([fungus,fungus_len,contig_len,alen,matchs,pid,bcv,refmap.qstart,refmap.qend,refmap.reference_start,refmap.reference_end])
        else:
            fail[contig].append([fungus,fungus_len,contig_len,alen,matchs,pid,bcv,refmap.qstart,refmap.qend,refmap.reference_start,refmap.reference_end])
    return contig_to_map,fail

def compute_cov_breadth(bam_file,pid_min,ref_len):

    # check mapping id is at 99%
    # check that length of alignment it at least 10% of contig
    contig_to_map,fail = filter_bam(bam_file,pid_min,pid_min)

    # assign contig to dmag, ignore multimap issue, it will add some noise but in the end if something is there, we will see it from breadth of cov
    ref_to_contigs = defaultdict(list)
    for contig,hits in contig_to_map.items():
        # infos needed, 
            # which contigs are mapped over multiple hits on the same ref
            # what position does this correspond to on the ref
        refs_to_hits = defaultdict(list)
        for hit in hits:
            refs_to_hits[hit[0]].append(hit)
        for ref,hits in refs_to_hits.items():
            alignemnt = [[hit[0]]+hit[-4:] for hit in hits]
            ref_to_contigs[ref].append([contig,alignemnt]) 

    ref_intervals = defaultdict(list)
    for ref,contigs in ref_to_contigs.items():
        for contig in contigs:
            for alignemnt in contig[1]:
                ref_intervals[alignemnt[0]].append(alignemnt[-2:])


    # compute breadth of cov
    breadth = 0
    cov = 0
    read_len = np.mean([el[2] for val in contig_to_map.values() for el in val])
    for ref,intervals in ref_intervals.items():
        breadth += sum([end-start for start,end in union(intervals)])/float(ref_len)
        cov += sum([end-start for start,end in intervals])/float(ref_len)

    # here in the quite special case of reads mapped to that particular seq, cov is just total length of intervals divided by cov 

    return breadth,cov, read_len

rule cov_breadth_info:
    input: bam = expand("{{path}}/{sample}_mapped_sorted.bam",sample=SAMPLES),
           bai = expand("{{path}}/{sample}_mapped_sorted.bam.bai",sample=SAMPLES)
    output: "{path}/sample_cov_breadth.tsv"
    params: path = "{path}"
    run:
        path = params["path"]
        results = []
        for sample in SAMPLES:
            bam_file = "%s/%s_mapped_sorted.bam"%(path,sample) 
            breadth,cov, read_len = compute_cov_breadth(bam_file,0.98,4580)
            results.append([sample,breadth,cov, read_len])
        with open(output[0],"w") as handle:
            handle.write("sample\tbreadth\tcov\tmean_read_len\n")
            handle.writelines("%s\n"%"\t".join(map(str,line)) for line in results)



def load_matrix(file,sample_order=None,strain_order=None):
    with open(file) as handle :
        header = next(handle).rstrip().split("\t")[1:]
        strains = []
        matrix = []
        for line in handle : 
            splitlines = line.rstrip().split("\t")
            strains.append(splitlines[0])
            matrix.append(list(map(float,splitlines[1:])))
    matrix = np.array(matrix)
    if sample_order :
        reorder_samples = [header.index(sample) for sample in sample_order]
        matrix = matrix[:,reorder_samples]
        header = sample_order
    if strain_order:
        reorder_strain = [strains.index(strain) for strain in strain_order]
        matrix = matrix[reorder_strain,:]
        strains = strain_order
    return matrix,header,strains


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


# -------------- actual snakemake start ------------------


rule results:
    input: "%s/18S_bam_plot.pdf"%ROOT


rule bwa_index:
    output:  "%s.pac"%REF
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb=get_resource("mem"),
    singularity: "docker://quay.io/annacprice/bwasamtools:1.10"
    shell:   "bwa index {REF}"

# ---- map reads to the assembly contigs--------------------------------------------------
rule bwa_mem_to_bam:
    input:   index = "%s.pac"%REF,
    params: R1 = lambda w:SAMPLES[w.sample][0],
            R2 = lambda w:SAMPLES[w.sample][1]
    output:  "{path}/{sample}_mapped_sorted.bam"
    resources:
        slurm_partition = get_resource("partition",mult=10),
        mem_mb = get_resource("mem",mult=10)
    threads: 32
    singularity: "docker://quay.io/annacprice/bwasamtools:1.10"
    shell:   "bwa mem -t {threads} {REF} {params.R1} {params.R2} | samtools view  -b -F 4 -@{threads} - | samtools sort -@{threads} - > {output}"

rule index:
    input:   "{path}.bam"
    output: "{path}.bam.bai"
    resources:
        slurm_partition = get_resource("partition"),
        mem_mb=get_resource("mem"),
    singularity: "docker://quay.io/annacprice/bwasamtools:1.10"
    shell: "samtools index {input}"


rule cov_bam:
    input: expand("{{path}}/{sample}_mapped_sorted.bam",sample=SAMPLES)
    output: "{path}/18S_bam_plot.pdf"
    run:
        ref_name = "FRE_4583_nt"
        ref_len = 4580
        with PdfPages(output[0]) as pdf:
            for sample in SAMPLES:
                profile = np.zeros(ref_len)
                bam_file = "%s/%s_mapped_sorted.bam"%(ROOT,sample)
                samfile = pysam.AlignmentFile(bam_file, "rb")
                for index,col in enumerate(samfile.pileup(ref_name,0,ref_len,stepper="nofilter")):
                    profile[col.pos-1] = col.n
                plot(profile,sample,ref_len)
                pdf.savefig()
                plt.close()







