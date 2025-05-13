from ruffus import *
from ruffus.combinatorics import *
from cgatcore import pipeline as P

import sys
import os

PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
      "../pipeline.yml",
      "pipeline.yml"])

# Helper functions mapping tracks to conditions, etc
# determine the location of the input files (reads).
try:
    PARAMS["input"]
except KeyError:
    DATADIR = "."
else:
    if PARAMS["input"] == 0:
        DATADIR = "."
    elif PARAMS["input"] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS["input"]  # not recommended practise.

#########################################################################
#########################################################################
#########################################################################
# define input files
SEQUENCESUFFIXES = ("*.fastq.1.gz",
                    "*.fastq.gz",
                    "*.fa.gz",
                    )

SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                       for suffix_name in SEQUENCESUFFIXES])

SEQUENCEFILES_REGEX = regex(
    r"(.*\/)*(\S+)_[^2-3].(fastq.1.gz|fastq.gz|fa.gz)")

SEQUENCEFILES_REGEX_ALT = regex(
    r"(.*\/)*(\S+)_.*[1-8].(fastq.1.gz|fastq.gz|fa.gz)")

@active_if(PARAMS["fastq_lanes"]== "single")
@follows(mkdir("Alevin_r_6.30_allcb"))
@transform( SEQUENCEFILES, SEQUENCEFILES_REGEX, r"Alevin_r_6.30_allcb/\2_rna/.alevin_done")
def alevin_rna(infile, outfile):
    '''
    Docstring
    '''
    infile2 = infile.replace("_1.fastq.gz", "_2.fastq.gz")
    threads = PARAMS["alevin_threads"]
    outfolder = outfile.replace("/.alevin_done", "")
    statement = """
    salmon  alevin -l %(alevin_library_type)s -i %(alevin_index)s \
        -1 %(infile)s -2 %(infile2)s \
        -o %(outfolder)s -p %(threads)s  --tgMap %(alevin_tgmap)s \
        %(alevin_options)s  \
        2> %(outfolder)s/job.log
        > %(outfolder)s/job.err &&
    touch %(outfile)s
    """

    job_threads = PARAMS["alevin_threads"]

    P.run(statement, job_memory=PARAMS["alevin_rna_memory"],job_condaenv=PARAMS["alevin_conda_env"])


@active_if(PARAMS["fastq_lanes"]== "multiple")
@follows(mkdir("Alevin_r_6.30_allcb"))
@collate(SEQUENCEFILES, SEQUENCEFILES_REGEX_ALT, r"Alevin_r_6.30_allcb/\2_rna/.alevin_done")
def alevin_rna_alt(infiles, outfile):
    '''
    Docstring
    '''
    infiles2 = [file.replace(".fastq.1.gz", ".fastq.2.gz") for file in infiles]

    fastq1 = " ".join(infiles)
    fastq2 = " ".join(infiles2)

    threads = PARAMS["alevin_threads"]
    outfolder = outfile.replace("/.alevin_done", "")
    statement = """
    salmon alevin -l %(alevin_library_type)s -i %(alevin_index)s \
        -1 %(fastq1)s -2 %(fastq2)s \
        -o %(outfolder)s -p %(threads)s  --tgMap %(alevin_tgmap)s \
        %(alevin_options)s  \
        2> %(outfolder)s/job.log
        > %(outfolder)s/job.err &&
    touch %(outfile)s
    """

    job_threads = PARAMS["alevin_threads"]

    P.run(statement, job_memory=PARAMS["alevin_rna_memory"],job_condaenv=PARAMS["alevin_conda_env"])

@transform([alevin_rna,alevin_rna_alt], regex(r"Alevin_r_6.30_allcb/(.*)_rna/(.*)") , r"Alevin_r_6.30_allcb/\1_rna/alevin/.done")
def droplet_utils(infile, outfile):
    '''
    Docstring
    '''
    infile = infile.replace(".alevin_done", "alevin/quants_mat.gz")

    outfolder = outfile.replace("/.done", "")
    sample_name = outfolder.split("/")[1]
    #sample_name =  os.path.basename(outfolder).replace("_rna/alevin", "")
    filter_fdr = PARAMS["emptydrops_filter_fdr"]
    filter_empty = PARAMS["emptydrops_filter_empty"]
    filter_genes = PARAMS["emptydrops_filter_genes"]
    plot_empty_pval = PARAMS["emptydrops_plot_empty_pval"]
    lower = PARAMS["emptydrops_lower"]
    niters = PARAMS["emptydrops_niters"]
    test_ambient = PARAMS["emptydrops_test_ambient"]
    ignore = PARAMS["emptydrops_ignore"]
    retain = PARAMS["emptydrops_retain"]
    working_dir = PARAMS["emptydrops_figures_dir"]
    working_dir = working_dir.replace("/figures_dir","")
    figures_dir = PARAMS["emptydrops_figures_dir"] + "/alevin_dropletutils"
    barcode_plot_function = PARAMS["emptydrops_barcode_plot_function"]
    R_PATH = PARAMS["R_path"]


    '''
    Run the R script Alevin_droplet_write10x_emptydrops.R to get cellranger format raw and filtered counts
    after calling empty drops function
    '''
    statement = """Rscript %(R_PATH)s/Alevin_droplet_write10x_emptydrops.R
        --input_file=%(working_dir)s/%(infile)s
        --output_dir=%(outfolder)s
        --sample_label=%(sample_name)s
        --filter_fdr=%(filter_fdr)s
        --filter_empty=%(filter_empty)s
        --filter_genes=%(filter_genes)s
        --plot_empty_pval=%(plot_empty_pval)s
        --lower=%(lower)s
        --niters=%(niters)s
        --test_ambient=%(test_ambient)s
        --ignore=%(ignore)s
        --retain=%(retain)s
        --figures_dir=%(figures_dir)s
        --barcode_plot_function=%(R_PATH)s/%(barcode_plot_function)s
    2> %(working_dir)s/%(outfolder)s/job.log
    > %(working_dir)s/%(outfolder)s/job.err &&
    touch %(outfile)s
    """

    P.run(statement, job_memory=PARAMS["emptydrops_memory"],job_cores=PARAMS["emptydrops_cores"],job_condaenv=PARAMS["emptydrops_conda_env"])


@follows(droplet_utils)
def full():
    pass



if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
