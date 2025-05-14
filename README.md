# Dmel-adult-central-brain-atlas



Included here are the analysis scripts used in the companion manuscripts:

<br>

**A High-Resolution Atlas of the Brain Predicts Lineage and Birth Order Underly Neuronal Identity.**
\
_Aaron M. Allen, Megan C. Neville, Tetsuya Nojima, Faredin Alejevski, Devika Agarwal, David Sims, and Stephen F. Goodwin_

and

**A Role for Exaptation in Sculpting Sexually Dimorphic Brains from Shared Neural Lineages.**
\
_Aaron M. Allen, Megan C. Neville, Tetsuya Nojima, Faredin Alejevski, and Stephen F. Goodwin_

<br>

These papers investigate transcriptional hetergeneity and sex differences therein in adult _Drosophila melanogaster_ central brain neurons.

<br>
<br>

## Overview


- `alevin_pipeline`: pipeline for running alevin to align fastq files to transcriptome.
    - written by Devika Agarwal with modification by Aaron M. Allen
- `general_analyses`: other analyses scripts and code to generate figure elements.
    - written by Aaron M. Allen  
- `scenic_pipeline`: pipeline for gene regulatory network analysis was performed using SCENIC.
    - written by Devika Agarwal and Lucy Garner, published in Garner et al., 2023
    - modification for this manuscript by Aaron M. Allen
- `seurat_workflow`: workflow for processing post alignment, predominantly in Seurat.
    - written by Aaron M. Allen
- `yamls`: yaml files for conda/mamba environments.
    - `alevin.yaml` - for the `alevin_pipeline`.
    - `pyscenic_py.yaml` - for the `scenic_pipeline`.
    - `seurat.yaml` - for the `seurat_workflow` and `general_analyses`.
    - `solo.yaml` - for running Solo doublet detection script in the `seurat_workflow`.






