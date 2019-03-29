---
title: Docs
template: blank
---

[Source](http://atlas.gs.washington.edu/mouse-atac/docs/ "Permalink to Mouse ATAC Atlas")

# Mouse ATAC Atlas

## Data Availability and Github

In several portions of this tutorial we may call out files on our [downloads page][1]. Any R/python scripts referenced in this tutorial are available on our [Github page][2].

## Generating Peak by Cell Matrices

Many of the initial steps of processing raw sci-ATAC-seq libraries used for this study are similar to our past efforts. Code and documentation on processing sequencing data can be found on the [Fly ATAC Github][3]. Documentation on various downstream steps can be found in our [Fly ATAC documentation][4].

The above documentation should allow you to obtain a binarized matrix of peaks by cells for a given dataset. Alternatively, you can download the datasets we provide [on our downloads page][1] as a starting point for trying out some analyses.

Below we document some of the most important and generalizable methods that are unique to this study.

## Loading Matrix Files

Through our [downloads page][1] we provide two major formats for all matrix data, matrix-market and RDS.

### Matrix Market Format

Matrix market format is the format used by the `cellranger` pipeline from 10X genomics, which may be familiar to many of you. This format is simply a text file that allows reconstruction of a sparse matrix, along with the peak (or gene) and cell names that specify the row and column names of the matrix, respectively.

This format is readable in `R`: 
    
    
    library(Matrix)
    my_matrix = readMM("my_matrix.mtx.gz")
    row.names(my_matrix) = read.delim('my_matrix.peaks.txt', header=FALSE)$V1
    colnames(my_matrix) = read.delim('my_matrix.cells.txt', header=FALSE)$V1

as well as `python` (based on code from [10X Genomics Documentation][5]). Note the code below requires `scipy` to be installed in your python environment:
    
    
    import csv
    import scipy.io
    
    mat = scipy.io.mmread("my_matrix.mtx.gz")
    peaks = [row[0] for row in csv.reader(open("my_matrix.peaks.txt"), delimiter="t")]
    cells = [row[0] for row in csv.reader(open("my_matrix.cells.txt"), delimiter="t")]

In principle this text format is readable in many other languages as well, although all of our tutorials are written in `R`.

### RDS Format

We also provide these matrices in RDS format, which is easily readable in `R`, but has less cross-compatability with other languages. Row and column names are already assigned within the matrix.
    
    
    library(Matrix)
    my_matrix = readRDS("my_matrix.rds")

## Nucleosome Banding Scores

One of the new quality control steps that we have used in this manuscript is scoring of individual cells for evidence of a banded insert size distribution that would be indicative of a high-quality ATAC-seq library. We do this by first extracting an insert size distribution for each cell and then calculating a score using an FFT of this distribution.

### Performance

The current implementation of the scripts below is a proof of concept and could use further performance optimization to reduce runtimes for datasets with many cells. 

As such, we recommend using a whitelist of cell barcodes to limit computations to non-background cells (cells that meet your reads per cell threshold and other quality control steps). This file, a text file with one cell barcode per line, can be provided with the `\--barcodes` argument in both scripts below.

One bottleneck may be the FFT step, which should be readily optimizable using alternative implementations in R or python.

### Extracting Insert-Size Distributions

The first step is to extract a histogram of insert sizes across all reads associated with each cell. The input to this script is a BAM file of the aligned reads and the output is a per-cell histogram of insert sizes in TSV format.

It assumes that the read names come in the format `cellid:other_text`, which allows the cell ID to be extracted. The script is simple and can be modified as needed if your BAM does not conform to this format.

### Requirements

Requires `pysam` to be installed in your python environment.
    
    
    python get_insert_size_distribution_per_cell.py mybam.bam insert_sizes.txt --barcodes mybarcodes.txt

### Generating Banding Scores

Using `insert_sizes.txt` as input, you can then generate banding scores. 

### Requirements

Script requires `dplyr` and `argparse` to be installed in your R environment.
    
    
    Rscript calculate_nucleosome_banding_scores.R insert_sizes.txt banding_scores.txt --barcodes mybarcodes.txt

The output file in this case is `banding_scores.txt`, a TSV separated file with `cell` and `nucleosome_score` as columns. Lower scores indicate less apparent banding. We recommend plotting a distribution of all scores and setting a cutoff. In our case there was an long left tail that could be used to establish a pass-fail threshold for inclusion in downstream analysis. 

## Monocle Support

### Upcoming Integration

When we started this project no existing tools (to our knowledge) supported both sc-RNA-seq data and sc-ATAC-seq data. sc-ATAC-seq data requires different preprocessing strategies before input into downstream steps such as dimensionality reduction and differential accessibility testing. However, as many functions, such as clustering and dimensionality reduction algorithms themselves, should work for both datatypes, we didn't want to make an toolkit exclusive to sc-ATAC-seq data.

We are aiming to add full support full support for sc-ATAC-seq data to our sc-RNA-seq toolkit, Monocle, by the end of the summer. We also expect that as sc-ATAC-seq becomes more widely used, other popular tools will come to support it as well.

Therefore, while we provide tutorial code below as examples of useful code that mirrors what we have used in our manuscript, we hope that some of it, as well as much additional functionality, will be enabled by existing tools soon.

## Dimensionality Reduction

Our dimensionality reduction approach (commonly referred to in the literature as LSI) has been [documented previously][6] and also used by other groups. However, we provide some additional discussion of it here in the context of this dataset and provide some sample code.

The input to all dimensionality reduction in our manuscript was the binarized count matrix that includes only qc filtered cells available in our downloads [here][7] (`atac_matrix.binary.qc_filtered.rds` for RDS formatted version).

The subset of peaks that we used as input to TFIDF are available in our downloads [here][7] (`atac_matrix.tfidf.qc_filtered.peaks.txt`). The output of TFIDF that we obtained using these sites and the cells in the binary matrix above is also available in our downloads [here][7] (`atac_matrix.tfidf.qc_filtered.rds` for RDS formatted version).

Note that this TFIDF matrix, or subsets of it corresponding to particular clusters, tissues, etc. was used as input to all dimensionality reduction steps (different from the provided example code).

For datasets of your own, we recommend choosing a set of peaks that is present in at least 1-3% of cells or so depending on the exact subset to avoid very lowly sampled peaks. We find that this in addition to avoiding inclusion of cells with very low reads per cell (based on the overall distribution) often improves results.

As a starting point, we also provide an example function on our Github page that given a matrix will do TFIDF, PCA, and t-SNE for you and return the resulting PCA and TSNE coordinates. Note that this function takes the binarized matrix and a `site_frequency_threshold` argument (default 0.03 or site observed in at least 3% of cells). The code contains several hard-coded parameters as inputs to PCA and t-SNE in the example, but you may modify to suit your needs:
    
    
    source('dim_reduction/dim_reduction.R')
    binarized_matrix = readRDS('atac_matrix.binary.qc_filtered.rds')
    results.dim_reduction = atac_dim_reduction(binarized_matrix, site_frequency_threshold=0.02)

This function outputs a list with two items `pca_coords` and `tsne_coords`, which contain the PCA and t-SNE coordinates as dataframes where the cell IDs are included as the rownames.

## Clustering

In our manuscript, we performed clustering in t-SNE space using an older version of Seurat. We expect that many users might instead want to cluster in PCA space (although we expect the results to be broadly similar for this dataset) and use the most recent versions of Seurat, so provide an adapted approach here.

Below we show how to take PCA or t-SNE coordinates and feed them into Seurat to perform Louvain clustering. Again, we imagine that this will become more streamlined once existing tools, such as Monocle, support ATAC-seq data.
    
    
    source('clustering/clustering.R')
    # We provide a utility function to take the results from the dimensionality reduction
    # performed above and put them in a Seurat object, although this code is simple
    # and you may modify it to suit your needs
    seurat_obj = make_seurat(binary_matrix, tsne_coords=results.dim_reduction$tsne_coords, pca_coords=results.dim_reduction$pca_coords)
    
    # You may then use Seurat as normal with desired parameters for clustering
    seurat_obj = FindClusters(seurat_obj, reduction.type='pca', save.SNN=TRUE)
    seurat_obj.tsne_clustering = FindClusters(seurat_obj, reduction.type='tsne', save.SNN=TRUE)
    

## Cicero and Calculating Activity Scores

[Cicero][8], a tool used in various parts of our analyses is now available ([also read the paper here][9]). In addition to outputting maps of coaccessibility, it can also calculate gene activity scores as presented in our manuscript ([see documentation][10]). Please see the paper and documentation for detailed information on Cicero methods and how to start using the tool.

For detailed documentation of Cicero and activity score calculations, we point you to the [Cicero website][8].

## KNN-based scRNA-seq Integration

In our manuscript we make comparisons between gene level scores derived from our data (activity scores) and gene expression measurements from several sc-RNA-seq studies. As mentioned above, you may calculate activity scores using Cicero, or alternatively, this aproach can take other gene by cell matrices as input. For example, we have used aggregate read counts in 25-50kb regions around TSSs, aggregate binarized counts in peaks around TSSs, or even binarized counts in peaks overlapping promoters as inputs with reasonable success as well.

We are in the process of wrapping up this functionality into an R package of its own. We may also have an example function available between now and when this package is made available.

## Differentially Accessible and Specific Peaks 

In our study we calculate differential accessibility of peaks (see documentation on how this is done on our [Fly ATAC page][11]) against a random set of 2000 cells from the remainder of the dataset to make computation more feasible.

### Calculating locus-by-cluster specificity scores

### Performance

Warning: The current script is not particularly optimized for speed and so running it on the full set of sites and clusters can take a few hours to finish.

As mentioned in the manuscript, this analysis relies heavily on the framework first described in [Cabili et al.][12] and implemented in [Monocle][13] \- specificity scores were computed by directly adapting the Monocle code for the `calculateMarkerSpecificity()` function to our ATAC dataset. We start with a matrix of locus-by-cluster proportions (i.e. how many cells in a cluster have each site accessible). We then normalize this matrix (to account for differences in read depth between different clusters of cells) and then calculate Jensen-Shannon specificity scores on this normalized matrix. Finally, we square the scores and multiple them again by the normalized proportion matrix (to emphasize the sites that are both restricted to a cluster _AND_ fairly common in that cluster. Finally, we try to determine an appropriate cutoff for calling a locus "specific" to a cluster. To do so, we use the list of loci that are not differentially accessible in that cluster as a null distribution for determining a specificity score threshold at which no more than 10% of the called "specific" loci belong to the null class.

The script provided here requires 3 files as inputs and requires that you specify 4 outputs. The required inputs are: (1) a list (saved in RDS format) of the results from the differential accessibility tests done for each cluster. Each item in the list (1 for each cluster) is a matrix where the rownames are the locus IDs, the first column (`beta`) contains the "beta" values for the differential accessibility tests, the second column (`pval`) contains the p values, and the third column (`qval`) contains the q values. (2) a matrix of proportions representing the fraction of cells each locus is accesible in for each cluster. Each column is a cluster and each row is a locus (column names and row names should be specified). Again this should be saved as an RDS file. And finally (3), a simple text matrix of the median read depths for cells in each cluster. For this matrix, the first column is the cluster name and the second column is the read depth measure. The 4 output file names you will be asked to specify are for a text file of sites that are specific to each cluster (`-ST`), an rds file for the same information (`-SR`), a text file for all specificity scores for all clusters by all loci (`-AT`), and an RDS file for all specificity scores for all clusters by all loci (`-AR`). 

Below, we present some code to generate the specificity scores we used for the paper:

### Requirements

Requires `Matrix`, `limma`, `reshape2`, and `argparse` to be installed in your R environment.
    
    
    wget http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Specificity_Scores/specificity_data.tgz
            tar -xzf ./specificity_data.tar.gz
            Rscript specificity_calculator.R -pvals ./specificity_data/peak_pvals_for_vignette.rds -props ./specificity_data/peak_proportions_for_vignette.rds -DP ./specificity_data/peak_depths_for_vignette.txt -ST ./specific_peaks.txt -SR ./specific_peaks.rds -AT ./all_peak_specificity_scores.txt -AR ./all_peak_specificity_scores.rds

If you wanted to run this on gene activity scores instead of peak accesibility, the command would be the same, but you would just specify the gene activity files included in tarball instead of the peak accessibility files.
    
    
    Rscript specificity_calculator.R -pvals ./specificity_data/gene_pvals_for_vignette.rds -props ./specificity_data/gene_proportions_for_vignette.rds -DP ./specificity_data/gene_depths_for_vignette.txt -ST ./specific_genes.txt -SR ./specific_genes.rds -AT ./all_gene_specificity_scores.txt -AR ./all_gene_specificity_scores.rds

## LD Score Regression

In our manuscript we use a tool called LD score regression (LDSC) to compute enrichments for trait heritability within differentially accessible peaks across different cell types.

LDSC already has extensive documentation provided on [Github][14] and their [wiki][15]. We take a set of peaks for each cell type and then run LDSC on set individually, annotating which SNPs used by LDSC fall in the cell type's peaks. One model is generated per cell type in this fashion and then used to calculate heritability enrichments for many phenotypes. In our case, we only annotate SNPs as falling in DA peaks if they lift over from the human genome to mouse, but this would be unncessary for human data.

As LDSC has been widely used outside of our study, we expect others will likely develop their own tools for generating LDSC models. However, here we provide modified versions of scripts we used to run LDSC internally. Note that we provide these mostly as a reference to see how we did things and toserve as an example for your own analyses -- we do not plan to support LDSC directly, which already has extensive documentation and user support.

### Warning

We found that generating many LDSC models required us to parallelize. To simplify the script, we show it using serial execution, but this may not be a tractable option in your case and some modifications may be required. LDSC authors may also be able to provide advice on generating many models without requiring parallelization.

### Requirements

There are a few requirements to use LDSC and these scripts specifically. LDSC must be installed (see [Github][14]) and you must have downloaded the requisite files required to use LDSC. Python 2 is required for LDSC. bedtools must be installed. You must be able to point to the path of the liftOver utility and a relevant chain file as described below.
    
    
    python ldscore_regression/ldscore_peaks.py 
        --output_directory outputs 
        --sample_sheet samplesheet.txt 
        --sumstats sumstats.txt 
        --master_peaks all_peaks.bed 
        --ldsc_path ldsc.py 
        --lift_over_path liftOver 
        --liftover_chain hg19ToMm9.over.chain.gz 
        --baseline_prefix ldsc/1000G_Phase3_baseline_v1.1/baseline. 
        --annot_template_prefix ldsc/1000G_Phase3_cell_type_groups/cell_type_group.1 
        --bfile_prefix ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC 
        --hapmap_snps_prefix ldsc/hapmap3_snps/hm 
        --ld_score_prefix ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC 
        --frqfile_prefix ldsc/1000G_Phase3_frq/1000G.EUR.QC

Note that this script relies on being in the same directory as the other scripts provided in the `ldscore_regression` folder on our Github. Documentation for the script arguments is provided below for reference:

| Argument                    | Description                                                                                                                                                                                                                                    |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `\--output_directory`:      | Folder to output results to.                                                                                                                                                                                                                   |
| `\--sample_sheet`           | A TSV file with columns `sample_id` and `sites`, which are a name for a sample (say a cluster ID, for example) and a BED file with sites to use for enrichment testing (cell type differentially accessible peaks, for example).               |
| `\--sumstats`               | A TSV file with columns `phenotype` and `sumstats` that are a name for the phenotype and the path to a `.sumstats.gz` file respectively.                                                                                                       |
| `\--master_peaks`           | A BED file with all peaks in the entire dataset vs. just the specific ones. Not really necessary at this point, but the script does require this argument.                                                                                     |
| `\--liftover_path`          | path to liftOver utility as downloaded from UCSC.                                                                                                                                                                                              |
| `\--liftover_chain`         | A liftover chain file such as `hg19ToMm9.over.chain.gz` or other chains available via [UCSC][16].                                                                                                                                              |
| `\--ldsc_path`              | path to ldsc.py from LD score regression installation.                                                                                                                                                                                         |
| `\--score_sumstats_options` | Optional string containing options for use with LDSC in scoring sumstats files. Must start with a space and be quoted. Example: " --chisq-max 99999999999"                                                                                     |
| `\--baseline_prefix`        | Prefix of baseline model for use with LDSC, such as `ldsc/1000G_Phase3_baseline_v1.1/baseline.`, for example. We had to rerun the basline model ourselves rather than use the one provided for download, but it may be different in your case. |
| `\--annot_template_prefix`  | Prefix for `.annot.gz` files provided with LD score regression. `ldsc/1000G_Phase3_cell_type_groups/cell_type_group.1` in our case.                                                                                                            |
| `\--bfile_prefix`           | Prefix for files used in `\--bfile` argument to LD score regression. `ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC` in our case.                                                                                                                   |
| `\--hapmap_snps_prefix`     | Prefix for files used in `\--hapmap_snps` argument to LD score regression. `ldsc/hapmap3_snps/hm` in our case.                                                                                                                                 |
| `\--ld_score_prefix`        | Prefix for `\--w-ld-chr` argument to LD score regression. `ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC` in our case.                                                                                                                |
| `\--frqfile_prefix`         | Prefix for `\--frqfile-chr` argument used with LD score regression. `ldsc/1000G_Phase3_frq/1000G.EUR.QC` in our case.                                                                                                                          |

As written, the script outputs several directories `make_template_files`, `train_models`, and `score_sumstats` that contain files that are inputs or outputs of LD score regression commands. It also contains a file `results_gathered.txt`, which contains various statistics aggregated from all LD score regression results. See LD score regression documentation for more information on reported statistics.

[1]: http://atlas.gs.washington.edu/data/
[2]: https://github.com/shendurelab/mouse-atac
[3]: https://github.com/shendurelab/fly-atac
[4]: http://atlas.gs.washington.edu/fly-atac/docs/#use-case-1-id-clades-with-lsi
[5]: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
[6]: http://atlas.gs.washington.edu/fly-atac/docs/#use-case-2-re-cluster-cells-with-t-sne
[7]: http://127.0.0.1:4000/mouse-atac/data/#atac-matrices
[8]: https://cole-trapnell-lab.github.io/cicero-release/
[9]: https://www.cell.com/molecular-cell/fulltext/S1097-2765(18)30547-1
[10]: https://cole-trapnell-lab.github.io/cicero-release/docs/#cicero-gene-activity-scores
[11]: http://atlas.gs.washington.edu/fly-atac/docs/#use-case-3-differential-accessibility
[12]: https://www.ncbi.nlm.nih.gov/pubmed/21890647
[13]: http://cole-trapnell-lab.github.io/monocle-release/
[14]: https://github.com/bulik/ldsc
[15]: https://github.com/bulik/ldsc/wiki
[16]: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/
