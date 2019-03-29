
[Source](http://atlas.gs.washington.edu/mouse-atac/data/ "Permalink to Mouse ATAC Atlas")

# Mouse ATAC Atlas

## Release Updates

* `08/14/2018` adds `UCSC Trackhub and BigWigs` section which provides documentation on new bigWig files and a UCSC Trackhub for visualization.
* `08/10/2018` adds another metadata file option, `cell_metadata.tissue_freq_filtered.txt`, which removes cells in each tissue belonging to a `cell_label` that accounts for less than 0.5% of the cells in that tissue. These very low frequency labels are often not cell types expected their respective tissues and could be due to slight imperfections in clustering, for example.
* `08/06/2018`: 
    * Switches any references to Colon to LargeIntestine and Intestine to SmallIntestine. This impacts any file that lists a tissue name, but all files were previously consistent with one another and usable.
    * Updates specificity score section including added tutorial input files and minor updates to existing gene-level specificity scores.
    * Adds a `BAM Files` for BAM downloads.
* `08/02/2018`: Initial release.

## ATAC Matrices

Similar to sc-RNA-seq, sci-ATAC-seq data is typically analyzed in sparse `peak` (row) by `cell (column)` matrices. The first set we provide are binarized counts. The second set has rare peaks filtered out and is then normalized with TFIDF to allow for input to PCA/TSNE, for example. Note that only cells in our final QC filtered set are included.

See [tutorials][1] for examples of how to read these formats into `R` or `python` along with documentation on lots of other downstream analysis.

### Matrix Formats

In general, we provide two formats for all matrices:

| File Type | Description |  
| --------- | ----------- |  
| `.mtx.gz` |  

Same as format generated by 10X Genomics `cellranger` pipeline (matrix market format). `.mtx` files are provided with `.txt` files containing peak and cell IDs that correspond to the rows and columns of the matrix, respectively.

 |  
| `.rds` | 

RDS format that can be read into R directly with `readRDS` function.

 | 

## Activity Score Matrices

We also report "gene activity scores", where a single number is calculated based on a weighted combination of proximal and distal sites for each gene (see manuscript for details; both quantitative and binarized calculations provided below). Unlike the ATAC matrices above, these are in `gene (row)` by `cell (column)` format.

### Important: Size Factor Normalization

Quantitative scores provided here are not normalized by size factors, so you may apply size factor normalization to these values if needed.

## Metadata

For all cells and peaks used in our QC filtered set, we report tables of metadata including information about tissue source, cell type assignment, TSNE coordinates, cluster assignments, etc. for cells, and intersections with genes (TSS only) for peaks.

| Name                                        | Last modified | Size    | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                  |  
| ------------------------------------------- | ------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |  
| [cell_metadata.txt][2]                      | 08/06/18      | 13.3MB  | Metadata for cells in TSV format, including several features such as TSNE coordinates, cluster assignments, and cell type assignments. [format details][3]                                                                                                                                                                                                                                                                                                   |  
| [cell_metadata.tissue_freq_filtered.txt][4] | 08/10/18      | 13.0MB  | Same as `cell_metadata.txt`, but removes cells in each tissue belonging to a `cell_label` that accounts for less than 0.5% of the cells in that tissue. These very low frequency labels are often not cell types expected their respective tissues and could be due to slight imperfections in clustering, for example. Provided matrices would need to be subsetted to match this set of cells if using this metadata.                                      |  
| [peak_promoter_intersections.txt][5]        | 08/02/18      | 6.9MB   | Metadata for peaks-intersected TSS pairs in TSV format [format details][6]                                                                                                                                                                                                                                                                                                                                                                                   |  
| [cell_type_assignments.xlsx][7]             | 08/06/18      | 134.1KB | Excel document with three tabs `expected cell types`, `cell type markers`, and `Cell type assignments` that contain a pairs of tissues and expected cell types, a list of positive markers for each cell type, and the table of cell type assignments with extra details about assignment criteria when applicable, respectively. This is meant to document justification for cell type assignments provided in `cell_metadata.txt` above. [format details][8]  

`expected cell types`

* `tissue`: the tissue in the cell type tissue pair
* `expected_cell_type`: the cell type in the tissue cell type pair
`cell type markers`

* `Cell type`: the cell type from above table
* `subset_cluster`: a list of positive markers for that cell type
`Cell type assignments`

* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space
* `pipeline_cell_type`: broad automatic assignment made by classifier (see manuscript).
* `manual_cell_type`: a broad assigned cell type where applicable
* `cell_label`: assigned cell type (same as in cell metadata file)
* `notes`: extra notes about assignment criteria where applicable

 | 

## Differential Accessibility

We report results from differential accessibility (DA) tests performed between each cluster of cells (final iterative clusters) and a set of 2K sampled cells. See manuscript for details. These are reported using both the binarized ATAC matrix (contains all peaks) and the binarized gene activity score matrix (contains a single entry per gene).

### DA Test Format

The following columns are provided in each DA test file. Note that files contain one entry per combination of `cluster`, `subset_cluster`, and `peak/gene_short_name`.

| Column                 | Description                                                                                                                                                                                               |  
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |  
| `status`               | status of successful test completion reported by monocle (OK for all in this case)                                                                                                                        |  
| `family`               | the distribution family used by monocle (binomialff in this case)                                                                                                                                         |  
| `pval`                 | Uncorrected P-value returned by monocle.                                                                                                                                                                  |  
| `beta`                 | Beta derived from the model returned by monocle. The coefficient of the term noting cluster membership. Negative values indicate less accessibility within the specified cluster than the 2K sampled set. |  
| `qval`                 | Q-value returned by monocle corrected for all tests performed for the specified cluster.                                                                                                                  |  
| `peak/gene_short_name` | For ATAC matrix this is the peak ID of the peak being tested (chr_start_stop; header `peak`). For activity score matrix, this is the common gene name of the gene being tested (header `gene`).           |  
| `cluster`              | Cluster assignment in initial t-SNE.                                                                                                                                                                      |  
| `subset_cluster`       | Cluster assignment in iterative t-SNE space.                                                                                                                                                              |  

## Specificity Scores

We report specificity scores to rank elements by their restricted accessiblity within each of our clusters (see manuscript for details). Only sites that had significant specificity scores at our empirically determined false discovery rate threshold are reported.

These are provided in both Excel format and text format to allow browsing of results.

## Comparisons to scRNA-seq Datasets

In our manuscript we also examine similarity between our sci-ATAC-seq dataset and several sc-RNA-seq datasets. We do this using a cluster-level correlation-based approach and a cell-by-cell KNN-based approach. Both used activity scores as calculated by Cicero as input (see above). Here we provide the cluster-level correlations for each dataset/tissue and the cell-by-cell KNN results for each dataset we have compared to.

| Name                 | Last modified | Size  | Description                                                                                                                                                                                                                                                                                                        |  
| -------------------- | ------------- | ----- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |  
| [knn_results.txt][9] | 08/06/2018    | 8.0MB | KNN results listing assignments of each cell in our study in TSV format. Some cells will be missing and designated as "NA" due to the filtering of very low frequency cell type assignments within each tissue and thresholds on the number of non-zero peaks per cell as described in methods. [format details][10]  

* `cell`: cell barcode (combined and corrected)
* `tissue`: the tissue in the cell type tissue pair
* `cell_label`: assigned cell type in our study
* `assigned_cell_label_microwell`: Label assigned from MCA dataset for applicable cells (Han et al.). See columns below for explanation of "NA" values.
* `assigned_cell_label_tm`: Label assigned from Tabula Muris Consortium dataset for applicable cells. See columns below for explanation of "NA" values.
* `label_status_microwell`: In cases where labels assigned from microwell dataset are "NA" this column provides a reason why. Entries with non-NA labels assigned with have "assigned" in this column. NA valuesin in the `assigned_cell_label_microwell` column can be NA due to tissues not overlapping between the two datasets (not_compared), a cell not being included in a comparison due it not meeting our thresholds for number of non-zero sites (filtered_by_depth), or KNN neighbors not providing a clear majority label (low_support).
* `label_status_tm`: same as label_status_microwell but referring to the assigned Tabula Muris label.

 |  
| [correlation_results.txt][11] |  08/06/2018 |  201.9KB |  Cluster-level Spearman correlations as a metric of similarity between each of our chromatin profiles and cell type assignments made in each of the same sc-RNA-seq studies in TSV format. [format details][12]

* `id`: combined ID for major + iterative cluster assignment
* `tissue`: the tissue in the cell type tissue pair
* `cell_label`: assigned cell type in our study for the specified cluster
* `scrna_seq_dataset`: The sc-RNA-seq study used to compute the correlation. One of "tabula_muris" or "microwell".
* `scrna_seq_label`: Label assigned from sc-RNA-seq dataset.
* `correlation`: Spearman correlation between the activity scores and expression values for genes with variable gene expression as described in manuscript.

 |  
| [microwell_truncated_labels.txt][13] |  08/02/2018 |  6.6KB |  As described in our manuscript, to facilitate comparisons to the MCA dataset, we made minor modifications to the labels provided in this study (see methods for details). This file provides original labels used in the MCA paper and then the set of labels that we used in TSV format. [format details][14]

* `mca_label`: the original label for this cell provided in the MCA data release
* `mca_label.modified`: the modified label that we used for comparisons in the manuscript

 | 

## Cicero Maps

We have also run Cicero ([Pliner _et al._][15]), which connects regulatory elements to their target genes using coaccessibility as a measure of connectedness, as measured by sci-ATAC-seq. We have generated Cicero maps for each cluster in the dataset. Maps and peak sets are combined into single files with columns to indicate the cluster and subset_cluster entries correspond to.

| Name                          | Last modified | Size  | Description                                                                                                                                                                   |  
| ----------------------------- | ------------- | ----- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |  
| [master_cicero_conns.txt][16] | 01/11/2018    | 3.1GB | Cicero was run on the peak sets provided below (accessible in at least 1% of cells for a given cluster). This is the resulting map provided in TSV format. [format details][17]  

* `Peak1`: Peak ID (see below) for first peak in connection.
* `Peak2`: Peak ID (see below) for second peak in connection.
* `coaccess`: coaccessibility score for connection. Higher values indicate stronger signal. See [Pliner _et al._][18] for details.
* `peak1.isproximal`: Peak1 is proximal if it has overlap with 5kb upstream and 1kb downstream of any TSS (set to "Yes" if proximal).
* `peak1.tss.gene_name`: Common gene name of the proximal TSS overlapping Peak1
* `peak1.tss.gene_id`: Ensemble gene id of the proximal TSS overlapping Peak1
* `peak2.isproximal`: same as for peak1
* `peak2.tss.gene_name`: same as for peak1
* `peak2.tss.gene_id`: same as for peak1
* `conn_type`: one of proximal_proximal, distal_proximal, distal_ distal to describe the combination of proximal and distal status of Peak1 and Peak2.
* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space
* `subcluster`: combined cluster ID (cluster + subset cluster)

 |  
| [master_cicero_conns.rds][19] |  01/11/2018 |  286.3MB |  Same as above in RDS format. |  
| [open_sites_0.01.txt][20] |  01/12/2018 |  455.3MB |  Peaks that served as input to Cicero in TSV format (accessible in >=1% of cells in a given cluster). [format details][21]

* `site`: ID of peak (chr_start_end)
* `num_cells_expressed`: number of cells in cluster with at least one read overlapping this site
* `fraction_cells_expressed`: fraction of cells in cluster with at least one read overlapping this site
* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space

 |  
| [open_sites_0.01.rds][22] |  01/11/2018 |  91.5MB |  Same as above in RDS format. | 

## Basset Results

We have also trained convolutional neural network (CNN) models with Basset ([Github][23]; [Kelley _et al._][24]) to find motifs that distinguish our clusters from one another. Here we provide results relevant to interpretation of these motifs as well as the actual models generated by Basset so they may be used in downstream analyses. To interpret the "filters" in the first layer of the CNN model, we utilize common tools for interpreting PWMs such as TomTom and MEME from the [MEME suite][25] in conjunction with the [Hocomoco PWM database][26]. 

| Name                        | Last modified | Size   | Description                                                                                             |  
| --------------------------- | ------------- | ------ | ------------------------------------------------------------------------------------------------------- |  
| [training_params.txt][27]   | 01/11/2018    | 286.0B | Parameters used to train the CNN using Basset. See Basset documentation.                                |  
| [trained_model.th][28]      | 01/11/2018    | 88.0MB | CNN model trained using Basset. See Basset documentation.                                               |  
| [train_test_data.h5][29]    | 01/11/2018    | 4.3GB  | Data used to train and test the CNN in the `hdf5` format. See Bassett documentation.                    |  
|                             |               |        |                                                                                                         |  
| [filter_influences.txt][30] | 01/11/2018    | 1.2MB  | Filter influences as calculated using the CNN model in our manuscript in TSV format. [format details][31]  

* `filter`: name of filter PWM
* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space
* `infl`: Influence score. See manuscript.
* `filter_id`: Numeric ID for the filter

 |  
| [filter_sd_mean.txt][32] |  01/11/2018 |  16.5KB |  Contains summary statistics for filter activity over the test data. [format details][33]

* `filter`: name of filter PWM
* `ic`: information content of the filter over all clusters
* `mean`: mean influence of the filter over all clusters
* `sd`: standard deviation for influence of the filter over all clusters

 |  
| [filters_meme.txt][34] |  01/11/2018 |  343.2KB |  PWMs for all the first layer filters in the trained CNN (used for interpretation of filters). See MEME documentation. |  
| [tomtom_hits.txt][35] |  01/11/2018 |  247.8KB |  Filters matched to known motifs using TomTom with value format details 

* `#Query ID`: name of filter PWM
* `Target ID`: motif name in Hocomoco database
* `Optimal offset`: See TomTom documentation
* `Overlap`: See TomTom documentation
* `Query consensus sequence`: See TomTom documentation
* `Target consensus sequence`: See TomTom documentation
* `Orientation`: See TomTom documentation

 |  
|  |   |   |   |  
| [cells_by_motifs.txt][36] |  01/11/2018 |  3.4GB |  Aggregated motif scores for cells. [format details][37]

* `filter`: name of filter PWM used to scan accessible sites using TomTom
* `cell`: cell barcode (combined and corrected)
* `motif_score`: the total number of sites that had this matched this PWM
* `total_motifs`: total number of motifs detected in all the accessible sites for this cell
* `motif_activity`: shows the relative frequency of this motif relative with all the motifs scored for this cell (motif score/total motifs)

 |  
| [cells_by_motifs.rds][38] |  01/11/2018 |  674.4MB |  Same as above in RDS format. | 

## Annotation Enrichments

To aid in interpretation we have found it helpful to calculate gene set enrichments using peaks that are DA and have positive betas (are open). We report these enrichments for a number of different gene sets.

### Specificity Score Format

The following columns are provided to report annotation enrichments for each cluster:

| Column             | Description                                                       |  
| ------------------ | ----------------------------------------------------------------- |  
| `cluster`          | Cluster assignment in initial t-SNE.                              |  
| `subset_cluster`   | Cluster assignment in the iterative t-SNE space.                  |  
| `term`             | Term for the gene set reported in the original gene set.          |  
| `name`             | A cleaned version of the `term` column used in figures.           |  
| `p.value`          | P-value of the reported enrichment.                               |  
| `Adjusted.p.value` | P-value of the reported enrichment adjusted for multiple testing. |  
| `fold_change`      | Fold change of the reported enrichment.                           |  
| `gene_coverage`    | Fraction of the gene set covered in enrichment test.              |  

## GWAS _h_2 Enrichments

As described in the manuscript, we report enrichments in heritability (_h_2) in DA peaks with positive betas for each cluster across many human traits as measured by GWAS.

These enrichments are calculated using a tool called partitioned LD score regression (LDSC; [Finucane _et al._][39]; [Github][40]). We also report the trained LDSC models and baseline model which could be used to calculate enrichments for any other trait given the appropriate summary statistics.

| Name                       | Last modified | Size    | Description                                                                                                                                       |  
| -------------------------- | ------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |  
| [gwas_metadata.xlsx][41]   | 08/02/2018    | 27.2KB  | A table with information about all GWAS (other than UKBB) used in our analysis and where summary statistics can be accessed. [format details][42] |  
| [gwas.all_results.txt][43] | 08/02/2018    | 538.6KB | Heritability enrichments for GWAS examined in manuscript (not UKBB) for each cluster in TSV format. [format details][44]  

* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space
* `gwas`: the phenotype studied in GWAS
* `h2`: the heritability estimated for the trait in question (`gwas` column)
* `h2_se`: standard error of `h2`
* `snp_count`: The total SNPs reported to be used by LDSC in log files.
* `prop_snps`: The proportion of `snp_count` that fall into DA peaks for the specified cluster
* `prop_h2`: The proportion of total `h2` that the DA peaks account for.
* `coefficient`: The coefficient reported by the LDSC model for the DA peak membership term.
* `scaled_coefficient`: `coefficient` column divided by the per SNP heritability. Roughly interpretable as an enrichment
* `enrichment`: an enrichment calculated independent of the baseline enrichment
* `coefficient_zscore`: the z-score reported for the unscaled `coefficient`
* `coefficient_std_error`: standard error for `coefficient`
* `coefficient_pval`: P-value for the coefficient being non-zero (enriched over baseline)
* `coefficient_qval`: `pvalue` adjusted for multiple testing

 |  
| [gwas.clustered_matrix.txt][45] |  08/02/2018 |  8.0KB |  A matrix of -log10(qval) for any significant enrichments from above file in TSV format. [format details][46]

Each row is a trait (as indicated by `trait` column) and each additional column is a cluster of cells. Entries appear in order according to the hierarchical clustering that appears in the manuscript. Zero entries indicate enrichments that are not significant.

 |  
|  |   |   |   |  
| [gwas.tissues.all_results.txt][47] |  08/06/2018 |  100.1KB |  Heritability enrichments for GWAS examined in manuscript (not UKBB) for peaks called from each tissue in TSV format. [format details][48]

* `tissue`: the tissue (replicates are not collapsed)
* `gwas`: the phenotype studied in GWAS
* `h2`: the heritability estimated for the trait in question (`gwas` column)
* `h2_se`: standard error of `h2`
* `snp_count`: The total SNPs reported to be used by LDSC in log files.
* `prop_snps`: The proportion of `snp_count` that fall into DA peaks for the specified cluster
* `prop_h2`: The proportion of total `h2` that the DA peaks account for.
* `coefficient`: The coefficient reported by the LDSC model for the DA peak membership term.
* `scaled_coefficient`: `coefficient` column divided by the per SNP heritability. Roughly interpretable as an enrichment
* `enrichment`: an enrichment calculated independent of the baseline enrichment
* `coefficient_zscore`: the z-score reported for the unscaled `coefficient`
* `coefficient_std_error`: standard error for `coefficient`
* `coefficient_pval`: P-value for the coefficient being non-zero (enriched over baseline)
* `coefficient_qval`: `pvalue` adjusted for multiple testing

 |  
| [gwas.tissues.clustered_matrix.txt][49] |  08/06/2018 |  9.6KB |  A matrix of -log10(qval) for any significant enrichments from above file in TSV format. [format details][50]

Each row is a trait (as indicated by `trait` column) and each additional column a tissue. Entries appear in order according to the hierarchical clustering that appears in the manuscript. Zero entries indicate enrichments that are not significant.

 |  
|  |   |   |   |  
| [ukbb.all_results.txt][51] |  08/02/2018 |  9.4MB |  Heritability enrichments for UKBB GWAS examined in manuscript for each cluster in TSV format. [format details][52]

* `cluster`: cluster assignment in initial t-SNE
* `subset_cluster`: cluster assignment in iterative t-SNE space
* `field`: the phenotype description provided by the Neale Lab
* `field_cleaned`: a (sometimes) shortened version of `field` to remove redundant text
* `field_code`: the non-descriptive ID for each `field` value
* `effective_n`: the effective sample size for this trait
* `h2`: the heritability estimated for the trait in question (`gwas` column)
* `h2_se`: standard error of `h2`
* `snp_count`: The total SNPs reported to be used by LDSC in log files.
* `prop_snps`: The proportion of `snp_count` that fall into DA peaks for the specified cluster
* `prop_h2`: The proportion of total `h2` that the DA peaks account for.
* `coefficient`: The coefficient reported by the LDSC model for the DA peak membership term.
* `scaled_coefficient`: `coefficient` column divided by the per SNP heritability. Roughly interpretable as an enrichment
* `enrichment`: an enrichment calculated independent of the baseline enrichment
* `coefficient_zscore`: the z-score reported for the unscaled `coefficient`
* `coefficient_std_error`: standard error for `coefficient`
* `coefficient_pval`: P-value for the coefficient being non-zero (enriched over baseline)
* `coefficient_qval`: `pvalue` adjusted for multiple testing

 |  
| [ukbb.clustered_matrix.txt][53] |  08/02/2018 |  81.9KB |  A matrix of -log10(qval) for any significant enrichments from above file in TSV format. [format details][54]

Each row is a trait (as indicated by `trait` column) and each additional column a cluster in TSV format. Entries appear in order according to the hierarchical clustering that appears in the manuscript. Zero entries indicate enrichments that are not significant.

 |  
| [ukbb.subset.clustered_matrix.txt][55] |  08/02/2018 |  45.0KB |  A matrix of -log10(qval) for a subset of phenotypes (as shown in manuscript) in TSV format. [format details][56]

Each row is a trait (as indicated by `trait` column) and each additional column a cluster. Entries appear in order according to the hierarchical clustering that appears in the manuscript. Zero entries indicate enrichments that are not significant.

 |  
|  |   |   |   |  
| [ld_score_regression_models.tar.gz][57] |  08/02/2018 |  9.7GB |  Set of models for each cluster and a baseline model for comparison trained with LDSC. [format details][58]

When unpacked with `tar -xzvf ld_score_regression_models.tar.gz`, will contain `cluster_models` and `baseline_model` subdirectories. These contain the models trained using DA peaks from each cluster and the baseline model to compare against, respectively.

Note that `cluster_models` will contain one file per cluster per chromosome, so there will be a large number of files in this subdirectory.

These files may then be used in conjuction with the final step of LDSC (see [Github][40] for usage and format descriptions) to calculate enrichments of heritability within the DA peaks for a given cluster for summary statistics from any GWAS.

 | 

## BAM Files

While we provide some raw data on [GEO (GSE111586)][59], we also provide BAM files of the sequences aligned to `mm9` here in case users would like to use them for their own pipelines or methods development. Below we provide one file per tissue, named by their `tissue.replicate` ID as specified in the `tissue.replicate` column of the `cell_metadata.txt` file in the [metadata section][60] above. This means there will be two files for tissues where we performed a replicate and a single file for all other tissues (in addition to BAM index files).

### Important: BAM Details

* Each read is assigned to a cell ID (the sequence specified in the `cell` column of the same metadata file mentioned above). This is encoded in the read name as `cellid:otherinfo`, so the sequence before the colon is the corrected cell barcode sequence for the read.
* Reads are already deduplicated.
* There will be cell IDs that do not appear in our final set of cells, as data is a superset of what ultimately passes our QC steps.
* Files may not download correctly in Chrome (and other web browsers), but they can easily be downloaded with `wget` or `curl`, by right clicking and copying the link address. For example:
    
    
    wget http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/bams/BoneMarrow_62016.bam

## UCSC Trackhub and Bigwigs

We also provide bigWig files and a UCSC trackhub to visualize aggregated pseudo-bulk ATAC-seq profiles for the cells from each cluster. Note that for the smallest clusters, the data will appear fairly sparse even in aggregate at any single locus. In general, we prefer methods for assessing differential accessibility or specificity computationally over visual inspection, although viewing tracks is often useful to get a sense for the data at a given locus.

### UCSC Trackhub

You may access our [UCSC trackhub here][61]. By default the hub will contain a track at the top called `_All_Peak_Calls`, which annotates regions that we called within LSI clusters for each tissue (see Methods). These peaks were used as our features for all downstream analysis.

The trackhub will also contain a track for each cluster in the dataset named according to the convention `cell_label-id`, where `cell_label` and `id` are defined in the same way as they are in our `cell_metadata.txt` file [above][60]. Spaces and periods in cell labels have been removed or replaced as necessary.

### Important

We have generated these tracks using [DeepTools3][62] and `CPM` or `Counts Per Million` normalization with using the `bamCoverage` command with arguments `-bs 1 --normalizeUsing CPM --skipNAs`. The default range displayed on UCSC is 0 to 4 for all tracks, which we find generally works well. However, it is possible that it may need to be adjusted in some cases.

### BigWig Files

In case you would like access to the files used to make the trackhub above, we provide them for download below. Each file is named in the same manner as described above with a `.bw` or `.bb` extension.

[1]: http://atlas.gs.washington.edu/mouse-atac/docs/
[2]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_metadata.txt
[3]: http://atlas.gs.washington.edu#cell-metadata-format
[4]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_metadata.tissue_freq_filtered.txt
[5]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/peak_promoter_intersections.txt
[6]: http://atlas.gs.washington.edu#peak-metadata-format
[7]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_type_assignments.xlsx
[8]: http://atlas.gs.washington.edu#marker-sheet-format
[9]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/scrna_seq_comparisons/knn_results.txt
[10]: http://atlas.gs.washington.edu#knn-format
[11]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/scrna_seq_comparisons/correlation_results.txt
[12]: http://atlas.gs.washington.edu#correlation-scrnaseq-format
[13]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/scrna_seq_comparisons/microwell_truncated_labels.txt
[14]: http://atlas.gs.washington.edu#mca-mapping-format
[15]: https://cole-trapnell-lab.github.io/cicero-release/
[16]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Cicero_Maps/master_cicero_conns.txt
[17]: http://atlas.gs.washington.edu#cicero-map-format
[18]: https://www.biorxiv.org/content/early/2017/07/20/166066
[19]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Cicero_Maps/master_cicero_conns.rds
[20]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Cicero_Maps/open_sites_0.01.txt
[21]: http://atlas.gs.washington.edu#cicero-peak-format
[22]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Cicero_Maps/open_sites_0.01.rds
[23]: https://github.com/davek44/Basset
[24]: http://genome.cshlp.org/content/early/2016/05/03/gr.200535.115.abstract
[25]: http://meme-suite.org/
[26]: http://atlas.gs.washington.edu/hocomoco11.autosome.ru
[27]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/training_params.txt
[28]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/trained_model.th
[29]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/train_test_data.h5
[30]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/filter_influences.txt
[31]: http://atlas.gs.washington.edu#influences-format
[32]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/filter_sd_mean.txt
[33]: http://atlas.gs.washington.edu#filter-stats-format
[34]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/filters_meme.txt
[35]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/tomtom_hits.txt
[36]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/cells_by_motifs.txt
[37]: http://atlas.gs.washington.edu#cell-motifs-format
[38]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/Basset/cells_by_motifs.rds
[39]: https://www.nature.com/articles/ng.3404
[40]: https://github.com/bulik/ldsc
[41]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/gwas_metadata.xlsx
[42]: http://atlas.gs.washington.edu#gwas-metadata-format
[43]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/gwas.all_results.txt
[44]: http://atlas.gs.washington.edu#gwas-format
[45]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/gwas.clustered_matrix.txt
[46]: http://atlas.gs.washington.edu#gwas-matrix-format
[47]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/gwas.tissues.all_results.txt
[48]: http://atlas.gs.washington.edu#gwas-tissues-format
[49]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/gwas.tissues.clustered_matrix.txt
[50]: http://atlas.gs.washington.edu#gwas-tissues-matrix-format
[51]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/ukbb.all_results.txt
[52]: http://atlas.gs.washington.edu#ukbb-format
[53]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/ukbb.clustered_matrix.txt
[54]: http://atlas.gs.washington.edu#ukbb-matrix-format
[55]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/ukbb.subset.clustered_matrix.txt
[56]: http://atlas.gs.washington.edu#ukbb-matrix-subset-format
[57]: http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/gwas_h2_enrichments/ld_score_regression_models.tar.gz
[58]: http://atlas.gs.washington.edu#ldsc_model-format
[59]: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111586
[60]: http://atlas.gs.washington.edu#metadata
[61]: http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&hubUrl=http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/bigwigs/trackhub/hub.txt
[62]: https://deeptools.readthedocs.io/en/develop/

  