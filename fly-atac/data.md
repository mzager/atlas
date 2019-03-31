---
title: Docs
template: blank
description: A Chromatin Cell Atlas of the Developing Fly Embryo
keywords: SOMA, Brotman Baty
---

[Source](http://atlas.gs.washington.edu/fly-atac/data/ "Permalink to A Chromatin Cell Atlas of the Developing Fly Embryo")

## Published Data

On this page, we make some of the raw data files available for conducting your own analyses.

## Tutorial Data

Here are links to all of the data used in the tutorial (code is presented for downloading and processing these files through the course of the tutorial too). 

| Name                                     | Modified          | Size | Description                                                                                                                                                                                                                             |
| ---------------------------------------- | ----------------- | ---- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [irlba_1.0.3.tar.gz][1]                  | 18-Nov-2017 17:20 | 177K | Source code for `irlba` version 1.0.3, the version used for latent semantic indexing in Cusanovich, Reddington, Garfield _et al._                                                                                                       |
| [monocle_2.3.5.tar.gz][2]                | 18-Nov-2017 17:20 | 16M  | Source code for `Monocle` version 2.3.5, the version used for pseudotemporal ordering in Cusanovich, Reddington, Garfield _et al._                                                                                                      |
| [DDRTree_0.1.4.tar.gz][3]                | 16-Jan-2018 11:23 | 12K  | Source code for `DDRTree` version 0.1.4, the version used for pseudotemporal ordering in Cusanovich, Reddington, Garfield _et al._                                                                                                      |
| [monocle_patch.R][4]                     | 18-Nov-2017 17:20 | 6.6K | A small update to Monocle that causes the `differentialGeneTest` to return beta values.                                                                                                                                                 |
| [6to8.2kbmatrix.sparse.binary.rds][5]    | 18-Nov-2017 17:20 | 113M | A binary matrix recording whether insertions were observed in 2kb windows throughout the genome in all 6-8 hour cells. This file is provided as a sparse matrix in `.rds` format so that it can be easily loaded into `R`.              |
| [6to8.summitmatrix.sparse.binary.rds][6] | 18-Nov-2017 17:19 | 69M  | A binary matrix recording whether insertions were observed in a master list of potential regulatory elements in all 6-8 hour cells. This file is provided as a sparse matrix in `.rds` format so that it can be easily loaded into `R`. |
| [6to8.xycalls.txt][7]                    | 18-Nov-2017 17:19 | 300K | A file assigning a sex to each cell from the 6-8 hour time point on the basis of reads mapping to the sex chromosomes. A `1` indicates male and a `2` indicates female.                                                                 |
| [6to8.read.report.txt][8]                | 18-Nov-2017 17:19 | 429K | A file reporting the total number of unique reads recovered from each cell. This was used in the logistic regression analysis to control for cells having varying complexity.                                                           |
| [2to4_files.tar.gz][9]                   | 18-Nov-2017 17:19 | 44M  | A tarball of all the files required for pseudotemporal ordering of cells from the 2-4 hour time point.                                                                                                                                  |

[1]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/irlba_1.0.3.tar.gz
[2]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/monocle_2.5.3.tar.gz
[3]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/DDRTree_0.1.4.tar.gz
[4]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/monocle_patch.R
[5]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/6to8.2kbmatrix.sparse.binary.rds
[6]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/6to8.summitmatrix.sparse.binary.rds
[7]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/6to8.xycalls.txt
[8]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/6to8.read.report.txt
[9]: http://krishna.gs.washington.edu/content/members/cusanovich/fly_embryogenesis/updated_data/vignette/2to4_files.tar.gz

  