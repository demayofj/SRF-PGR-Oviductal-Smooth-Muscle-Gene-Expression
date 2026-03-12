# Serum Response Factor Cooperates with Progesterone Receptor to Regulate Oviductal Smooth Muscle Gene Expression During Preimplantation Embryo Transport

This work identifies SRF as an important regulator of the myosalpinx and described the genes it physically interacts with and modulates in the mouse oviduct during the key stage of embryo transport. We have also shown evidence of cooperation between SRF and PGR in regulating the oviductal smooth muscle genes. These findings provide insight into the gene regulatory requirements for oviductal embryo transport and identify potential contributors to disorders of the human fallopian tube such as ectopic pregnancy.

Below are instructions for reproducing the results seen in this work.

## External tools

The following external tools were used in preparing these results:

* Command-line tools: **Trimmomatic (v0.39), HISAT2 (v2.2.1), SAMTools (v1.18), MACS3 (v3.0.0), HOMER (v4.11), DeepTools (v3.5.4)**
* R packages: **DESeq2**
* Online tools: **CIBERSORTx**

Please install the above tools before proceeding. For help with using the command-line tools (used to create count matrices), please see the [DeMayo Bulk Omics](https://github.com/Petronian/demayo-bulk-omics) repository.

**If you use `conda` or `mamba`, the `environment.yml` file in this repository should specify all the packages you need to run the code here.** You can quickly install all dependencies by executing `conda/mamba env create -n [NAME] -f environment.yml`.

## RNA sequencing

Please see the corresponding section in the manuscript for creating count matrices.

Once the count matrices are created, please see the `notebooks/RNA-Sequencing' for interactive R notebooks demonstrating how to create and analyze DEG lists for the Pgr cre/+ Srf d/d and Myh11 cre Pgr d/d datasets.

## ChIP sequencing

Please see the corresponding section in the manuscript for identifying ChIP peaks.

Once the peaks are identified, please see the `notebooks/ChIP-Sequencing' for interactive R notebooks demonstrating how to analyze these lists for the Pgr cre/+ Srf d/d and Myh11 cre/+ Pgr d/d datasets.

## Deconvolution analysis

Deconvolution analysis was performed using [CIBERSORTx](https://cibersortx.stanford.edu/) to quantitatively infer the expected changes in cell populations based on bulk gene expression. Please visit the CIBERSORTx website for more information about its proper usage.

This step requires two single cell count matrices: a reference matrix and a test matrix. Single-cell RNA-seq data from estrus phase oviducts published by McGlade et al. were utilized as the reference matrix for the deconvolution pipeline since they delineated oviduct cells into muscle, secretory epithelial, ciliated epithelial, and fibroblast populations. The test matrix represents the Pgr cre/+ Srf d/d data analyzed during this experiment.

Please see the `notebooks/Deconvolution` folder for information on how the single-cell datasets were prepared and used within CIBERSORT, as well as how the results were subsequently analyzed.

## Analysis of published transcriptomic datasets

This section requires [Ingenuity Pathway Analysis (IPA)](https://digitalinsights.qiagen.com/products-overview/discovery-insights-portfolio/analysis-and-visualization/qiagen-ipa/) software. Please see the software's website for more information about its use.

IPA takes in DEG lists of various samples and outputs information about differential regulation of gene pathways that are present within the samples represented by the DEG lists. Please see the **RNA sequencing** section above for information about generating DEG lists; an additional list representing Pgr knockout (PRKO) data was downloaded from Akison et al. and is present in Table S1.

Additionally, SRF and PGR expression mapping in human fallopian tube single cell RNA-seq was done using [this online tool](https://www.simplesinglecell.org/winuthayanon/humanFT/); please see the webpage for more information.
