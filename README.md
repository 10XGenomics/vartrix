# VarTrix

VarTrix is a tool for extracting single-cell variant information from single-cell sequencing datasets. VarTrix is **not** a variant calling tool. Users must provide a variant call set that they have previously generated through other means. VarTrix is useful for evaluating heterogeneity within a single-cell sample, which means that the types of variants that will be useful to a user are either *somatic* or *contained within a copy number variant (CNV) event*.

VarTrix uses Smith-Waterman alignment to evaluate reads that map to a known variant locus and assign single cells to these variants. This process works on both single-cell RNA sequencing datasets as well as single-cell DNA sequencing datasets.

VarTrix works with any properly formatted sequence resolved VCF. VarTrix works with SNVs, insertions and deletions.

At this point, all multi-allelic sites are ignored. They will still be a column in the final matrix to maintain ordering relative to the input VCF, but all values will be empty.

## Use cases
VarTrix is useful for evaluating heterogeneity of single cell datasets, including tumor samples and cell lines. 

### scRNA-seq
Allele specific expression in tumor samples can lead to strong correlations between the presence of a specific expressed variant and expression-based clustering. Overlaying the variant information with expression clustering can lead to new insight in to the cause of a specific cancer and to the accumulation of mutations that may lead to relapse or drug resistance.

#### Input variants
Generating variants from scRNA-seq datasets is challenging. Noise inherent in reverse transcription leads to a high false positive rate. We recommend looking at the Broad Institute's GATK and Mutect2 best practices guide for (calling variants in RNAseq)[https://software.broadinstitute.org/gatk/documentation/article.php?id=3891]. An alternative approach is to determine somatic variants using WGS data generated from the same sample as the scRNA-seq library.

### scDNA-seq
Assignment of variants in scDNA data can improve understanding of tumor and cell line heterogeneity. Copy number expansion in tumor cells or chromothripsis in cell lines can lead to different allele fractions of germline variants being associated with subclonal populations. Somatic variants in tumor cells can be associated with subclonal populations and associated with subclones that lead to relapse. Similar to scRNA-seq datasets, variant assignment to specific cells can be overlaid with copy number based clustering.

## Installation

VarTrix has automatically generated downloadable binaries for generic linux and Mac OSX under the [releases page](https://github.com/10XGenomics/vartrix/releases). The linux binaries will work on [any of our supported OSes](https://support.10xgenomics.com/os-support). 

## Compiling from source
VarTrix is standard Rust executable project, that works with stable Rust >=1.13. Install Rust through the standard channels, then type `cargo build --release`. The executable will appear at `target/release/vartrix`. As usual it's important to use a release build to get good performance.

## Usage

`--vcf (-v)`: Input VCF formatted variants to be assigned. REQUIRED.

`--bam (-b)`: Input CellRanger BAM. This BAM must have the `CB` tag to define the barcodes of cell barcodes. Must also have an index file. REQUIRED.

`--fasta (-f)`: A FASTA file for the reference genome used in the BAM. Must have a index file. REQUIRED.

`--cell-barcodes (-c)`: A cell barcodes file as produced by CellRanger that defines which barcodes were called as cells. One barcode per line. In CellRanger runs, this can be found in the sub-folder `outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv` where `${refGenome}` is the name of the reference genome used in your CellRanger run. REQUIRED.

`--out-matrix (-o)`: The path to write a Market Matrix format matrix out to. This is the same sparse matrix format used by CellRanger, and can be loaded into external tools like Seraut. REQUIRED.

`--out-variants`: The path to write a neat formatting of the variants to for loading into external tools. This file represents the column labels for `--out-matrix` in the format of `$chromosome_$pos`.

`--padding`: The amount of padding around the variant to use when constructing the reference and alternative haplotype for alignment. This should be no shorter than your read length. DEFAULT: 100bp.

`--scoring-method (-s)`: The scoring method to be used in the output matrix. In the default `consensus` mode, the matrix will have a `1` if all reads at the position support the ref allele, a `2` if one or more reads support the alt allele, and a `3` if one or more reads support both the alt and the ref allele. In the `alt_frac` mode, the output matrix will have the fraction of alternate allele reads seen at this position. In the `coverage` mode, two matrices are produced. The matrix sent to `--out-matrix` is the number of alt reads seen, and the matrix sent to `--ref-matrix` is the number of ref reads seen. DEFAULT: consensus.

`--ref-matrix`: If `--scoring-method` is set to `coverage`, this must also be set. This is the path that the reference coverage matrix will be written to.

`--threads`: The number of parallel threads to use.

`--log-level`: One of `info`, `error` or `debug`. Increasing levels of logging. `Debug` mode is extremely verbose and will report on the fate of every single read. DEFAULT: error.

`--mapq`: The minimum mapping quality of reads to be considered. Default: 0.

`--primary-alignments`: Boolean flag -- consider only primary alignments? Default: false.

`--no-duplicates`: Boolean flag -- ignore alignments marked as duplicates? Take care when turning this on with scRNA-seq data, as duplicates are marked in that pipeline for every extra read sharing the same UMI/CB pair, which will result in most variant data being lost. Default: false.


## Log level considerations
The default logging level will only report on errors. The next log level, `info`, will report on basic information like the number of variants and barcodes seen, as well as reporting on sites that are problematic in consensus mode. In `debug` mode, the constructed haplotypes and alignments for every single read will be reported. For large datasets, this can produce an extremely large log file.

### Problematic sites
With the log level set to `info` or higher, upon the final scoring step, VarTrix will report on barcode/variant pairs that are inconsistent for potential manual inspection. This situation arises when multiple reads for a given barcode/variant combination have equal alignment scores to both the ref and alt haplotype. The most common cause for this is that this location is a multi-allelic site that was not reported as such in the VCF. This is most often seen in cancer samples with large copy number expansions. In these cases, VarTrix will not consider these reads when populating the matrix.

## Troubleshooting
If you have any bug reports, feature requests, or questions about VarTrix please submit a GitHub issue. We will do our best to respond in a timely fashion and address your concerns. However, VarTrix is not an officially supported 10x Genomics software tool and so we cannot promise that we will be able to address all concerns.

If any uncaught errors happen during execution, VarTrix uses the `human_panic` library which will package the full backtrace into a temporary file. Please put the contents of this file in your GitHub issue.

## Loading data into Seraut

Below is some example code for using the output of VarTrix with Seraut to enable highlighting of variants on expression clusters.

```
library(Seurat)
library(Matrix)
library(stringr)
```

### To get the variants in the correct format
```
vawk '{print $1,$2}' my.vcf > SNV.loci.txt
sed -i 's/\s/:/g' SNV.loci.txt 
```

Where `my.vcf` is the VCF that you used as input to VarTrix.


### Read in the matrix, barcodes, and variants
```
# Read in the sparse genotype matrix
snv_matrix <- readMM("matrix.mtx")

# convert the matrix to a dataframe
snv_matrix <- as.data.frame(as.matrix(snv_matrix))

#read in the cell barcodes output by Cellranger
barcodes <- read.table("filtered_matrix_mex/barcodes.tsv", header = F)

# read in SNV loci
# Should be constructed a single column. For example

# chr1:1234-1235
# chr2:2345-2346

# Construct the final table to add to the Seurat object
snps <- read.table("SNV.loci.txt", header = F)

row.names(snv_matrix) <- barcodes$V1

colnames(snv_matrix) <- snps$V1
```

Where `matrix.mtx` is the output consensus matrix from VarTrix.

Pull a unique variant for from snps (for example chr1:1624866)

```
# Construct the data.frame
gt_chr1_1624866 <- data.frame(gt_chr1_1624866$`chr1:1624866`)

row.names(gt_chr1_1624866) <- barcodes$V1

colnames(gt_chr1_1624866) <- "chr1:1624866"

# Make the encoding more readable
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "0", "No Call")
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "1", "ref/ref")
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "2", "alt/alt")
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "3", "alt/ref")
```

Example output
```
                    chr1:1624866
AAACGAAAGAAATTCG	No Call			
AAACGAAAGCTGGAGT	No Call			
AAACGAAAGTATCTGC	No Call			
AAACGAACACCGAAAG	ref/ref			
AAACGAAGTGCAAGAC	No Call			
AAACGAAGTGTCCAGC	alt/ref
```

```
seurat.data <- Read10X("/path/to/filtered_matrix_mex")
seurat_obj <- CreateSeuratObject(raw.data = seurat.data, min.cells = 1, project = "seurat_obj")

seurat_obj <- AddMetaData(object = seurat_obj, metadata = gt_chr1_1624866)
```

Now process your data to the point of generating a tSNE. If you have never done this before, consult the [Seraut tutorial](https://satijalab.org/seurat/get_started.html).


### Plot the tSNE with variants layered
```
TSNEPlot(object = seurat_obj, do.label = T, colors.use = c("azure2","black", "yellow","red"), 
         pt.size = 2, group.by = "chr1:1624866", label.size = 0.0, plot.order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" ),
         plot.title = "chr1:1624866", do.return = T)
```

## License
VarTrix is licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.
