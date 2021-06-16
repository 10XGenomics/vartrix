# VarTrix

VarTrix is a software tool for extracting single cell variant information from 10x Genomics single cell data. VarTrix will take a set of previously defined variant calls and use that to identify those variants in the single cell data. VarTrix does not perform variant calling. VarTrix is useful for evaluating heterogeneity within a sample, which means that the types of variants that will be useful are either *somatic* or *contained within a copy number variant (CNV) event*. 
 
## Overview of how it works

![overview](https://github.com/10xgenomics/vartrix/blob/master/VarTrix_WorkFlow.png)

VarTrix uses Smith-Waterman alignment to evaluate reads that map to each known input variant locus and assign single cells to these variants. This process works on both 10x single cell gene expression datasets as well as 10x single cell DNA datasets.

VarTrix works with any properly formatted sequence resolved VCF. VarTrix works with SNVs, insertions and deletions. 

At this point, all multi-allelic sites are ignored. They will still be a row in the final matrix to maintain ordering relative to the input VCF, but all values will be empty.

## Use cases
VarTrix is useful for evaluating heterogeneity of 10x single cell datasets, including ones from tumor samples and cell lines. VarTrix can be used to evaluate either *somatic variants*  or *variants contained within a copy number variant (CNV) event*.

### 10x Genomics Single Cell Gene Expression Data
Allele specific expression in tumor samples can lead to strong correlations between the presence of a specific expressed variant and expression-based clustering. Overlaying the variant information with expression clustering can lead to new insights about specific diseases and to the accumulation of mutations that may lead to different phenotypes such as relapse or drug resistance.

### 10x Genomics Single Cell DNA Data
Assignment of variants in scDNA data can improve understanding of tumor and cell line heterogeneity. Copy number expansion in tumor cells or chromothripsis in cell lines can lead to different allele fractions of germline variants being associated with subclonal populations. Somatic variants in tumor cells can be associated with subclonal populations and associated with subclones that lead to relapse. Similar to single cell gene expression datasets, variant assignment to specific cells can be overlaid with copy number based clustering.

## Support
VarTrix is provided as an open-source tool for use by the community. Although we cannot provide full support for the software please submit a [GitHub Issue](https://github.com/10XGenomics/vartrix/issues) if you have any problems, questions or comments. We would also be happy to consider Pull Requests that fix bugs or provide enhancements.

## Installation

VarTrix has automatically generated downloadable binaries for generic linux and Mac OSX under the [releases page](https://github.com/10XGenomics/vartrix/releases). The linux binaries are expected to work on [our supported Operating Systems](https://support.10xgenomics.com/os-support). 

Download the appropriate binary and run

    tar xvzf vartrix-v1.0-x86_64-apple-darwin.tar.gz

or

	tar xvzf vartrix-v1.0-x86_64-linux.tar.gz
	
to decompress it. This will produce a new directory containing the executable binary.  

## Compiling from source
VarTrix is a standard Rust executable project, that works with stable Rust >=1.13. Building from source requires that clang and clang developer libraries
are installed on your system. See [here](https://rust-lang.github.io/rust-bindgen/requirements.html) for details.

If you need to compile from source, [install Rust](https://www.rust-lang.org/en-US/install.html), then type `cargo build --release` from within the directory containing the VarTrix source code. The executable will appear at `target/release/vartrix`. As usual it's important to use a release build to get good performance.

## Testing
If you have compiled VarTrix from source, you can run the tiny test dataset by typing the command `cargo test --release` from wthin the directory containing the VarTrix source code.

The test data files in the `test/` folder also provide a simple example of the inputs and outputs for VarTrix.

## Inputs
VarTrix requires a pre-called variant set in VCF format, an associated set of alignments in BAM or CRAM format, and a genome FASTA file. All sequence names must match between the files. VarTrix also requires a cell barcodes file produced by Cell Ranger, for single cell gene expression data, or Cell Ranger DNA, for single cell DNA data.

### Generating Input Variants
Pre-called variants to be used as input to VarTrix can be generated in many different ways such as gathering calls from existing variant databases or performing variant calling on bulk or single cell genome or transcriptome data. It is important to note that generating variants from bulk or single cell RNA-seq datasets is challenging. Noise inherent in reverse transcription leads to a high false positive rate. We recommend looking at the Broad Institute's GATK and Mutect2 best practices guide for [calling variants in RNAseq](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891). An alternative approach is to determine somatic variants using WGS data generated from the same sample as the scRNA-seq library.

## Outputs
VarTrix produces genome matrices in the same [Market Exchange](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices) format that Cell Ranger uses. This is a sparse matrix format that can be read by common packages. Column labels are the cell barcodes included in the cell barcode input file (specified with `--cell-barcodes`). If the `--out-variants` option is used, VarTrix will produce an additional text file of row names where each variant is named as `$chromosome_$pos`. 

The matrix will contain information about each variant for each cell barcode. The exact output is determined by the [parameters that are used](#usage).  

## Usage
From the directory containing the `vartrix` binary, run VarTrix as:

    ./vartrix -v <path_to_input_vcf> -b <path_to_cellranger_bam> -f <path_to_fasta_file> -c <path_to_cell_barcodes_file> -o <path_for_output_matrix>

### Required arguments 

`--vcf (-v)`: Input VCF formatted variants to be assigned.

`--bam (-b)`: Input Cell Ranger BAM. This BAM must have the `CB` tag to define the barcodes of cell barcodes. Must also have an index file.

`--fasta (-f)`: A FASTA file for the reference genome used in the BAM. Must have a index file.

`--cell-barcodes (-c)`: A cell barcodes file as produced by Cell Ranger that defines which barcodes were called as cells. One barcode per line. In Cell Ranger runs, this can be found in the sub-folder `outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv` where `${refGenome}` is the name of the reference genome used in your Cell Ranger run. This file can be used as column labels for the output matrix.

`--out-matrix (-o)`: The path to write a Market Matrix format matrix out to. This is the same sparse matrix format used by Cell Ranger, and can be loaded into external tools like Seurat.

### Options

`--out-variants`: The path to write a neat formatting of the variants to for loading into external tools. This file represents the row labels for `--out-matrix` in the format of `$chromosome_$pos`.

`--padding`: The amount of padding in base pairs around the variant to use when constructing the reference and alternative haplotype for alignment. This should be no shorter than your read length. DEFAULT: `100`.

`--scoring-method (-s)`: The scoring method to be used in the output matrix. In the default `consensus` mode, the matrix will have a `0` if no reads cover the variant, `1` if all reads at the position support the ref allele, a `2` if one or more reads support the alt allele, and a `3` if one or more reads support both the alt and the ref allele. In the `alt_frac` mode, the output matrix will have the fraction of alternate allele reads seen at this position. In the `coverage` mode, two matrices are produced. The matrix sent to `--out-matrix` is the number of alt reads seen, and the matrix sent to `--ref-matrix` is the number of ref reads seen. DEFAULT: `consensus`.

`--ref-matrix`: If `--scoring-method` is set to `coverage`, this must also be set. This is the path that the reference coverage matrix will be written to.

`--umi`: Boolean flag -- consider UMIs when calculating coverage? In this mode, the BAM must have a `UB` tag reporting the UMI for a given alignment. A consensus of the most common base for a given UMI will be used. This will mostly effect the `coverage` and `alt_frac` modes. DEFAULT: false.

`--threads`: The number of parallel threads to use. DEFAULT: 1.

`--log-level`: One of `info`, `error` or `debug`. Increasing levels of logging. `Debug` mode is extremely verbose and will report on the fate of every single read. DEFAULT: `error`.

`--mapq`: The minimum mapping quality of reads to be considered. Default: `0`.

`--primary-alignments`: Boolean flag -- consider only primary alignments? Default: `false`.

`--no-duplicates`: Boolean flag -- ignore alignments marked as duplicates? Take care when turning this on with scRNA-seq data, as duplicates are marked in that pipeline for every extra read sharing the same UMI/CB pair, which will result in most variant data being lost. Default: `false`.

## Scoring method
In all output modes, each alignment is evaluated by looking at which haplotype provides the higher alignment score. Each alignment is then assigned a value of `{-1, 1, 2}` for being ambiguous, reference supporting and alternate supporting respectively. If the `--umi` flag is set, then a consensus is taken across each UMI for a given cell. Each UMI must exceed a hardcoded 90% threshold for supporting the same allele or that UMI is set to `-1`. Only positive values will be considered when populating the final output table.

In `--umi` mode, the coverage table reports the number of unambiguous UMIs that support reference or alternate allele per cell. 

## Log level considerations
The `--log-level` option sets the amount of information recorded in the output log. 

The default logging level (`error`) will only report on errors. The next log level, `info`, will report on basic information like the number of variants and barcodes seen, as well as reporting on sites that are problematic (see below). In `debug` mode, the constructed haplotypes and alignments for every single read will be reported. For large datasets, this can produce an extremely large log file.

### Problematic sites
With the log level set to `info` or higher, upon the final scoring step, VarTrix will report on barcode/variant pairs that are inconsistent for potential manual inspection. This situation arises when multiple reads for a given barcode/variant combination have equal alignment scores to both the ref and alt haplotype. The most common cause for this is that this location is a multi-allelic site that was not reported as such in the VCF. This is most often seen in cancer samples with large copy number expansions. In these cases, VarTrix will not consider these reads when populating the matrix.

## Troubleshooting
If any uncaught errors happen during execution, VarTrix uses the [human_panic library](https://github.com/rust-clique/human-panic)  which will package the full backtrace into a temporary file. Please include this if you submit a Github issue asking for help to resolve a crash.

## Using VarTrix and Seurat to overlay variant information with gene expression clusters
You can use the [Seurat](https://satijalab.org/seurat/) package to combine with output of VarTrix with other analyses of your single cell expressioon data. 

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

#read in the cell barcodes output by Cell Ranger
barcodes <- read.table("filtered_matrix_mex/barcodes.tsv", header = F)

# read in SNV loci
# Should be constructed a single column. For example

# chr1:1234-1235
# chr2:2345-2346

# Construct the final table to add to the Seurat object
snps <- read.table("SNV.loci.txt", header = F)

colnames(snv_matrix) <- barcodes$V1

row.names(snv_matrix) <- snps$V1
```

Where `matrix.mtx` is the output consensus matrix from VarTrix.

Pull a unique variant for from snps (for example chr1:1624866)

```
# Construct the data.frame
gt_chr1_1624866 <- data.frame(gt_chr1_1624866$`chr1:1624866`)

row.names(gt_chr1_1624866) <- barcodes$V1

colnames(gt_chr1_1624866) <- "chr1:1624866"

# Make the encoding more readable
# See documentation of the consensus scoring mode (`-s consensus`)
# No reads detected
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "0", "No Call")
# Only ref detected
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "1", "ref/ref")
# Only alt detected
gt_chr1_1624866$`chr1:1624866` <- str_replace(as.character(gt_chr1_1624866$`chr1:1624866`), "2", "alt/alt")
# Both alleles detected
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

Now process your data to the point of generating a tSNE. If you have never done this before, consult the [Seurat tutorial](https://satijalab.org/seurat/get_started.html).


### Plot the tSNE with variants layered
```
TSNEPlot(object = seurat_obj, do.label = T, colors.use = c("azure2","black", "yellow","red"), 
         pt.size = 2, group.by = "chr1:1624866", label.size = 0.0, plot.order = c("alt/alt" ,"alt/ref", "ref/ref","No Call" ),
         plot.title = "chr1:1624866", do.return = T)
```

## License
VarTrix is licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.
