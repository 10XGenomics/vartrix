## VarTrix

VarTrix is a tool for extracting single-cell variant information from single-cell sequencing datasets. VarTrix uses Smith-Waterman alignment to evaluate reads that map to a known variant locus and assign single cells to these variants. This process works on both single-cell RNA sequencing datasets as well as single-cell DNA sequencing datasets.

At this point, all multi-allelic sites are ignored. They will still be a column in the final matrix to maintain ordering, but all values will be empty.

### Usage

`--vcf (-v)`: Input VCF formatted variants to be assigned. REQUIRED.

`--bam (-b)`: Input CellRanger BAM. This BAM must have the `CB` tag to define the barcodes of cell barcodes. REQUIRED.

`--fasta (-f)`: A FASTA file for the reference genome used in the BAM. Must have a index file. REQUIRED.

`--cell-barcodes (-c)`: A cell barcodes file as produced by CellRanger that defines which barcodes were called as cells. One barcode per line. In CellRanger runs, this can be found in the sub-folder `outs/filtered_gene_bc_matrices_mex/${refGenome}/barcodes.tsv`. REQUIRED.

`--out-matrix (-o)`: The path to write a Market Matrix format matrix out to. This is the same sparse matrix format used by CellRanger, and can be loaded into external tools like Seraut. REQUIRED.

`--out-variants`: The path to write a neat formatting of the variants to for loading into external tools. This file represents the column labels for `--out-matrix` in the format of `$chromosome_$pos`.

`--padding`: The amount of padding around the variant to use when constructing the reference and alternative haplotype for alignment. This should be no shorter than your read length. DEFAULT: 100bp.

`--scoring-method (-s)`: The scoring method to be used in the output matrix. In the default `binary` mode, the matrix will have a `1` if all reads at the position support the ref allele, and a `2` if one or more reads support the alt allele. In the `alt_frac` mode, the output matrix will have the fraction of alternate allele reads seen at this position. In the `coverage` mode, two matrices are produced. The matrix sent to `--out-matrix` is the number of alt reads seen, and the matrix sent to `--ref-matrix` is the number of ref reads seen. DEFAULT: binary.

`--ref-matrix`: If `--scoring-method` is set to `coverage`, this must also be set. This is the path that the reference coverage matrix will be written to.

`--threads`: The number of parallel threads to use.

`--log-level`: One of `info`, `error` or `debug`. Increasing levels of logging. `Debug` mode is extremely verbose and will report on the fate of every single read. DEFAULT: error.

`--mapq`: The minimum mapping quality of reads to be considered. Default: 0.

`--primary-alignments`: Boolean flag -- consider only primary alignments? Default: false.

`--no-duplicates`: Boolean flag -- ignore alignments marked as duplicates? Take care when turning this on with scRNA-seq data, as duplicates are marked in that pipeline for every extra read sharing the same UMI/CB pair, which will result in most variant data being lost. Default: false.


### License
VarTrix is distributed under the MIT license.
