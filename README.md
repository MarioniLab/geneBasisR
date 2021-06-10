# geneBasisR

`geneBasis` is an approach to select an optimal gene library (based on scRNA-seq data) as a function of designated number of genes. 
Essentially, we attempt to identify all sources of heterogeneity driving transcriptional relationships between cells, represented in the form of k-NN graph. For greater details on the method, read our paper. 


<p align="center">
  <img src="geneBasis_cartoon.png" width="500">
</p>


### Installation

```
## Install development version
devtools::install_github("MarioniLab/geneBasisR") 
```

### Main usage

Main functions of the package are `gene_search` and `evaluate_library`.

1. `gene_search` takes as inputs scRNA-seq data and number of genes to select.

Requirements for scRNA-seq data format:

a) SingleCellExperiment object, containing assay 'logcounts'. Henceforth this SingleCellExperiment object will be referred to as sce.

b) Rownames of sce correspond to unique gene identifiers.

c) If colData(sce) contains field 'cell', this field should correspond to unique identifiers of cell entries (if not, we use colnames(sce) instead).

2. `evaluate_library` takes as inputs scRNA-seq data (as a SingleCellExperiment objects, using logcounts) and charcter vector of gene names, and estimates the quality of the selected library at cell type, cell and gene levels. Note that this is independent from `gene_search` meaning that you can plug any selection you want and assess how ~complete it is.


Explore vignette and tutorials to get a grasp on the package and its functions.


### Tutorials

1. [Extended vignette of library design and its evaluation for mouse embryo, E8.5](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/ef2d83ae4eaf607c447037dc8981b18d4e7821af/geneBasis_mouseEmbryo_extended.html)

2. Add-on: illustration of **geneBasis** within an individual cell type + suggestion for how to pre-select relevant genes:
[Vignette of library design within brain cells, mouse embryo, E8.5](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/ef2d83ae4eaf607c447037dc8981b18d4e7821af/geneBasis_mouseEmbryo_within_celltype.html)

3. [Vignette of library design for spleen dataset](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/ef2d83ae4eaf607c447037dc8981b18d4e7821af/geneBasis_spleen.html) . Here we introduce how to create working sce object from raw .txt data and introduce a workflow to compare two independent selections.


