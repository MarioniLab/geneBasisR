# geneBasisR

`geneBasisR` is a package that:

a. Selects an optimal targeted gene panel (based on scRNA-seq data) as a function of designated number of genes. 

b. Provides evaluation of any selected gene panel on -cell type/-cell and -gene levels. 

For greater details on the method, please read our paper: . Also, explore vignette and tutorials to get a better grasp on the package and its functions.

### Installation

```
## Install development version
library(devtools)
devtools::install_github("MarioniLab/geneBasisR") 
```


### Gene panel selection

`gene_search` is the main function of the package and it selects the gene panel of designated size. The schematic below illustrates the steps of the algorithm.

<p align="center">
  <img src="geneBasis_cartoon.png" width="500">
</p>

Essential to specify arguments of `gene_search` are counts matrix (stored in SingleCellExperiment object or henceforth *sce*) and *n_genes_total* specifying the size of the panel. 

Requirements for scRNA-seq data format:

a. SingleCellExperiment object, containing assay 'logcounts'. Henceforth this SingleCellExperiment object will be referred to as sce.

b. Rownames of sce correspond to unique gene identifiers.

c. If colData(sce) contains field 'cell', this field should correspond to unique identifiers of cell entries (if not, we use colnames(sce) for cell IDs instead).

Few notes:

a) If the initial gene screening (e.g. HVG selection) has not been performed, use `retain_informative_genes` prior to `gene_search`.

b) `gene_search` works in iterative fashion and adds genes one by one. The practicality of this is if the initially chosen n_genes_total returned the selection that seems to be insufficient, the selected panel can be plugged back in (specified *genes_base*) to avoid the repetition and discover additional to the selection genes.


```
library(geneBasisR)

# sce - SingleCellExperiment object, wehre normalized counts stored in 'logcounts' assay
# discard definetely uninteresting genes
sce = retain_informative_genes(sce)

# run gene selection
genes = gene_search(sce, n_genes_total = 50)

```

#### Gene panel evaluation

We evaluate gene panels on next levels:

- cell type: for each cell type, we estimate how often cells from the cell type are assigned with the correct cell type based on their neighbors in the 'selection' graph.

- cell: for each cell, we compare normalized distances between neighbors in 'true' and 'selection' graphs. 

- gene: for each gene, we assess imputation accuracy based on the average expression values across cell's neighbors in the 'selection' graph.

The wraper function that performs evaluation is `evaluate_library` takes as inputs scRNA-seq data (as a SingleCellExperiment objects, using logcounts) and character vector of gene names.

Few notes:

a) This function is independent from `gene_search` meaning that you can plug any selection you want and assess how ~complete it is.

b) For assessment of accuracy of cell type mappings, you need to provide cell type labels. We require that it is stored in colData(sce). The default name for the field is 'celltype' - in case, it differes, please specify it in the argument *celltype.id*.


### Tutorials

1. [Extended vignette of library design and its evaluation for mouse embryo, E8.5](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/ef2d83ae4eaf607c447037dc8981b18d4e7821af/geneBasis_mouseEmbryo_extended.html)

2. Add-on: illustration of **geneBasis** within an individual cell type + suggestion for how to pre-select relevant genes:
[Vignette of library design within brain cells, mouse embryo, E8.5](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/03c70e494ae2f36136f6f5cfefdb60e7d5f76fe4/geneBasis_mouseEmbryo_within_celltype.html)

3. [Vignette of library design for spleen dataset](https://rawcdn.githack.com/MarioniLab/geneBasis_tutorials/ef2d83ae4eaf607c447037dc8981b18d4e7821af/geneBasis_spleen.html) . Here we introduce how to create working sce object from raw .txt data and introduce a workflow to compare two independent selections.


