library(Augur)
library(Seurat)
library(tidyverse)
library(Matrix)
library(ggpubr)
library(SCINA)
library(CellID)
library(harmony)


download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150861/suppl/GSE150861%5Fcell%5Fbatch%2Ecsv%2Egz","data/batch.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150861/suppl/GSE150861%5Fbarcodes%2Etsv%2Egz", "data/cell.csv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150861/suppl/GSE150861%5Ffeatures%2Etsv%2Egz", "data/features.tsv.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150861/suppl/GSE150861%5Fmatrix%2Emtx%2Egz", "data/matrix.mtx.gz")
download.file("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz", "data/panglaodb.tsv.gz")

download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_filtered_gene_bc_matrices.tar.gz", "data/pbmc8k.tar.gz")
untar("data/pbmc8k.tar.gz")
control <- Read10X("filtered_gene_bc_matrices/GRCh38/")
control <- control[rownames(control) %in% features,]

Immune_ref <- read_rds("../CellIDPaperScript/Immune Cells/data/cite_seq/seurat_cbmc_filtered.rds")
batch <- read_csv("data/batch.csv.gz") %>%  rename(cell = X1) %>%  mutate(batch = str_remove(batch, "-rep.*$")) %>% separate(col = batch, into = c("patient", "day")) %>%  mutate(condition = paste0(patient, "_", day)) %>%  column_to_rownames("cell")
features <- read_tsv("data/features.tsv.gz", col_names = FALSE)[[2]]
cell <- read_csv("data/cell.csv.gz", col_names = FALSE)[[1]]
matrix <- readMM("data/matrix.mtx.gz") %>% set_colnames(cell) %>% set_rownames(features)
panglao <- read_tsv("data/panglaodb.tsv.gz")
immune_geneset <- panglao[panglao$organ == 'Immune system',] %>%  group_by(`cell type`) %>%  summarise(geneset = list(`official gene symbol`)) %$%  set_names(geneset,`cell type`) %>%  head(-1)

SeuratAll  <- CreateSeuratObject(matrix, min.cells = 10, meta.data = batch, min.features = 500)
SeuratAll <- FindVariableFeatures(SeuratAll, nfeatures = 3000)
SeuratAll <- NormalizeData(SeuratAll)
SeuratAll <- ScaleData(SeuratAll)
SeuratAll <- RunPCA(SeuratAll, npcs = 50)
SeuratAll <- RunTSNE(SeuratAll, dims = 1:30, num_threads =16)
SeuratAll <- RunUMAP(SeuratAll, dims = 1:30)
SeuratAll <- RunMCA(SeuratAll)
HGTMatrix <- RunCellHGT(SeuratAll, pathways = c(CITE_signatures,cc.genes))
predictions <- rownames(HGTMatrix)[apply(HGTMatrix,2, which.max)]
predictions <- ifelse(apply(HGTMatrix,2, max)>2, predictions, "NA")
SeuratAll$predictions <- predictions
ggarrange(plotlist = lapply(c("day", "patient", "predictions"), function(x) DimPlot(SeuratAll, group.by = x)), labels = "AUTO")
DimPlot(SeuratAll, group.by = "predictions", reduction = "harmonyUMAP", label =T)
SeuratAll <- RunHarmony(SeuratAll, c("day", "patient"))

SeuratAll <- RunTSNE(SeuratAll, dims = 1:30, num_threads =16, reduction = "harmony", reduction.name = "harmonyTSNE")
SeuratAll <- RunUMAP(SeuratAll, dims = 1:30, reduction = "harmony", reduction.name = "harmonyUMAP")
DimPlot(SeuratAll, group.by = "predictions", reduction = "harmonyUMAP", label =T)
FeaturePlot(SeuratAll, "CD34", reduction = "harmonyUMAP", order =T)
stickbug <- SCINA(SeuratAll@assays$RNA@data, signatures = immune_geneset,rm_overlap = F)
SeuratAll$scina <- stickbug$cell_labels

table(SeuratAll$scina)
DimPlot(SeuratAll, group.by = "scina", reduction = "harmonyUMAP", label =T) %>%  plotly::ggplotly()
DimPlot(SeuratAll, group.by = "predictions", reduction = "harmonyUMAP", label =T)

CITE_signatures <- GetGroupGeneSet(Immune_ref)

GO <- fgsea::gmtPathways("https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018")
GO <- lapply(GO, function(x) x[x !=''])
HGTMatrixFunc <- RunCellHGT(SeuratAll, pathways = GO, minSize = 5)
SeuratAll@assays$FUNC <- CreateAssayObject(HGTMatrixFunc)
GO$`regulation of inflammatory response (GO:0050727)`
FeaturePlot(SeuratAll, "regulation of inflammatory response (GO:0050727)", min.cutoff = 2, reduction = "harmonyUMAP", order = T)
FeaturePlot(SeuratAll, "ATF3", min.cutoff = 0, reduction = "harmonyUMAP", order = T)
SeuratAll <- SeuratWrappers::RunFastMNN(SeuratAll, nThreads = 32)
