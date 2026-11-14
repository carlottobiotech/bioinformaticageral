# 1. Abrir janela para escolher onde criar a pasta do projeto
pasta_base <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::selectDirectory(caption = "Escolha onde criar (ou abrir) a pasta do projeto")
} else {
  choose.dir(caption = "Escolha onde criar (ou abrir) a pasta do projeto")
}

# Interrompe de forma limpa caso clique em Cancelar
if (is.null(pasta_base) || is.na(pasta_base) || pasta_base == "") {
  stop("Nenhuma pasta foi selecionada. Operacao cancelada.")
}

# Evita criar Projeto_scRNA_BALF dentro dela mesma na segunda execucao
if (basename(pasta_base) == "Projeto_scRNA_BALF") {
  dir_projeto <- pasta_base
} else {
  dir_projeto <- file.path(pasta_base, "Projeto_scRNA_BALF")
}

dir.create(file.path(dir_projeto, "dados"), recursive = TRUE, showWarnings = FALSE)
setwd(dir_projeto)

# Conferir o diretorio definido
getwd()

# 2. Instalar apenas pacotes que ainda nao estejam instalados
pacotes <- c("Seurat", "tidyverse", "patchwork", "hdf5r", "writexl")
faltando <- pacotes[!(pacotes %in% installed.packages()[, "Package"])]
if (length(faltando) > 0) {
  install.packages(faltando)
}

# Carregar pacotes
library(Seurat)
library(tidyverse) # ja carrega dplyr e ggplot2
library(patchwork)
library(hdf5r)
library(writexl)

# 3. Listar e importar os ficheiros .h5 da pasta 'dados'
arquivos <- list.files(path = "dados/", pattern = "\\.h5$", full.names = FALSE)
print(arquivos)

if (length(arquivos) == 0) {
  stop(
    "Pasta criada com sucesso, mas ainda esta vazia. Coloque os 7 arquivos .h5 em: ",
    file.path(dir_projeto, "dados"),
    " e rode o bloco novamente."
  )
}

lista_inicial <- list()
nomes_amostras <- c()

for (arq in arquivos) {
  nome <- gsub("\\.h5$", "", arq)
  nomes_amostras <- c(nomes_amostras, nome)
  
  cts <- Read10X_h5(filename = paste0("dados/", arq), use.names = TRUE, unique.features = TRUE)
  if (is.list(cts)) {
    cts <- cts$`Gene Expression`
  }
  
  lista_inicial[[nome]] <- CreateSeuratObject(
    counts = cts,
    project = nome,
    min.cells = 3,
    min.features = 200
  )
}

# 3. Unificacao de Matrizes (Merge) e Metadados
merged_seurat <- merge(
  x = lista_inicial[[1]],
  y = lista_inicial[2:length(lista_inicial)],
  add.cell.ids = nomes_amostras,
  project = "BALF_Pneumonia"
)

barcodes_originais <- rownames(merged_seurat@meta.data)
merged_seurat$sample <- barcodes_originais

merged_seurat@meta.data <- separate(
  merged_seurat@meta.data,
  col = "sample",
  into = c("Patient", "Type", "Barcode"),
  sep = "_"
  
  # 4. Controlo de Qualidade (QC) e Filtragem
  merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
  
  # Conferir se o padrao "^MT-" encontrou os genes mitocondriais
  summary(merged_seurat$mitoPercent)
  
  merged_seurat_filtered <- subset(
    merged_seurat,
    subset = nCount_RNA > 800 &
      nFeature_RNA > 500 &
      nFeature_RNA < 5000 &
      mitoPercent < 10
  )
  
  # Conferir quantos genes e celulas restaram
  merged_seurat_filtered
  
  # 5. Normalizacao e Escalonamento
  merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)
  
  merged_seurat_filtered <- FindVariableFeatures(
    merged_seurat_filtered,
    selection.method = "vst",
    nfeatures = 2000
  )
  
  merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
)

rownames(merged_seurat@meta.data) <- barcodes_originais

# Conferir se a separacao das colunas Patient e Type funcionou
head(merged_seurat@meta.data)
table(merged_seurat$Patient, merged_seurat$Type)

# 6. Reducao de Dimensionalidade e Clusterizacao
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)

merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)

merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)

merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)

# Conferir os clusters gerados (esperado: 0 a 21)
table(Idents(merged_seurat_filtered))

# 7. Quantificacao de Marcadores por Cluster e Exportacao
DefaultAssay(merged_seurat_filtered) <- "RNA"

merged_seurat_filtered <- JoinLayers(merged_seurat_filtered)

marcadores_clusters <- FindAllMarkers(
  merged_seurat_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Exportar tabela completa em CSV
write.csv(marcadores_clusters, file = "top_markers_clusters.csv", row.names = FALSE)

# Criacao do escudo (.gitignore) para evitar upload de ficheiros pesados
writeLines(
  c(
    "dados/*.h5",
    "*.rds",
    "*.RData",
    "*.Rhistory",
    ".Rproj.user/"
  ),
  con = ".gitignore"
)

# Conferir o conteudo do arquivo criado
readLines(".gitignore")

# 8.1 Graficos de Violino: Controle de Qualidade (Bruto vs. Filtrado)
p_vln_bruto <- VlnPlot(
  merged_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Type",
  ncol = 2
)
ggsave(
  "Grafico de Violino BRUTO sem QC.png",
  plot = p_vln_bruto,
  width = 10,
  height = 6,
  dpi = 300
)

p_vln_qc <- VlnPlot(
  merged_seurat_filtered,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Type",
  ncol = 2
)
ggsave(
  "Grafico de Violono de Controle de qualidade.png",
  plot = p_vln_qc,
  width = 10,
  height = 6,
  dpi = 300
)

# Conferir se os dois arquivos foram salvos na pasta do projeto
list.files(pattern = "\\.png$")

# Filtrar os 2 principais marcadores por cluster (44 genes no total, clusters 0 a 21) e exportar para Excel
top2_markers_excel <- marcadores_clusters %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write_xlsx(top2_markers_excel, path = "top_markers_clusters_formatado.xlsx")

# Conferir resultado e arquivos gerados na pasta do projeto
nrow(top2_markers_excel)
head(top2_markers_excel)
list.files(pattern = "top_markers_clusters")

# 8.2 Graficos Base: PCA e UMAP por Condicao e Cluster
p_pca <- DimPlot(merged_seurat_filtered, reduction = "pca", group.by = "Type") +
  ggtitle("PCA - Análise de Componentes Principais")
ggsave("PCA.png", plot = p_pca, width = 8, height = 5, dpi = 300)

p_umap_clusters <- DimPlot(merged_seurat_filtered, reduction = "umap", label = TRUE) +
  ggtitle("Clusters Celulares (BALF)")

p_umap_condicao <- DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Type") +
  ggtitle("Condição Clínica")

ggsave(
  "UMAP_Clusters_Condicao.png",
  plot = p_umap_clusters + p_umap_condicao,
  width = 14,
  height = 6,
  dpi = 300
)

# Conferir os arquivos PNG gerados na pasta do projeto
list.files(pattern = "\\.png$")

# 8.3 Sobreposicao de Expressao no UMAP (FeaturePlot)
p_feat_cd68 <- FeaturePlot(
  merged_seurat_filtered,
  features = "GRCh38------CD68",
  reduction = "umap"
)

ggsave(
  "featureplot sobrepondo clusters com a expressão de genes.png",
  plot = p_feat_cd68,
  width = 8,
  height = 5,
  dpi = 300
)

# Conferir se o arquivo foi salvo
list.files(pattern = "\\.png$")

# 8.4 Perfil de Marcadores Imunitarios por Condicao (DotPlot)
marcadores_chave <- c(
  "GRCh38------CD68",
  "GRCh38------KLRD1",
  "GRCh38------CD8A",
  "GRCh38------CCL19"
)

p_dot <- DotPlot(
  merged_seurat_filtered,
  features = marcadores_chave,
  group.by = "Type"
) +
  RotatedAxis() +
  ggtitle("Expressão de Marcadores por Condição (Dot Plot)")

ggsave("DOTPLOT.png", plot = p_dot, width = 8, height = 5, dpi = 300)

# Conferir se o DOTPLOT.png foi salvo na pasta do projeto
list.files(pattern = "\\.png$")

# 8.5 Ordenacao Sequencial dos Top Marcadores (0 a 21) e Escalonamento
top2_markers <- marcadores_clusters %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

top2_markers$cluster <- factor(top2_markers$cluster, levels = 0:21)

top2_markers <- top2_markers %>% arrange(cluster)

features_top2 <- unique(top2_markers$gene)

merged_seurat_filtered <- ScaleData(
  merged_seurat_filtered,
  features = features_top2
)

# Conferir quantos genes unicos ficaram na ordem dos clusters
length(features_top2)
head(features_top2, 10)

# 8.6 Top Marcadores por Cluster: Heatmaps Alinhados (0 a 10 e 11 a 21)
# Subdivisao limpa removendo niveis vazios para evitar sobreposicao de rotulos nas bordas
sub_0_10 <- subset(merged_seurat_filtered, idents = as.character(0:10))
Idents(sub_0_10) <- droplevels(Idents(sub_0_10))

sub_11_21 <- subset(merged_seurat_filtered, idents = as.character(11:21))
Idents(sub_11_21) <- droplevels(Idents(sub_11_21))

# Parte 1: Clusters 0 a 10
png("Heatmap_Clusters_0_a_10_Final.png", width = 1400, height = 800, res = 120)
print(
  DoHeatmap(sub_0_10, features = features_top2, angle = 0, hjust = 0.5) +
    scale_y_discrete(labels = function(x) gsub("^GRCh38------", "", x)) +
    NoLegend() +
    ggtitle("Top Marcadores - Clusters 0 a 10") +
    theme(
      text = element_text(size = 10),
      axis.text.y = element_text(size = 8, face = "italic")
    )
)
dev.off()

# Parte 2: Clusters 11 a 21 (Continuacao exata na mesma sequencia de genes)
png("Heatmap_Clusters_11_a_21_Final.png", width = 1400, height = 800, res = 120)
print(
  DoHeatmap(sub_11_21, features = features_top2, angle = 0, hjust = 0.5) +
    scale_y_discrete(labels = function(x) gsub("^GRCh38------", "", x)) +
    NoLegend() +
    ggtitle("Top Marcadores - Clusters 11 a 21") +
    theme(
      text = element_text(size = 10),
      axis.text.y = element_text(size = 8, face = "italic")
    )
)
dev.off()

# Conferir todos os arquivos finais gerados na pasta do projeto
list.files(pattern = "\\.(png|csv|xlsx)$")