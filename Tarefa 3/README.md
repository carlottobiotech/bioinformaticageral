# Projeto de Análise de scRNA-Seq: COVID-19 vs. *Pseudomonas aeruginosa* (BALF)

Repositório destinado ao pipeline completo de transcriptômica de célula única (Tarefa 3 da disciplina de Bioinformática). O objetivo é caracterizar a heterogeneidade celular e avaliar a resposta imunológica no microambiente alveolar humano frente à infecção viral severa (SARS-CoV-2) e à pneumonia bacteriana (*Pseudomonas aeruginosa*).

## 1. Metadados do Projeto

* **Base de Dados:** NCBI GEO - `GSE155249`
* **Tipo de Dados:** Single-Cell RNA-Seq (*Expression profiling by high throughput sequencing* - 10x Genomics).
* **Organismo:** *Homo sapiens* (humano / genoma de referência híbrido `GRCh38` + `SARS-CoV-2`).
* **Desenho Experimental:** 3 condições clínicas analisadas a partir de amostras de fluido de lavagem broncoalveolar (BALF), totalizando 7 matrizes de expressão pré-processadas (`.h5`): Controle Saudável (`PacA_Controle`), Pneumonia por *Pseudomonas aeruginosa* (`PacB_Pseudomonas`) e COVID-19 (`Pac1_COVID`, `Pac3A_COVID`, `Pac3B_COVID`, `Pac4_COVID`, `Pac5_COVID`).
* **Amostras:** Matrizes de contagem brutas no formato HDF5 (`.h5`) baixadas via NCBI GEO e processadas localmente em ambiente R/RStudio (Seurat v5). Os arquivos `.h5` não são versionados neste repositório devido ao tamanho (bloqueados via `.gitignore`).

## 2. Preparação do Ambiente e Aquisição de Dados

Seleção interativa do diretório do projeto, criação automática da pasta `dados/`, instalação condicional dos pacotes necessários no RStudio e leitura automatizada das matrizes `.h5`:

```r
# 1. Abrir janela para escolher onde criar (ou abrir) a pasta do projeto
pasta_base <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  rstudioapi::selectDirectory(caption = "Escolha onde criar (ou abrir) a pasta do projeto")
} else {
  choose.dir(caption = "Escolha onde criar (ou abrir) a pasta do projeto")
}

if (is.null(pasta_base) || is.na(pasta_base) || pasta_base == "") {
  stop("Nenhuma pasta foi selecionada. Operacao cancelada.")
}

if (basename(pasta_base) == "Projeto_scRNA_BALF") {
  dir_projeto <- pasta_base
} else {
  dir_projeto <- file.path(pasta_base, "Projeto_scRNA_BALF")
}

dir.create(file.path(dir_projeto, "dados"), recursive = TRUE, showWarnings = FALSE)
setwd(dir_projeto)

# 2. Instalar apenas pacotes que ainda nao estejam instalados
pacotes <- c("Seurat", "tidyverse", "patchwork", "hdf5r", "writexl")
faltando <- pacotes[!(pacotes %in% installed.packages()[, "Package"])]
if (length(faltando) > 0) {
  install.packages(faltando)
}

# Carregar pacotes
library(Seurat)
library(tidyverse) # inclui dplyr e ggplot2
library(patchwork)
library(hdf5r)
library(writexl)

# 3. Listar e importar os arquivos .h5 da pasta 'dados'
arquivos <- list.files(path = "dados/", pattern = "\\.h5\$", full.names = FALSE)

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
  nome <- gsub("\\.h5\$", "", arq)
  nomes_amostras <- c(nomes_amostras, nome)
  
  cts <- Read10X_h5(filename = paste0("dados/", arq), use.names = TRUE, unique.features = TRUE)
  if (is.list(cts)) {
    cts <- cts\$`Gene Expression`
  }
  
  lista_inicial[[nome]] <- CreateSeuratObject(
    counts = cts,
    project = nome,
    min.cells = 3,
    min.features = 200
  )
}
```

## 3. Unificação de Matrizes (Merge) e Metadados

Fusão dos objetos individuais em um objeto Seurat único (`merged_seurat`) e extração dos metadados clínicos (`Patient`, `Type` e `Barcode`), preservando os identificadores originais das células:

```r
merged_seurat <- merge(
  x = lista_inicial[[1]],
  y = lista_inicial[2:length(lista_inicial)],
  add.cell.ids = nomes_amostras,
  project = "BALF_Pneumonia"
)

barcodes_originais <- rownames(merged_seurat@meta.data)
merged_seurat\$sample <- barcodes_originais

merged_seurat@meta.data <- separate(
  merged_seurat@meta.data,
  col = "sample",
  into = c("Patient", "Type", "Barcode"),
  sep = "_"
)

rownames(merged_seurat@meta.data) <- barcodes_originais
```

## 4. Controle de Qualidade (QC) e Filtragem

Cálculo da porcentagem de transcritos mitocondriais (`^MT-`) e remoção de gotículas vazias, células apoptóticas e *doublets*, retendo **35.729 células** de alta qualidade e **20.002 genes**:

```r
merged_seurat\$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

merged_seurat_filtered <- subset(
  merged_seurat,
  subset = nCount_RNA > 800 &
           nFeature_RNA > 500 &
           nFeature_RNA < 5000 &
           mitoPercent < 10
)
```

## 5. Normalização e Escalonamento

Normalização logarítmica das contagens, identificação dos 2.000 genes altamente variáveis (`vst`) e escalonamento linear dos dados:

```r
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)

merged_seurat_filtered <- FindVariableFeatures(
  merged_seurat_filtered,
  selection.method = "vst",
  nfeatures = 2000
)

merged_seurat_filtered <- ScaleData(merged_seurat_filtered)
```

## 6. Redução de Dimensionalidade e Clusterização

Execução da Análise de Componentes Principais (PCA), construção do grafo de vizinhança (KNN), agrupamento celular não supervisionado (**22 clusters: 0 a 21**) e projeção não linear UMAP:

```r
merged_seurat_filtered <- RunPCA(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:20)
```

## 7. Quantificação de Marcadores por Cluster e Exportação

Unificação das camadas de expressão do Seurat v5 (`JoinLayers`), detecção global de biomarcadores positivos para cada população celular (`FindAllMarkers`) e exportação da tabela completa (`.csv`) e do resumo com o Top 2 marcadores por cluster (`top_markers_clusters_formatado.xlsx`):

```r
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

# Filtrar os 2 principais marcadores por cluster (44 genes no total, clusters 0 a 21) e exportar para Excel
top2_markers_excel <- marcadores_clusters %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write_xlsx(top2_markers_excel, path = "top_markers_clusters_formatado.xlsx")
```

## 8. Análise de Expressão Diferencial e Visualização de Subpopulações

A modelagem estatística de marcadores, o mapeamento topológico e a geração dos gráficos finais foram realizados em R utilizando os pacotes `Seurat`, `ggplot2`, `patchwork` e `dplyr`.

### Visualizações Geradas

* **Violin Plots (Bruto vs. Pós-QC):** Avaliação da distribuição de genes (`nFeature_RNA`) e contagens (`nCount_RNA`), confirmando a remoção de *doublets* (>6.000 genes na COVID-19) e ruído basal.
* **PCA Plot:** Verificação da variância linear entre Controle, COVID-19 e *Pseudomonas*, evidenciando um eixo inflamatório exclusivo nas infecções agudas.
* **UMAP & FeaturePlot (`CD68`):** Projeção topológica dos 22 clusters celulares e sobreposição da expressão de `CD68`, delimitando o compartimento majoritário de macrófagos alveolares frente aos linfócitos infiltrantes.
* **Dot Plot por Condição:** Comparação quantitativa de marcadores-chave (`CD68`, `KLRD1`, `CD8A`, `CCL19`), demonstrando a hegemonia mieloide no Controle (~80%) e o pico citotóxico antiviral (T CD8+ e NK) na COVID-19.
* **Heatmaps Divididos (Clusters 0–10 e 11–21):** Agrupamento do Top 2 marcadores por cluster. Como a plotagem única dos 22 clusters gerava forte compressão horizontal das colunas menores e sobreposição ilegível de textos e números, o heatmap foi dividido em dois painéis sequenciais alinhados na mesma ordem gênica (0 a 10 e 11 a 21), destacando desde subpopulações de macrófagos (`PLTP`, `SPP1`, `FABP4`) até o cluster 11 enriquecido em transcritos virais (`SARS-CoV-2i-S`, `SARS-CoV-2i-ORF8`).

### Proteção do Repositório (`.gitignore`)

Criação do arquivo `.gitignore` para evitar o upload de arquivos pesados (`.h5` e `.rds`) e arquivos temporários do RStudio:

```r
# Opcao via R
writeLines(
  c("dados/*.h5", "*.rds", "*.RData", "*.Rhistory", ".Rproj.user/"),
  con = ".gitignore"
)
```

```bash
# Opcao via Terminal (Git Bash / Linux / macOS)
cat << 'EOF' > .gitignore
dados/*.h5
*.rds
*.RData
*.Rhistory
.Rproj.user/
EOF
```

### Scripts de Geração dos Gráficos Finais

#### 8.1 Gráficos de Violino: Controle de Qualidade (Bruto vs. Filtrado)

```r
p_vln_bruto <- VlnPlot(
  merged_seurat,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Type",
  ncol = 2
)
ggsave("Grafico de Violino BRUTO sem QC.png", plot = p_vln_bruto, width = 10, height = 6, dpi = 300)

p_vln_qc <- VlnPlot(
  merged_seurat_filtered,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "Type",
  ncol = 2
)
ggsave("Grafico de Violino de Controle de qualidade.png", plot = p_vln_qc, width = 10, height = 6, dpi = 300)
```

#### 8.2 Gráficos Base: PCA e UMAP por Condição e Cluster

```r
p_pca <- DimPlot(merged_seurat_filtered, reduction = "pca", group.by = "Type") +
  ggtitle("PCA - Análise de Componentes Principais")
ggsave("PCA.png", plot = p_pca, width = 8, height = 5, dpi = 300)

p_umap_clusters <- DimPlot(merged_seurat_filtered, reduction = "umap", label = TRUE) +
  ggtitle("Clusters Celulares (BALF)")

p_umap_condicao <- DimPlot(merged_seurat_filtered, reduction = "umap", group.by = "Type") +
  ggtitle("Condição Clínica")

ggsave("UMAP_Clusters_Condicao.png", plot = p_umap_clusters + p_umap_condicao, width = 14, height = 6, dpi = 300)
```

#### 8.3 Sobreposição de Expressão no UMAP (FeaturePlot)

```r
p_feat_cd68 <- FeaturePlot(
  merged_seurat_filtered,
  features = "GRCh38------CD68",
  reduction = "umap"
)
ggsave("featureplot sobrepondo clusters com a expressão de genes.png", plot = p_feat_cd68, width = 8, height = 5, dpi = 300)
```

#### 8.4 Perfil de Marcadores Imunológicos por Condição (DotPlot)

```r
marcadores_chave <- c(
  "GRCh38------CD68",
  "GRCh38------KLRD1",
  "GRCh38------CD8A",
  "GRCh38------CCL19"
)

p_dot <- DotPlot(merged_seurat_filtered, features = marcadores_chave, group.by = "Type") +
  RotatedAxis() +
  ggtitle("Expressão de Marcadores por Condição (Dot Plot)")
ggsave("DOTPLOT.png", plot = p_dot, width = 8, height = 5, dpi = 300)
```

#### 8.5 Ordenação Sequencial dos Top Marcadores (0 a 21)

```r
top2_markers <- marcadores_clusters %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

top2_markers$cluster <- factor(top2_markers$cluster, levels = 0:21)
top2_markers <- top2_markers %>% arrange(cluster)

features_top2 <- unique(top2_markers\$gene)

merged_seurat_filtered <- ScaleData(merged_seurat_filtered, features = features_top2)
```

#### 8.6 Top Marcadores por Cluster: Heatmaps Alinhados (0 a 10 e 11 a 21)

```r
# Subdivisao limpa removendo niveis vazios para evitar sobreposicao de rotulos nas bordas
sub_0_10 <- subset(merged_seurat_filtered, idents = as.character(0:10))
Idents(sub_0_10) <- droplevels(Idents(sub_0_10))

sub_11_21 <- subset(merged_seurat_filtered, idents = as.character(11:21))
Idents(sub_11_21) <- droplevels(Idents(sub_11_21))

# Parte 1: Clusters 0 a 10 (removendo o prefixo "GRCh38------" nos rotulos do eixo Y)
png("Heatmap_Clusters_0_a_10_Final.png", width = 1400, height = 800, res = 120)
print(
  DoHeatmap(sub_0_10, features = features_top2, angle = 0, hjust = 0.5) +
    scale_y_discrete(labels = function(x) gsub("^GRCh38------", "", x)) +
    NoLegend() +
    ggtitle("Top Marcadores - Clusters 0 a 10") +
    theme(text = element_text(size = 10), axis.text.y = element_text(size = 8, face = "italic"))
)
dev.off()

# Parte 2: Clusters 11 a 21 (continuacao exata na mesma sequencia de genes)
png("Heatmap_Clusters_11_a_21_Final.png", width = 1400, height = 800, res = 120)
print(
  DoHeatmap(sub_11_21, features = features_top2, angle = 0, hjust = 0.5) +
    scale_y_discrete(labels = function(x) gsub("^GRCh38------", "", x)) +
    NoLegend() +
    ggtitle("Top Marcadores - Clusters 11 a 21") +
    theme(text = element_text(size = 10), axis.text.y = element_text(size = 8, face = "italic"))
)
dev.off()
```
