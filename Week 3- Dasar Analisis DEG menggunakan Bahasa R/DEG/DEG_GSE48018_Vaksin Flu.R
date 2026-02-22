# Dataset: GSE48018 (Influenza vaccination vs Baseline)
# Illumina HumanHT-12 V3.0 expression beadchip GPL6947

# INSTALL & LOAD PACKAGE

install.packages("BiocManager")

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
library(GEOquery)
library(limma)

BiocManager::install("illuminaHumanv4.db", ask = FALSE, update = FALSE)
library(AnnotationDbi)
library(illuminaHumanv4.db)

install.packages(c("pheatmap", "ggplot2", "dplyr"))
library(pheatmap)
library(ggplot2)
library(dplyr)

install.packages("umap")
library(umap)

# PENGAMBILAN DATA DARI GEO 

options(timeout = 6000)
gset <- getGEO("GSE48018", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

## gset akan terdiri dari matrix ekspresi gen, dan banyak data lain seperti metadata sampel, metadata probe, dll
## Timeout diatur jadi 6000 detik agar tidak mudah error saat download (defaut = 60)

#PRE-PROCESSING DATA EKSPRESI 

ex <- exprs(gset)

##Sekarang ex hanya berisi data matriks ekspresi gen dari gset 

qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))

## Ini adalah bagian dari transformasi log2

LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)

## Jika TRUE, maka perlu dilakukan transformasi log2

if (LogTransform) {
  # Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# DEFINISI KELOMPOK SAMPEL 

group_info <- pData(gset)[["treatment:ch1"]]

## group_info hanya ambil data metadata sampel kolom ‘treatment:ch1’ dari gset 

groups <- make.names(group_info)

## Spasi diubah jadi titik, dsb, supaya namanya valid

gset$group <- factor(groups)

## factor berfungsi untuk bilang bahwa objek ini (groups) adalah kelompok

nama_grup <- levels(gset$group)
print(nama_grup)

## Levels berfungsi untuk mendeskripsikan nama kelompok yg sudah generate oleh factor

# DESIGN MATRIX (KERANGKA STATISTIK) 

design <- model.matrix(~0 + gset$group)

colnames(design) <- levels(gset$group)

grup_vaksin <- "trivalent.influenza.vaccination"
grup_baseline <- "baseline"

contrast_formula <- paste(grup_vaksin, "-", grup_baseline)

# ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)

fit <- lmFit(ex, design)

contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

# Mengambil hasil akhir DEG

topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)

# ANOTASI NAMA GEN 

probe_ids <- rownames(topTableResults)

# Mapping probe
gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

# Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

## Konsepnya mirip join di SQL


## Boxplot, Density plot, dan UMAP yang akan dibuat dibawah adalah bagian dari QC untuk memastikan bahwa datanya wajar

# BOXPLOT DISTRIBUSI NILAI EKSPRESI

group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

## Data non-cancer normalnya terlihat seragam antar-kelompok


# DENSITY PLOT DISTRIBUSI NILAI EKSPRESI

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )

## Data non-cancer normalnya terlihat tumpah tindah antar-kelompok

# UMAP DISTRIBUSI NILAI EKSPRESI

umap_input <- t(ex)

umap_result <- umap(umap_input)

umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Baseline vs Vaccinated",
    x = "UMAP 1",
    y = "UMAP 2"
  )

## Data non-cancer normalnya tidak terlihat pemisahan antar-kelompok



# VISUALISASI VOLCANO PLOT 

## Dari sini mulai masuk analisis DEG

volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)

volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 0.1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -0.1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG Vaksin Flu")


# VISUALISASI HEATMAP 


## Agar lebih tepat sasaran, hanya ambil 50 gen paling signifikan berdasarkan adj.P.Val
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$PROBEID, ]

gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID,      
  top50$SYMBOL        
)

rownames(mat_heatmap) <- gene_label

# Data cleaning
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(
  Group = gset$group
)

rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",                 
  annotation_col = annotation_col,
  show_colnames = FALSE,         
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

# SAVE
write.csv(topTableResults, "GSE48018_DEG.csv")