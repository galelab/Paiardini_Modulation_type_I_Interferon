#Run Rv4.1.0 Fruitfly
options(mc.cores = 8)
library(SetRank)
# devtools::install("/share/lwhitmo/GeneSets_SetRank/GeneSets/GeneSets.Homo.sapiens")
library(GeneSets.Homo.sapiens)
library(limma)
library(edgeR)
library(stats)
library(factoextra)
library(umap)
library(Rtsne)
library(ggplot2)
library(ExpressionNormalizationWorkflow)
library(stringr)
library(pvca)
library(mclust)
library(fossil)
library(amap)
library(dendextend)
library(data.table)
library(rrcov)
library("Hmisc")
library(corrplot)
source("./heatmap3LW_function.r")

# TITLE: COVTEN-ANALYSIS FIXED FOR SAMPLE SWAPS
# AUTHOR: LEANNE WHITMORE
######################
# FILES/OBJECTS
######################
rhesus2human <- read.csv(
    file = "rhesus2human.csv",
    header = TRUE,
    stringsAsFactors = FALSE
)

# When set to true code will run enrichment analysis (takes a bit of time so have added this option to easily turn off)
SetRankRun <- TRUE
if (isTRUE(SetRankRun)) {
    # converters
    symbol2EntrezID <- createIDConverter(
        "Homo.sapiens", "SYMBOL",
        "ENTREZID"
    )
    IDConverter <- createIDConverter(
        "Homo.sapiens", "ENTREZID",
        "SYMBOL"
    )
}
####################
# FUNCTIONS
####################

get_upper_tri <- function(cormat) {
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
}

get_lower_tri <- function(cormat) {
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}

theme_minimal_LW <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face = "italic"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

scale_fill_Publication <- function(..) {
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}

scale_colour_Publication <- function(..) {
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ..)
}

generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

vizualize_DE_genes_bp <- function(results, plot_file) {
    print("STATUS: Generating bar plot of number of DE genes...")
    results_t <- t(summary(results))
    results_t <- results_t[, -2]

    for (i in 1:(length(row.names(results_t)))) {
        results_t[i, 1] <- results_t[i, 1] * -1
    }

    DE <- as.data.frame(results_t)
    DE <- setnames(DE,
        old = c("Var1", "Var2", "Freq"),
        new = c("Time_Point", "group", "DE_genes")
    )

    # Create plot
    ggplot(DE, aes(
        x = Time_Point, y = DE_genes, fill = group,
        label = DE$DE_genes
    )) +
        geom_bar(stat = "identity", position = "identity") +
        # geom_text(size = 5, position = position_stack(vjust = 0) )+
        # theme_light() +
        theme_minimal() +
        scale_fill_manual(values = c("#0808c4", "#da9618")) +
        # xlab("Time point")
        ylab("Number of Differentially Expressed Genes") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            axis.text.y = element_text(size = 15)
        )
    ggsave(plot_file, width = 6, height = 4, units="in", dpi = 300)
}

generate_boxplots_voom <- function(data, labels, filename, figres, maintitle, ylabtitle) {
    png(filename, width = 10, height = 8, units = "in", res = figres)
    # par(mar=c(1,1,1,1))
    minvalue <- min(data)
    maxvalue <- max(data)
    boxplot(data,
        labels = labels, ylim = c(minvalue - 1, maxvalue + 1),
        ylab = ylabtitle, main = maintitle, cex.axis = .6, las = 2,
        frame = FALSE
    )
    dev.off()
}

generate_density_plot <- function(data, labels, filename, figres) {
    png(filename, res = figres)
    par(xpd = TRUE)
    if (length(labels) > 10) {
        plotDensities(data, legend = FALSE)
    } else {
        plotDensities(data,
            legend = "topright",
            inset = c(-0.2, 0), levels(labels)
        )
    }
    dev.off()
}

normalize_data <- function(CM2, targetfile, group, control = FALSE) {

    # order target and count matrix so they are the same (THIS IS IMPORTANT)
    CM2 <- CM2[, rownames(targetfile)]

    # CHECK IF ORDER IS THE SAME
    if (all.equal(colnames(CM2), rownames(targetfile)) != TRUE) {
        print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
        print(rownames(targetfile))
        print(colnames(CM2))
    }

    # normalize
    CM2 <- DGEList(counts = CM2)
    CM2 <- calcNormFactors(CM2, method = "TMM") # TMM normalization
    png(file.path(norm_results, "mean_variance_norm.png"))
    Pi.CPM <- voom(counts = CM2, normalize.method = "none", plot = T, span = 0.1, save.plot=T)
    dev.off()
    write.csv(Pi.CPM$E, file.path(norm_results, paste0("1.norm_matrix_", group, ".csv")))

    sig_HGNC <- merge(rhesus2human, Pi.CPM$E,
        by.x = "Gene.stable.ID",
        by.y = "row.names",
        all.X = T, all.Y = T
    )

    sig_HGNC <- sig_HGNC[, !(names(sig_HGNC) %in% c("Gene.stable.ID"))]
    sig_HGNC <- avereps(sig_HGNC,
        ID = sig_HGNC$HGNC.symbol
    )
    rownames(sig_HGNC) <- sig_HGNC[, "HGNC.symbol"]
    sig_HGNC <- sig_HGNC[, !(colnames(sig_HGNC) %in% c("HGNC.symbol"))]
    sig_HGNC <- as.matrix(data.frame(sig_HGNC))
    write.csv(sig_HGNC, file.path(norm_results, paste0("1.norm_matrix_HGNC_", group, ".csv")), quote = FALSE)
    return(Pi.CPM)
}

pca_fun <- function(exprs, labels, results_path,
                    base_file_name, target_columns,
                    figres = 100, size = 1, pca=FALSE) {
    # Run PCA/SVD reduction
    if (isFALSE(pca)) {
        pca <- prcomp(t(exprs))
    }
    E <- get_eig(pca)
    cx <- sweep(t(exprs), 2, colMeans(t(exprs)), "-")
    sv <- svd(cx)


    vizualize_pca(
        file.path(results_path, paste0("svd_", base_file_name)),
        sv$u, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, E, size
    )
    vizualize_pca(
        file.path(results_path, paste0("pca_", base_file_name)),
        pca$x, labels[, target_columns[1]],
        labels[, target_columns[2]],
        figres, E, size
    )
    vizualize_scree_plot(
        file.path(
            results_path,
            paste0("scree_", base_file_name)
        ), pca, figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    is_pc1_0 <- loadingscores$PC1 > 0
    is_pc2_0 <- loadingscores$PC2 > 0

    loadingscores <- loadingscores[is_pc1_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC1)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc1", base_file_name, ".txt")),
        loadingscores["PC1"], figres
    )

    loadingscores <- as.data.frame(pca$rotation)
    loadingscores <- loadingscores[is_pc2_0, ]
    loadingscores <- loadingscores[with(loadingscores, order(-PC2)), ]
    save_loading_scores(
        file.path(results_path, paste0("loadingscores_pc2", base_file_name, ".txt")),
        loadingscores["PC2"], figres
    )
    return(pca)
}

umap_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns,
                     figres = 100, size = 1, UMAP=FALSE) {
    # Runs default paramaters of umap
    if (isFALSE(UMAP)) {
        UMAP <- umap(t(exprs))
    }
    vizualize_umap(
        file.path(results_path, paste0("umap_", base_file_name)),
        UMAP$layout, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size
    )

    return(UMAP)
}

vizualize_umap <- function(plot_file, U, class1, class2, figres, size) {
    # Vizualize umap reduction
    library(Polychrome)
    minx <- min(U[, 1])
    maxx <- max(U[, 1])
    miny <- min(U[, 2])
    maxy <- max(U[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right")
        } else {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right")
        } else if (length(levels(factor(class1))) > 6) {
            qplot(U[, 1], U[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_pca <- function(plot_file, PCA, class1, class2, figres, E, size) {
    # Vizualize PCA  results
    library(Polychrome)
    minx <- min(PCA[, 1])
    maxx <- max(PCA[, 1])
    miny <- min(PCA[, 2])
    maxy <- max(PCA[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right")
        } else {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") + scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36))
        } else if (length(levels(factor(class1))) > 6) {
            qplot(PCA[, 1], PCA[, 2], color = factor(class2), shape = factor(class1), size = I(size)) +
                theme_Publication() +
                theme(legend.title = element_blank()) +
                xlab(paste0("PC1 ", round(E$variance.percent[1], digits = 2), "%")) +
                ylab(paste0("PC2 ", round(E$variance.percent[2], digits = 2), "%")) +
                theme(legend.position = "right") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = 300)
}

tsne_fun <- function(exprs, labels, results_path,
                     base_file_name, target_columns, figres = 300, size = 1, T = FALSE) {
    # Runs default paramaters of umap
    if (isFALSE(T)) {
        T <- Rtsne(t(exprs), perplexity = 1)
    }
    vizualize_tSNE(
        file.path(results_path, paste0("tsne_", base_file_name)),
        T$Y, labels[, target_columns[1]],
        labels[, target_columns[2]], figres, size
    )
    return (T)
}

vizualize_tSNE <- function(plot_file, T, class1, class2, figres, size) {
    # Vizualize tsne reduction
    library(Polychrome)
    minx <- min(T[, 1])
    maxx <- max(T[, 1])
    miny <- min(T[, 2])
    maxy <- max(T[, 2])
    if (length(levels(factor(class2))) <= 3) {
        if (length(levels(factor(class1))) <= 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right")
        } else {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                scale_color_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                scale_fill_manual(values = c("Prot" = "red", "NonProt" = "black", "Uninfected" = "gray")) +
                theme(legend.position = "right") + scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    } else {
        P36 <- createPalette(length(levels(factor(class2))), c("#ff0000", "#00ff00", "#0000ff"))
        if (length(levels(factor(class1))) <= 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlab("tSNE 1") +
                ylab("tSNE 2") +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right")
        } else if (length(levels(factor(class1))) > 6) {
            qplot(T[, 1], T[, 2], shape = factor(paste(class1)), color = factor(class2), size = I(size)) +
                theme_Publication() + theme(legend.title = element_blank()) +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                scale_color_manual(values = as.character(P36)) +
                scale_fill_manual(values = as.character(P36)) +
                xlim(minx, maxx) + ylim(miny, maxy) +
                theme(legend.position = "right") +
                scale_shape_manual(values = seq(1, length(levels(factor(class1)))))
        }
    }
    ggsave(plot_file, width = 6, height = 4, units = "in", dpi = figres)
}

vizualize_scree_plot <- function(plot_file, PCA, figres) {
    # Vizualize principle component variation results
    scree.plot <- fviz_eig(PCA, addlabels = TRUE, hjust = -0.3)
    png(plot_file, width = 7, height = 6, units = "in", res = figres)
    print(scree.plot)
    dev.off()
}

save_loading_scores <- function(write_file, df, figres) {
    # Save list of genes that have a positive effect on variation of principle
    # component 1 and 2 sorted from most influential
    write.table(df, file = write_file)
}

filter_read_counts_mean <- function(cm, filter_cutoff) {
    # Filter value was calculated by:
    #### Filters by row means usually set at 10 reads per gene across all samples

    A <- rowMeans(cm)
    isexpr <- A >= filter_cutoff
    cmfl <- cm[isexpr, ]
    return(cmfl)
}

rename_samples <- function(samples) {
    newsampleIDs <- c()
    for (i in samples) {
        # i <- str_remove(i, "_nohmrRNA_noglobin")
        # i <- str_remove(i, "_RNA1_Lib\\d+_GL\\d+_S\\d+")
        i <- str_remove(i, "_RNA\\d+\\S+")
        i <- str_remove(i, "_RNA\\S+")
        i <- str_remove(i, "G\\d+_CoVTEN-01_")

        newsampleIDs <- c(newsampleIDs, i)
    }
    return(newsampleIDs)
}

generate_design_matrix <- function(normmatrix, target, groups) {
    ti <- factor(target$TimePoint)
    Xid <- factor(target$Animal_ID)
    treat <- factor(target$Sample_Type)
    mm <- model.matrix(~ 0 + ti:treat + Xid)

    rownames(mm) <- colnames(normmatrix)
    colnames(mm) <- make.names(colnames(mm))
    mm <- mm[, colnames(mm)[order(tolower(colnames(mm[, ])))]]
    mm <- mm[, colSums(mm) > 0]

    excludeAll <- nonEstimable(mm)
    if (length(excludeAll) > 0) {
        message("WARNING: These samples are nonEstimatable, design matrix ", excludeAll)
    }

    if ("ti" %in% excludeAll) {
        return("interactions term non estimable")
    }
    mm <- mm[, !colnames(mm) %in% excludeAll]
    if (!is.fullrank(mm)) {
        return("not full rank")
    }
    return(mm)
}

generate_clusterdendograms <- function(hc, plotfilename1, adjvalue, labelsize = 0.7) {
    counter <- 0
    labelsf <- c()
    colors <- c()
    dend <- as.dendrogram(hc)
    dend_labels <- labels(dend)
    P36 <- createPalette(length(levels(factor(dend_labels))), c("#ff0000", "#00ff00", "#0000ff"))
    names(P36) <- unique(dend_labels)
    print(P36)
    for (i in dend_labels) {
        labelsf <- c(labelsf, i)
        colors <- c(colors, P36[[i]])
    }
    labels_colors(dend) <- colors
    labels_cex(dend) <- labelsize
    png(plotfilename1,
        units = "in", # bg = "transparent",
        width = 14.5, height = 5, res = 300
    )
    par(mar = c(6, 3, 2, 0.5), xpd = TRUE)
    plot(dend, xlab = "", main = "")
    mtext(paste0("Adj Rand index ", round(adjvalue, 3)))
    dev.off()
}

vizualize_genes_boxplots <- function(normmatrix, target, genelist, folder) {
    newcolnames <- paste(target$Time_Point, target$Animal_Outcome, sep="_")

    colnames(normmatrix) <- newcolnames
    normmatrix <- melt(normmatrix)

    for (gene in genelist) {
        gene_name <- rhesus2human[rhesus2human$Gene.stable.ID == gene, "HGNC.symbol"]
        if (isTRUE(identical(gene_name, character(0)))) {
            gene_name = gene
        }
        tempnormtmp <- normmatrix[normmatrix$Var1 == gene, ]
        tempnormtmp$Var2 <- factor(tempnormtmp$Var2,
            levels = c(
                "D0_Prot", "D0_NonProt", "D3_Prot", "D3_NonProt",
                "D7_Prot", "D7_NonProt", "D126_Prot", "D126_NonProt",
                "D129_Prot", "D129_NonProt", "D133_Prot", "D133_NonProt",
                "D455_Prot", "D455_NonProt"
            )
        )

        png(file.path(folder, paste0(gene_name, "_boxplot.png")), width = 8, height = 6, units = "in", res = 300)
        par(mar = c(8, 5, 5, 2))
        boxplot(tempnormtmp$value ~ tempnormtmp$Var2,
            main = gene_name, col = (c("red", "#6b6868")),
            xlab = "", ylab = "log2CPM", las = 2, frame = FALSE
        )
        dev.off()
    }
}

gene_enrichment <- function(genes, results_folder, cluster) {
    inputGenesTrans <- rhesus2human[rhesus2human$Gene.stable.ID %in% genes, ]
    inputGenesHGNC <- unique(unlist(inputGenesTrans$HGNC.symbol))
    inputGenes <- symbol2EntrezID(inputGenesHGNC)
    network <- setRankAnalysis(inputGenes, collection,
        use.ranks = FALSE,
        setPCutoff = 0.01,
        fdrCutoff = 0.05
    )

    generate_folder(results_folder)
    #### IMPORTANT OUTPUT INFORMATION###
    # SetRank value -  value reflects the prominence of a gene set in a gene set network (based on PageRank algorithm developed by google)
    # pSetRank value - significance of the SetRank value
    exportSingleResult(network, inputGenes,
        collection, paste0(results_folder, "/de_unranked_", cluster),
        IDConverter = IDConverter
    )
    # png(file.path(results_folder, paste0("de_unranked_", cluster,".png")), res=100)
    # plot(network, layout = layout.spring)
    # dev.off()
    return(network)
}

############################
# Load in data
# countmatrix - load count matricies for all ABL3
# target - load target file for all ABL3
##########################

# --Read in target files
message("STATUS: Load tables")
cm <- read.table("../count_matrix.txt", header = TRUE, sep = "\t", row.names = 1, as.is = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
target <- read.csv("../targetfile_LWedit_fixed4swaps.csv", sep = ",", row.names = 1, as.is = TRUE, check.names = FALSE, header = TRUE, stringsAsFactors = FALSE)

# --Rename sample names in count matrix
newsampleIDs <- rename_samples(colnames(cm))
colnames(cm) <- newsampleIDs
newsampleIDs <- rename_samples(rownames(target))
rownames(target) <- newsampleIDs

if (all.equal(colnames(cm), rownames(target)) != TRUE) {
    print("MASSIVE WARNING: RESULTS WILL BE WRONG IF THIS IS NOT EQUAL!!!!!!!!")
    cm <- cm[, rownames(target)]
    if (all.equal(colnames(cm), rownames(target)) == TRUE) {
        print("ISSUE HAS BEEN FIXED")
    } else {
        print("ISSUE IS NOT FIXED, PLEASE LOOK AT MANUALLY")
    }
}

#--Because of sample swap issues relabel rows in target file and count matrix with real information--#
rownames(target) <- paste(target$Animal_ID, target$Sample_Type, target$TimePoint, sep="_")
colnames(cm) <- rownames(target)

# --Write new renamed files to file
count_results <- "1.count_data/"
generate_folder(count_results)
write.csv(cm, file.path(count_results, "count_matrix_renamed.csv"))
write.csv(target, file.path(count_results, "target_renamed.csv"))

# --generate figure of all counts
generate_density_plot(
    cm, rownames(target), file.path(count_results, "de_intensities_raw_counts.png"), 100
)

# -- Normalize data
norm_results <- "1.norm_data/"
generate_folder(norm_results)
Pi.CPM <- normalize_data(cm, target, "all")
generate_boxplots_voom(Pi.CPM$E, target$Study_Group,
    file.path(norm_results, "boxplot_vnorm_all.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)

#-- filter out genes from each group that are below mean count of 5 across samples 
#Iteratively adjusted thresholds and decided that 5 was the best cutoff to get rid of the 
# mean-variance hook shown in the initial mean-variance plot in 1.norm_results/ (iterative results not saved just ran in R)
cmfl_counts <- filter_read_counts_mean(cm, 5)
write.csv(cmfl_counts, file.path(count_results, "count_matrix_renamed_fl.csv"))

norm_results <- "1.norm_data_fl/"
generate_folder(norm_results)
Pi.CPM <- normalize_data(cmfl_counts, target, "all")
generate_boxplots_voom(Pi.CPM$E, target$Study_Group,
    file.path(norm_results, "boxplot_vnorm_all.png"),
    100,
    maintitle = "Normalized count matrix",
    ylabtitle = "voom normalized expression"
)
saveRDS(Pi.CPM, file.path(norm_results,"normobject.rds"))

########################
# --Feature Reduction
########################
message("STATUS: Run Feature reduction")
feature_results <- "1.feature_red"
generate_folder(feature_results)

pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normtime.png",
    c("Sample_Type", "TimePoint"), 300, 3
)
pca <- pca_fun(
    Pi.CPM$E, target,
    feature_results, "_normAnimalID.png",
    c("Sample_Type", "Animal_ID"), 300, 3, pca=pca
)

umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_normtime.png",
   c("Sample_Type", "TimePoint"), 300, 3
)
umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_normAnimalID.png",
   c("Sample_Type", "Animal_ID"), 300, 3, UMAP=umap
)
umap <- umap_fun(
    Pi.CPM$E, target,
    feature_results, "_normAnimalID-time.png",
    c("TimePoint", "Animal_ID"), 300, 3,
    UMAP = umap
)

tsne <- tsne_fun(
    Pi.CPM$E, target,
    feature_results, "_normtime.png",
    c("Sample_Type", "TimePoint"), 300, 3,
)
tsne <- tsne_fun(
    Pi.CPM$E, target,
    feature_results, "_normAnimalID.png",
    c("Sample_Type", "Animal_ID"), 300, 3,
    T = tsne
)
tsne <- tsne_fun(
    Pi.CPM$E, target,
    feature_results, "_normAnimalID-time.png",
    c("TimePoint", "Animal_ID"), 300, 3,
    T = tsne
)

##########################
# --Hierarchal clustering
##########################
message("STATUS: Running hierarchal clustering...")
hi_cluster_results <- "1.hi_cluster_results"
generate_folder(hi_cluster_results)
d <- Pi.CPM$E
colnames(d) <- target$Animal_ID
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(target$Animal_ID)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(target$Animal_ID))))
generate_clusterdendograms(hc,
    paste0(hi_cluster_results, "/dendogram_AnimalID.png"), m, labelsize = 1)

colnames(d) <- paste(target$Animal_ID, target$TimePoint, sep = "-")
hc <- hcluster(t(d), method = "pearson", link = "average")
x <- cutree(hc, k = length(unique(target$Time_Point)))
m <- adj.rand.index(as.integer(x), as.integer(factor(as.character(target$TimePoint))))
generate_clusterdendograms(
    hc,
    paste0(hi_cluster_results, "/dendogram_normAnimalIDTime.png"), m
)

##############################
# --run lmfit de P. Edelfsen
##############################

mm_all <- generate_design_matrix(Pi.CPM$E, target, groups = c("all"))
Pi.lm <- lapply(rownames(Pi.CPM$E), function(.gene) { ### Here Pi.CPMO is normalized, with weights assigned by voom()
    # print(.gene)
    X <- Pi.CPM$E[.gene, ]
    W <- Pi.CPM$weights[rownames(Pi.CPM$E) == .gene, ]
    lm(X ~ 0 + mm_all, weights = W) ### Okay, this results in the same coefficients as the corresponding row of Pi.lmfit:
    # unname( Pi.lm$coefficients == Pi.lmfit$coefficients[.gene, ] )
    #
    # ### Calculate coefficients and vcov for contrasts
    # .coefs <- Pi.lm$coefficients %*% contr
    # .vcov <- t(contr) %*% vcov(Pi.lm) %*% contr
    # .stdev <- sqrt(diag(.vcov))
    # .tstats <- .coefs/.stdev
    # .pt <- pt(abs(.tstats), lower.tail = F, df = Pi.lm$df.residual) * 2 #2-sided t test
    # array( c(.coefs,.stdev,.tstats,.pt), dim = c(length(.coefs),4), dimnames = list( colnames(contr), c("coef","stdev","t","p") ))
})
names(Pi.lm) <- rownames(Pi.CPM$E)
# Pi.lmfit <- lmFit(Pi.CPMA, design = mm_all)
deresults_path <- "1.de_groups"
generate_folder(deresults_path)

contrastsmatrix <- c(
    "tiD1.treatIFN -tiBL.treatIFN",
    "ti2.treatIFN -tiBL.treatIFN",
    "ti4.treatIFN -tiBL.treatIFN",
    "ti5.treatIFN -tiBL.treatIFN",
    "ti7.treatIFN -tiBL.treatIFN",
    "tiD1.treatCTL -tiBL.treatCTL",
    "ti2.treatCTL -tiBL.treatCTL",
    "ti4.treatCTL -tiBL.treatCTL",
    "ti5.treatCTL -tiBL.treatCTL",
    "ti7.treatCTL -tiBL.treatCTL"
)

#################
# --de analysis
#################
contr <- makeContrasts(contrasts = contrastsmatrix, levels = mm_all)

contrasts.fit.lm <- function(fit, contrasts) {
    out <- sapply(fit, function(.lm) {
        ### Calculate coefficients and vcov for contrasts
        .coefs <- .lm$coefficients %*% contrasts
        .vcov <- t(contrasts) %*% vcov(.lm) %*% contrasts
        .stdev <- sqrt(diag(.vcov))
        .tstats <- .coefs / .stdev
        .pt <- pt(abs(.tstats), lower.tail = F, df = .lm$df.residual) * 2 # 2-sided t test
        array(c(.coefs, .stdev, .tstats, .pt), dim = c(length(.coefs), 4), dimnames = list(colnames(contr), c("coef", "stdev", "t", "p")))
    }, simplify = "array")
}

Pi.contrasts.lm <- contrasts.fit.lm(Pi.lm, contr) # This results in a 3-dimensional array, which is not entirely analogous to Pi.contrasts from above. The next two lines extract the coefficients and p-value matrices from the array, which is all you need to use the decideTests() function.
.coefs <- t(Pi.contrasts.lm[, "coef", ])
.tstats <- t(Pi.contrasts.lm[, "t", ])
.pt <- t(Pi.contrasts.lm[, "p", ])
.pt.adjusted <- apply(.pt, 2, p.adjust, method = "BH")


results.lmde <- decideTests(object = .pt, coefficients = .coefs, lfc = (.58), method = "separate", adjust.method = "BH", p.value = 0.05)
write.csv(.coefs, file = file.path(deresults_path, "coefficients.csv"), quote = F)
write.csv(.tstats, file = file.path(deresults_path, "t_stats.csv"), quote = F)
write.csv(.pt, file = file.path(deresults_path, "p_value.csv"), quote = F)
write.csv(.pt.adjusted, file = file.path(deresults_path, "p_value_adj.csv"), quote = F)


dataMatrixde <- .coefs
sigMask <- dataMatrixde * (results.lmde**2) # 1 if significant, 0 otherwise
ExpressMatrixde <- subset(dataMatrixde, rowSums(sigMask) != 0)

# filter for significant genes - up/down regulated
sigMask <- subset(sigMask, rowSums(sigMask) != 0)
# Pi.contrasts$genes <- data.frame(ID_REF=rownames(Pi.contrasts))

write.csv(ExpressMatrixde, file = file.path(deresults_path, "expression_matrix_de_lm.csv"), quote = F)
write.csv(results.lmde, file = file.path(deresults_path, "results_de_lm.csv"), quote = F)
write.csv(dataMatrixde, file = file.path(deresults_path, "full_expression_matrix_de_lm.csv"), quote = F)
ExpressMatrixde_HGNC <- merge(rhesus2human, ExpressMatrixde,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

ExpressMatrixde_HGNC <- ExpressMatrixde_HGNC[, !(names(ExpressMatrixde_HGNC) %in% c("Gene.stable.ID"))]
ExpressMatrixde_HGNC <- avereps(ExpressMatrixde_HGNC,
    ID = ExpressMatrixde_HGNC$HGNC.symbol
)
rownames(ExpressMatrixde_HGNC) <- ExpressMatrixde_HGNC[, "HGNC.symbol"]
ExpressMatrixde_HGNC <- ExpressMatrixde_HGNC[, !(colnames(ExpressMatrixde_HGNC) %in% c("HGNC.symbol"))]
ExpressMatrixde_HGNC <- as.matrix(data.frame(ExpressMatrixde_HGNC, check.names = FALSE ))
class(ExpressMatrixde_HGNC) <- "numeric"
write.csv(ExpressMatrixde_HGNC, file.path(deresults_path, "expression_matrix_de_lm_hgnc.csv"), quote = F)

dataMatrixde_HGNC <- merge(rhesus2human, dataMatrixde,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

dataMatrixde_HGNC <- dataMatrixde_HGNC[, !(names(dataMatrixde_HGNC) %in% c("Gene.stable.ID"))]
dataMatrixde_HGNC <- avereps(dataMatrixde_HGNC,
    ID = dataMatrixde_HGNC$HGNC.symbol
)
rownames(dataMatrixde_HGNC) <- dataMatrixde_HGNC[, "HGNC.symbol"]
dataMatrixde_HGNC <- dataMatrixde_HGNC[, !(colnames(dataMatrixde_HGNC) %in% c("HGNC.symbol"))]
dataMatrixde_HGNC <- as.matrix(data.frame(dataMatrixde_HGNC, check.names = FALSE))
class(dataMatrixde_HGNC) <- "numeric"
write.csv(dataMatrixde_HGNC, file.path(deresults_path, "full_expression_matrix_de_lm_hgnc.csv"), quote = F)

message(paste0("Dimensionality of DE genes ", dim(ExpressMatrixde)[1]))

# --Change column names and create matrix for coloring columns
new_colnames <- c()

for (i in colnames(ExpressMatrixde)) {
    i <- str_remove_all(i, "ti")
    i <- str_remove_all(i, "treat")
    i <- str_remove(i, " -BL.\\w+$")
    new_colnames <- c(new_colnames, i)
    # new_colnames <- c(new_colnames, v)

}
new_colnames_heatmapde <- new_colnames
colnames(ExpressMatrixde) <- new_colnames
colnames(results.lmde) <- new_colnames
colnames(dataMatrixde) <- new_colnames
colnames(dataMatrixde_HGNC) <- new_colnames
colnames(ExpressMatrixde_HGNC) <- new_colnames

colcolorlistgroup <- c(rep("pink", 5), rep("black", 5))
colcolormatrix <- as.matrix(colcolorlistgroup)
colnames(colcolormatrix) <- c("Group")

# --barplot
vizualize_DE_genes_bp(results.lmde, file.path(deresults_path, "barplot.png"))

# --heatmap
png(file.path(deresults_path, "heatmap.png"), width = 8, height = 10, units = "in", res = 300)
# par(mar = c(4, 4, -1, 2))
global_modulesde <- heatmap.L.4(ExpressMatrixde, figmargins=c(7,5),
    cutoff = 1, distmethod = "euclidean", cexcol = 2, labCol=new_colnames,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize=0.9,
    colsep = c(5), colcolorlist = colcolormatrix
)
dev.off()

####################################################
# make heatmap of dark red module, dark blue, black
#####################################################
##--darkred##
genes <- which(global_modulesde$modulesrows == "darkred")
df <- ExpressMatrixde[rownames(ExpressMatrixde) %in% names(genes), ]
df_trans <- merge(rhesus2human, df,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

write.csv(df_trans, file.path(deresults_path, "darkredmodulegenes.csv"))
dffinal <-  df_trans
genenames <- dffinal$HGNC.symbol
dffinal$HGNC.symbol <- NULL
# rownames(dffinal) <- dffinal$Gene.stable.ID
dffinal$Gene.stable.ID <- NULL
dffinal <- as.matrix(as.data.frame(dffinal))
png(file.path(deresults_path, "heatmap_darkredmodule.png"), width = 5, height = 7, units = "in", res = 300)
# par(mar = c(4, 4, -1, 2))
global_modulesdarkred <- heatmap.L.4(dffinal,
    figmargins = c(7, 6),
    cutoff = 1, distmethod = "euclidean", cexcol = 2, labCol = new_colnames, labRow=genenames, cexrow=1,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9,
    colsep = c(5), colcolorlist = colcolormatrix
)
dev.off()

##--darkblue--##
genes <- which(global_modulesde$modulesrows == "darkblue")
df <- ExpressMatrixde[rownames(ExpressMatrixde) %in% names(genes), ]
df_trans <- merge(rhesus2human, df,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

write.csv(df_trans, file.path(deresults_path, "darkbluemodulegenes.csv"))
dffinal <- df_trans
counter <- 1
for (i in dffinal$HGNC.symbol) {
    if (i == "") {
        print(i)
        print(dffinal$Gene.stable.ID[counter])
        dffinal$HGNC.symbol[counter] <- dffinal$Gene.stable.ID[counter]
    }
    counter <- counter + 1
}
genenames <- dffinal$HGNC.symbol
dffinal$HGNC.symbol <- NULL
# rownames(dffinal) <- dffinal$Gene.stable.ID
dffinal$Gene.stable.ID <- NULL
dffinal <- as.matrix(as.data.frame(dffinal))
png(file.path(deresults_path, "heatmap_darkbluemodule.png"), width = 5, height = 7, units = "in", res = 300)
# par(mar = c(4, 4, -1, 2))
global_modulesdarkblue <- heatmap.L.4(dffinal,
    figmargins = c(7, 11),
    cutoff = 1, distmethod = "euclidean", cexcol = 2, labCol = new_colnames, labRow = genenames, cexrow = 1,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9,
    colsep = c(5), colcolorlist = colcolormatrix
)
dev.off()

##--black--##

genes <- which(global_modulesde$modulesrows == "black")
df <- ExpressMatrixde[rownames(ExpressMatrixde) %in% names(genes), ]
df_trans <- merge(rhesus2human, df,
    by.x = "Gene.stable.ID",
    by.y = "row.names",
    all.X = T, all.Y = T
)

write.csv(df_trans, file.path(deresults_path, "blackmodulegenes.csv"))
dffinal <- df_trans
counter <- 1
for( i in dffinal$HGNC.symbol) {
    if(i=="") {
        print (i)
        print(dffinal$Gene.stable.ID[counter])
        dffinal$HGNC.symbol[counter] <- dffinal$Gene.stable.ID[counter]
    }
    counter = counter+1
}
genenames <- dffinal$HGNC.symbol
dffinal$HGNC.symbol <- NULL
# rownames(dffinal) <- dffinal$Gene.stable.ID
dffinal$Gene.stable.ID <- NULL
dffinal <- as.matrix(as.data.frame(dffinal))
png(file.path(deresults_path, "heatmap_blackmodule.png"), width = 5, height = 7, units = "in", res = 300)
# par(mar = c(4, 4, -1, 2))
global_modulesblack <- heatmap.L.4(dffinal,
    figmargins = c(7, 8),
    cutoff = 1, distmethod = "euclidean", cexcol = 2, labCol = new_colnames, labRow = genenames, cexrow = 0.7,
    clustermethod = "ward.D2", clusterdim = "row", ColSideColorsSize = 0.9,
    colsep = c(5), colcolorlist = colcolormatrix
)
dev.off()
###############################
# Enrichment Analysis on DE
###############################

# ---Build collection based on reference genes extracted (only have to build this once)
if (isTRUE(SetRankRun)) {

    # allGenesTrans <- rhesus2human[rhesus2human$Gene.stable.ID %in% rownames(Pi.CPM$E), ]
    # allGenesHGNC <- unique(unlist(allGenesTrans$HGNC.symbol))
    # referenceSet <- symbol2EntrezID(allGenesHGNC)

    # collection <- buildSetCollection(allDBs,
    #     referenceSet = referenceSet,
    #     maxSetSize = 500
    # )
    # saveRDS(collection, "./collectionallDBs.rds")
    collection <- readRDS("./collectionallDBs.rds")
}

if (isTRUE(SetRankRun)) {
    message("STATUS: Finding gene enrichments for each cluster in DE")
    for (cluster in unique(global_modulesde$modulesrows)) {
        print(paste0("STATUS: gene enrichments for module ", cluster))
        genes <- which(global_modules$modulesrows == cluster)
        network <- gene_enrichment(names(genes), file.path(deresults_path, "SetRank_results/"), cluster)
    }
}
if (isTRUE(SetRankRun)) {
    network <- gene_enrichment(rownames(ExpressMatrixde), file.path(deresults_path, "SetRank_results/"), "All")
}

######################
# Line Plots
#####################
Lineplots_results <- "1.DELinePlots"
generate_folder(file.path(deresults_path, Lineplots_results))

data_summary <- function(data, varname, groupnames) {
    require(plyr)
    summary_func <- function(x, col) {
        c(
            sum = sum(x[[col]], na.rm = TRUE),
            sd = sd(x[[col]], na.rm = TRUE)
        )
    }
    data_sum <- ddply(data, groupnames,
        .fun = summary_func,
        varname
    )
    return(data_sum)
}

lineplotsdraw <- function(path2file, membershipfile, pathofinterest, resultsfig) {
    df <- read.table(path2file, sep = "\t")
    dfgenes <- read.table(membershipfile, row.names = 1, header = TRUE)
    genes <- rownames(dfgenes[dfgenes[, pathofinterest] == "X", ])
    print(genes)
    write.csv(genes, str_replace(resultsfig, ".png", ".csv"))

    genesensembl <- rhesus2human[rhesus2human$HGNC.symbol %in% genes, "Gene.stable.ID"]
    expression <- Pi.CPM$E[rownames(Pi.CPM$E) %in% genesensembl, colnames(Pi.CPM$E) %in% rownames(target)]
    if (length(genes) >1 ) {
        colnames(expression) <- paste(target$TimePoint, target$Sample_Type, sep = ".")
        mexp <- melt(expression)
        time <- c()
        treatment <- c()
        for (x in mexp$Var2) {
            v <- strsplit(x, "\\.")
            time <- c(time, v[[1]][1])
            treatment <- c(treatment, v[[1]][2])
        }
        mexpdraw <- mexp
        mexpdraw$time <- time
        mexpdraw$treatment <- treatment

        df <- aggregate(mexp[, c("value")], list(mexp$Var2), median)
        dfsd <- aggregate(mexp[, c("value")], list(mexp$Var2), sd)
        df$sd <- dfsd$x
        print(head(df))
        dfgenemean <- aggregate(mexp[, c("value")], list(mexp$Var1,mexp$Var2 ), mean)

        # -- df processing #
        time <- c()
        treatment <- c()

        for (x in dfgenemean$Group.2) {
            v <- strsplit(x, "\\.")
            time <- c(time, v[[1]][1])
            treatment <- c(treatment, v[[1]][2])
        }
        dfgenemean$time <- time
        dfgenemean$treatment <- treatment
        dfgenemean$time <- factor(dfgenemean$time,
            levels = c(
                "BL", "D1", "2", "4", "5", "7"
            )
        )

        # -- dfmedian processing #
        time <- c()
        treatment <- c()

        for (x in df$Group.1) {
            v <- strsplit(x, "\\.")
            time <- c(time, v[[1]][1])
            treatment <- c(treatment, v[[1]][2])
        }
        df$time <- time
        df$treatment <- treatment
        df$time <- factor(df$time,
            levels = c(
                "BL", "D1", "2", "4", "5", "7"
            )
        )

        ggplot(data = dfgenemean, aes(x = time, y = x,  color=treatment)) +
            # geom_errorbar(aes(ymin = sum - sd, ymax = x + sd)) +
            scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            geom_point() +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        ggsave(str_replace(resultsfig, ".png", "_all.png"), width = 6, height = 3)


        ggplot(data = df, aes(x = time, y = x, group = treatment, color = treatment)) +
            geom_line() + labs(x="time", y="median", color="") +
            # geom_errorbar(aes(ymin = x - sd, ymax = x + sd)) +
            scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            geom_point() +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        ggsave(resultsfig, width = 4.5, height = 3)
    } else {
        message("Only one gene in path") 
        names(expression) <- paste(target$TimePoint, target$Sample_Type, sep = ".")
        mexp <- melt(expression)
        mexp$Var2 <- paste(target$TimePoint, target$Sample_Type, sep = ".")
        time <- c()
        treatment <- c()
        for (x in mexp$Var2) {
            v <- strsplit(x, "\\.")
            time <- c(time, v[[1]][1])
            treatment <- c(treatment, v[[1]][2])
        }
        mexpdraw <- mexp
        mexpdraw$time <- time
        mexpdraw$treatment <- treatment
        mexpdraw$time <- factor(mexpdraw$time,
            levels = c(
                "BL", "D1", "2", "4", "5", "7"
            )
        )

        df <- aggregate(mexpdraw[, c("value")], list(mexpdraw$Var2), median)
        dfsd <- aggregate(mexpdraw[, c("value")], list(mexpdraw$Var2), sd)
        df$sd <- dfsd$x
        time <- c()
        treatment <- c()

        for (x in df$Group.1) {
            v <- strsplit(x, "\\.")
            time <- c(time, v[[1]][1])
            treatment <- c(treatment, v[[1]][2])
        }
        df$time <- time
        df$treatment <- treatment
        df$time <- factor(df$time,
            levels = c(
                "BL", "D1", "2", "4", "5", "7"
            )
        )
        ggplot(data = df, aes(x = time, y = x, group = treatment, color = treatment)) +
            geom_line() +
            labs(x = "time", y = "median", color = "") +
            # geom_errorbar(aes(ymin = x - sd, ymax = x + sd)) +
            scale_color_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            scale_fill_manual(values = c("IFN" = "pink", "CTL" = "black")) +
            geom_point() +
            theme_Publication() +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        ggsave(resultsfig, width = 4.5, height = 3)
    }
}

lineplotsdraw( file.path(deresults_path, "SetRank_results/", "de_unranked_black_pathways.txt"),
               file.path(deresults_path, "SetRank_results/", "de_unranked_black_membership.txt"),
               "GO.0060337", file.path(deresults_path, Lineplots_results, "black_typeIinterferonsignaling.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_black_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_black_membership.txt"),
    "GO.0045087", file.path(deresults_path, Lineplots_results, "black_InnateImmuneresponse.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_skyblue_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_skyblue_membership.txt"),
    "GO.0002526", file.path(deresults_path, Lineplots_results, "skyblue_accuteinflammaatory.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_skyblue_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_skyblue_membership.txt"),
    "GO.0061515", file.path(deresults_path, Lineplots_results, "skyblue_myeloidcelldevelopment.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkred_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkred_membership.txt"),
    "GO.0034344", file.path(deresults_path, Lineplots_results, "darkred_reg_type_III_ifn.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkred_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkred_membership.txt"),
    "GO.0034155", file.path(deresults_path, Lineplots_results, "darkred_TLR7.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkblue_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_darkblue_membership.txt"),
    "GO.0045087", file.path(deresults_path, Lineplots_results, "darkblue_InnateImmuneresponse.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_yellow_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_yellow_membership.txt"),
    "GO.0001817", file.path(deresults_path, Lineplots_results, "yellow_regulationofcytokineproduction.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_yellow_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_yellow_membership.txt"),
    "hsa04668", file.path(deresults_path, Lineplots_results, "yellow_TNFsignaling.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_membership.txt"),
    "GO.0023024", file.path(deresults_path, Lineplots_results, "red_MHCcomplex.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_membership.txt"),
    "GO.0004888", file.path(deresults_path, Lineplots_results, "red_transmembrane.png")
)
lineplotsdraw(
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_pathways.txt"),
    file.path(deresults_path, "SetRank_results/", "de_unranked_red_membership.txt"),
    "GO.0032817", file.path(deresults_path, Lineplots_results, "red_regulationofNKproliferation.png")
)
