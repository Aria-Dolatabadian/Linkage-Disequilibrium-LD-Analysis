# =============================================================================
# Linkage Disequilibrium (LD) Analysis from XLSX Genotype Data
# =============================================================================
# Required packages: readxl, corrplot, ggplot2, reshape2
#
# Install if needed:
#   install.packages(c("readxl", "corrplot", "ggplot2", "reshape2"))
# =============================================================================

library(readxl)
library(corrplot)
library(ggplot2)
library(reshape2)

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA FROM XLSX
# ─────────────────────────────────────────────────────────────────────────────

xlsx_file <- "snp_ld_data.xlsx"   # <── set path to your file

snp_info  <- read_excel(xlsx_file, sheet = "SNP_Info")
geno_raw  <- read_excel(xlsx_file, sheet = "Genotypes")
snp_stats <- read_excel(xlsx_file, sheet = "Summary_Stats")

cat("=== Data Loaded ===\n")
cat(sprintf("Samples : %d\n", nrow(geno_raw)))
cat(sprintf("SNPs    : %d\n", ncol(geno_raw) - 1))   # minus Sample_ID column

# ─────────────────────────────────────────────────────────────────────────────
# 2. PREPARE GENOTYPE MATRIX
# ─────────────────────────────────────────────────────────────────────────────
# Dosage coding: 0 = hom-ref, 1 = het, 2 = hom-alt

sample_ids <- geno_raw$Sample_ID
snp_ids    <- colnames(geno_raw)[-1]          # drop Sample_ID column
geno_mat   <- as.matrix(geno_raw[, -1])       # numeric matrix: rows=samples, cols=SNPs
rownames(geno_mat) <- sample_ids
colnames(geno_mat) <- snp_ids

cat("\nGenotype matrix (first 5 samples, first 5 SNPs):\n")
print(geno_mat[1:5, 1:5])

# ─────────────────────────────────────────────────────────────────────────────
# 3. QUALITY CONTROL
# ─────────────────────────────────────────────────────────────────────────────

# Minor allele frequency per SNP
compute_maf <- function(g) {
  freq <- colMeans(g, na.rm = TRUE) / 2
  pmin(freq, 1 - freq)
}

# Missing rate per SNP
missing_rate <- colMeans(is.na(geno_mat))

# Hardy-Weinberg Equilibrium exact test (chi-square approximation)
hwe_test <- function(g_col) {
  n  <- sum(!is.na(g_col))
  n0 <- sum(g_col == 0, na.rm = TRUE)   # hom-ref
  n1 <- sum(g_col == 1, na.rm = TRUE)   # het
  n2 <- sum(g_col == 2, na.rm = TRUE)   # hom-alt
  q  <- (2 * n2 + n1) / (2 * n)
  p  <- 1 - q
  e0 <- p^2 * n; e1 <- 2*p*q*n; e2 <- q^2 * n
  chi2 <- sum(c((n0-e0)^2/max(e0,1e-9),
                (n1-e1)^2/max(e1,1e-9),
                (n2-e2)^2/max(e2,1e-9)))
  pchisq(chi2, df = 1, lower.tail = FALSE)
}

maf_vals  <- compute_maf(geno_mat)
hwe_pvals <- apply(geno_mat, 2, hwe_test)

qc_table <- data.frame(
  SNP          = snp_ids,
  MAF          = round(maf_vals, 4),
  Missing_Rate = round(missing_rate, 4),
  HWE_p        = round(hwe_pvals, 4),
  Pass_QC      = maf_vals >= 0.05 & missing_rate <= 0.05 & hwe_pvals >= 0.001
)

cat("\n=== Quality Control Summary ===\n")
print(qc_table)

# Keep only SNPs that pass QC
pass_snps <- qc_table$SNP[qc_table$Pass_QC]
geno_qc   <- geno_mat[, pass_snps, drop = FALSE]
cat(sprintf("\n%d / %d SNPs passed QC\n", length(pass_snps), length(snp_ids)))

# ─────────────────────────────────────────────────────────────────────────────
# 4. COMPUTE LD STATISTICS (r² and D')
# ─────────────────────────────────────────────────────────────────────────────

# Manual computation of r² and D' for a pair of biallelic SNPs (dosage coding)
compute_ld_pair <- function(g1, g2) {
  keep  <- !is.na(g1) & !is.na(g2)
  g1    <- g1[keep]; g2 <- g2[keep]
  n     <- length(g1)
  
  # Allele frequencies
  p1 <- mean(g1) / 2   # freq of alt allele at SNP1
  p2 <- mean(g2) / 2   # freq of alt allele at SNP2
  q1 <- 1 - p1;  q2 <- 1 - p2
  
  # Haplotype frequency estimation (EM-free shortcut under HWE)
  # p11 ≈ E[g1*g2] / 4  (prob both alt haplotypes together)
  p11 <- mean(g1 * g2) / 4
  
  D   <- p11 - p1 * p2
  
  # D' (Lewontin's)
  if (D >= 0) Dmax <- min(p1*q2, q1*p2) else Dmax <- min(p1*p2, q1*q2)
  Dprime <- if (Dmax > 0) D / Dmax else NA_real_
  
  # r²
  r2 <- D^2 / (p1 * q1 * p2 * q2)
  r2 <- if (is.nan(r2) || is.infinite(r2)) NA_real_ else r2
  
  # Correlation-based r (can be negative)
  r  <- cor(g1, g2, use = "complete.obs")
  
  list(r2 = r2, Dprime = Dprime, r = r, D = D)
}

n_snp <- ncol(geno_qc)
snp_names <- colnames(geno_qc)
r2_mat     <- matrix(NA, n_snp, n_snp, dimnames = list(snp_names, snp_names))
dprime_mat <- matrix(NA, n_snp, n_snp, dimnames = list(snp_names, snp_names))
diag(r2_mat) <- 1; diag(dprime_mat) <- 1

for (i in 1:(n_snp - 1)) {
  for (j in (i + 1):n_snp) {
    ld <- compute_ld_pair(geno_qc[, i], geno_qc[, j])
    r2_mat[i, j] <- r2_mat[j, i] <- ld$r2
    dprime_mat[i, j] <- dprime_mat[j, i] <- abs(ld$Dprime)
  }
}

cat("\n=== Pairwise r² Matrix ===\n")
print(round(r2_mat, 3))

cat("\n=== Pairwise |D'| Matrix ===\n")
print(round(dprime_mat, 3))

# ─────────────────────────────────────────────────────────────────────────────
# 5. LD BLOCKS (simple r² threshold method)
# ─────────────────────────────────────────────────────────────────────────────

ld_threshold <- 0.5

find_ld_blocks <- function(r2_mat, snp_df, threshold = 0.5) {
  n      <- nrow(r2_mat)
  blocks <- list()
  assigned <- rep(FALSE, n)
  
  for (i in 1:n) {
    if (assigned[i]) next
    block_snps <- i
    for (j in (i+1):n) {
      if (j > n) break
      if (snp_df$Chrom[j] != snp_df$Chrom[i]) break  # different chromosome
      if (!is.na(r2_mat[i, j]) && r2_mat[i, j] >= threshold) {
        block_snps <- c(block_snps, j)
        assigned[j] <- TRUE
      }
    }
    assigned[i] <- TRUE
    blocks[[length(blocks) + 1]] <- rownames(r2_mat)[block_snps]
  }
  blocks
}

snp_info_qc <- snp_info[snp_info$SNP_ID %in% snp_names, ]
ld_blocks   <- find_ld_blocks(r2_mat, snp_info_qc, threshold = ld_threshold)

cat(sprintf("\n=== LD Blocks (r² ≥ %.1f) ===\n", ld_threshold))
for (b in seq_along(ld_blocks)) {
  cat(sprintf("Block %d: %s\n", b, paste(ld_blocks[[b]], collapse = ", ")))
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. VISUALISATIONS
# ─────────────────────────────────────────────────────────────────────────────

# ── 6a. r² heatmap ────────────────────────────────────────────────────────

r2_melt <- melt(r2_mat, na.rm = TRUE)
colnames(r2_melt) <- c("SNP1", "SNP2", "r2")

p_r2 <- ggplot(r2_melt, aes(x = SNP1, y = SNP2, fill = r2)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r2)), size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "#FFFFFF",
    mid      = "#FFA07A",
    high     = "#8B0000",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "r²"
  ) +
  scale_x_discrete(limits = snp_names) +
  scale_y_discrete(limits = rev(snp_names)) +
  labs(title = "Pairwise Linkage Disequilibrium (r²)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

ggsave("ld_r2_heatmap.png", p_r2, width = 9, height = 7, dpi = 150)
cat("\nSaved: ld_r2_heatmap.png\n")

# ── 6b. D' heatmap ────────────────────────────────────────────────────────

dp_melt <- melt(dprime_mat, na.rm = TRUE)
colnames(dp_melt) <- c("SNP1", "SNP2", "Dprime")

p_dp <- ggplot(dp_melt, aes(x = SNP1, y = SNP2, fill = Dprime)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Dprime)), size = 3, color = "black") +
  scale_fill_gradient2(
    low      = "#FFFFFF",
    mid      = "#87CEEB",
    high     = "#00008B",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "|D'|"
  ) +
  scale_x_discrete(limits = snp_names) +
  scale_y_discrete(limits = rev(snp_names)) +
  labs(title = "Pairwise Linkage Disequilibrium (|D'|)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

ggsave("ld_dprime_heatmap.png", p_dp, width = 9, height = 7, dpi = 150)
cat("Saved: ld_dprime_heatmap.png\n")

# ── 6b-alt. corrplot LD visualization ─────────────────────────────────────
# corrplot gives a compact, publication-style triangular LD matrix

png("ld_corrplot.png", width = 800, height = 800, res = 120)
corrplot(
  r2_mat,
  method      = "color",
  type        = "upper",
  col         = colorRampPalette(c("white", "#FFA07A", "#8B0000"))(100),
  addCoef.col = "black",
  number.cex  = 0.75,
  tl.col      = "black",
  tl.srt      = 45,
  tl.cex      = 0.85,
  cl.lim      = c(0, 1),
  title       = "LD Matrix (r²) — corrplot",
  mar         = c(0, 0, 2, 0)
)
dev.off()
cat("Saved: ld_corrplot.png\n")

# ── 6c. MAF bar chart ─────────────────────────────────────────────────────

maf_df <- data.frame(SNP = snp_ids, MAF = maf_vals,
                     Pass = ifelse(snp_ids %in% pass_snps, "Pass", "Fail (MAF < 0.05)"))

p_maf <- ggplot(maf_df, aes(x = reorder(SNP, -MAF), y = MAF, fill = Pass)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 1, y = 0.07, label = "MAF threshold = 0.05",
           hjust = 0, size = 3.5, color = "red") +
  scale_fill_manual(values = c("Fail (MAF < 0.05)" = "#E74C3C", "Pass" = "#2ECC71"),
                    name = "QC Status") +
  labs(title = "Minor Allele Frequencies per SNP",
       x = "SNP", y = "MAF") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold", hjust = 0.5))

ggsave("maf_barplot.png", p_maf, width = 8, height = 5, dpi = 150)
cat("Saved: maf_barplot.png\n")

# ── 6d. LD decay plot ─────────────────────────────────────────────────────
# r² vs. physical distance for SNPs on the same chromosome

ld_decay <- data.frame()
for (i in 1:(n_snp - 1)) {
  for (j in (i + 1):n_snp) {
    sn_i <- snp_names[i]; sn_j <- snp_names[j]
    ri   <- snp_info[snp_info$SNP_ID == sn_i, ]
    rj   <- snp_info[snp_info$SNP_ID == sn_j, ]
    if (nrow(ri) == 0 || nrow(rj) == 0) next
    if (ri$Chrom != rj$Chrom) next
    dist_bp <- abs(rj$Position - ri$Position)
    ld_decay <- rbind(ld_decay, data.frame(
      dist_kb = dist_bp / 1000,
      r2      = r2_mat[i, j],
      Chrom   = paste0("Chr", ri$Chrom)
    ))
  }
}

if (nrow(ld_decay) > 0) {
  p_decay <- ggplot(ld_decay, aes(x = dist_kb, y = r2, color = Chrom)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_smooth(method = "loess", se = TRUE, span = 1.0, alpha = 0.15) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "LD Decay: r² vs. Physical Distance",
         x = "Distance (kb)", y = "r²", color = "Chromosome") +
    ylim(0, 1) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave("ld_decay.png", p_decay, width = 8, height = 5, dpi = 150)
  cat("Saved: ld_decay.png\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# 7. SUMMARY REPORT
# ─────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 55), "\n")
cat("            LINKAGE DISEQUILIBRIUM ANALYSIS REPORT\n")
cat(strrep("=", 55), "\n\n")

cat(sprintf("Dataset         : %d samples, %d SNPs\n", nrow(geno_mat), ncol(geno_mat)))
cat(sprintf("SNPs passing QC : %d\n\n", ncol(geno_qc)))

cat("Top 10 SNP pairs by r²:\n")
r2_pairs <- subset(r2_melt, as.character(SNP1) < as.character(SNP2))
r2_pairs  <- r2_pairs[order(-r2_pairs$r2), ]
print(head(r2_pairs, 10), row.names = FALSE)

cat("\nTop 10 SNP pairs by |D'|:\n")
dp_pairs <- subset(dp_melt, as.character(SNP1) < as.character(SNP2))
dp_pairs  <- dp_pairs[order(-dp_pairs$Dprime), ]
print(head(dp_pairs, 10), row.names = FALSE)

cat("\nOutputs:\n")
cat("  • ld_r2_heatmap.png     — pairwise r² heatmap (ggplot2)\n")
cat("  • ld_dprime_heatmap.png — pairwise |D'| heatmap (ggplot2)\n")
cat("  • ld_corrplot.png       — upper-triangle r² matrix (corrplot)\n")
cat("  • maf_barplot.png       — MAF bar chart with QC threshold\n")
cat("  • ld_decay.png          — LD decay by physical distance\n")
cat(strrep("=", 55), "\n")
