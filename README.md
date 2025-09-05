# SNP Density Ideogram Plotting with Centromeres (R & RIdeogram)

This document explains the workflow for generating SNP density ideogram plots with centromere annotation, using VCF data processed via VCFtools and plotted in R with the [RIdeogram](https://cran.r-project.org/package=RIdeogram) package.

---

## **1. Preparation of SNP Density Data**

Before running the R script, you need to generate the SNP density data files. Here are the steps:

### **a. Compress and Index the VCF Files**

If your VCF files are not already compressed and indexed, do:

```bash
# Compress and index the full dataset VCF
bgzip -c gwas_ready_filtered.vcf > gwas_ready_filtered.vcf.gz   # (skip if already gzipped)
tabix -p vcf gwas_ready_filtered.vcf.gz

# Compress and index the subset VCF
bgzip -c gwas_ready_filtered_20.vcf > gwas_ready_filtered_20.vcf.gz  # (skip if already gzipped)
tabix -p vcf gwas_ready_filtered_20.vcf.gz
```

### **b. Generate SNP Density Files with VCFtools**

Use [VCFtools](https://vcftools.github.io/) to calculate SNP density in 1Mb non-overlapping windows:

```bash
# For the first file (all samples)
vcftools --gzvcf gwas_ready_filtered.vcf.gz --SNPdensity 1000000 --out snp_density_all

# For the second file (subset of 20 samples)
vcftools --gzvcf gwas_ready_filtered_20.vcf.gz --SNPdensity 1000000 --out snp_density_20
```

This will produce `.snpden` files as input for the R script.

---

## **2. Plotting in R: SNP Density Ideograms**

### **a. Script Overview**

The script (`snp_density_RIdeogram_with_centomere.R`) reads SNP density files and a centromere annotation file, then plots ideograms, overlaying SNP density and centromere positions.

### **b. Environment**

This step was performed in a [conda](https://docs.conda.io/) environment called `r_env`.

**Environment setup example:**

```bash
conda create -n r_env r-base=4.2 data.table
conda activate r_env
# Install RIdeogram in R:
# install.packages("RIdeogram", repos="https://cloud.r-project.org")
```

### **c. Directory**

All analysis was performed in:

```
/home/mahmoudi/Bigdata/computing/Shima/b_carinata/snp_density_plot
```

---

## **3. Running the R Plotting Script**

Set your working directory in R or use full file paths. Execute the script as follows (after activating `r_env`):

```bash
cd /home/mahmoudi/Bigdata/computing/Shima/b_carinata/snp_density_plot
Rscript snp_density_RIdeogram_with_centomere.R
```

---

## **4. Script Details**

- **Inputs:**
  - `snp_density_all.snpden`: SNP density file for all samples.
  - `snp_density_20.snpden`: SNP density file for 20 samples.
  - `centromeres.tsv`: Tab-delimited file with columns `chr`, `arm`, `start`, `end` (centromere annotation).
- **Output:** 
  - SVG ideograms visualizing SNP density and centromeres (e.g., `snp_density_all.svg`, `snp_density_20.svg`).

---

## **5. R Script Reference**

```r
# --- SNP density ideograms with centromeres (RIdeogram, two separate plots) ---

options(stringsAsFactors = FALSE)

# Packages
need <- c("data.table", "RIdeogram")
for (p in need) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
library(data.table)
library(RIdeogram)

# --------- SET YOUR FILES HERE ---------
fileA <- "snp_density_all.snpden"
fileB <- "snp_density_20.snpden"
centromere_file <- "centromeres.tsv"   # your table (see expected columns below)
# --------------------------------------

# Your chromosome mapping (CM... -> B/C names)
chr_map <- c(
  "CM081008.1"="B1","CM081009.1"="B2","CM081010.1"="B3",
  "CM081011.1"="B4","CM081012.1"="B5","CM081013.1"="B6",
  "CM081014.1"="B7","CM081015.1"="B8","CM081016.1"="C1",
  "CM081017.1"="C2","CM081018.1"="C3","CM081019.1"="C4",
  "CM081020.1"="C5","CM081021.1"="C6","CM081022.1"="C7",
  "CM081023.1"="C8","CM081024.1"="C9"
)
desired_order <- c(paste0("B",1:8), paste0("C",1:9))

# ---------- Helpers ----------
read_snpden_upper <- function(path){
  stopifnot(file.exists(path))
  dt <- fread(path, showProgress = FALSE)
  if (nrow(dt) == 0) stop("File has 0 rows: ", path)
  
  orig <- names(dt)
  setnames(dt, names(dt), toupper(names(dt)))
  
  # Count column → Value
  cnt_candidates <- c("N_SNPS","N_SNP","SNP_COUNT")
  cnt_col <- intersect(cnt_candidates, names(dt))
  if (length(cnt_col) == 0)
    stop("No SNP count column in ", path, ". Columns: ", paste(orig, collapse=", "))
  setnames(dt, cnt_col[1], "COUNT")
  
  # Window coords → Start/End
  if ("BIN_START" %in% names(dt)) {
    START <- dt$BIN_START
    if ("BIN_END" %in% names(dt)) {
      END <- dt$BIN_END
    } else {
      bw  <- unique(na.omit(diff(sort(unique(dt$BIN_START)))))
      win <- if (length(bw) >= 1) bw[1] else 1e6
      END <- START + win - 1
    }
  } else if ("BIN" %in% names(dt)) {
    START <- dt$BIN
    bw  <- unique(na.omit(diff(sort(unique(dt$BIN)))))
    win <- if (length(bw) >= 1) bw[1] else 1e6
    END <- START + win - 1
  } else if ("POS" %in% names(dt)) {
    START <- dt$POS; END <- dt$POS
  } else {
    stop("No BIN_START/BIN_END/BIN/POS in ", path, ". Columns: ", paste(orig, collapse=", "))
  }
  if (!"CHROM" %in% names(dt)) stop("No CHROM column in ", path)
  
  out <- data.table(
    Chr_old = as.character(dt$CHROM),
    Start   = as.integer(START),
    End     = as.integer(END),
    Value   = suppressWarnings(as.numeric(gsub("[^0-9eE+\\.-]", "", as.character(dt$COUNT))))
  )
  out[!is.finite(Value) | is.na(Value), Value := 0]
  out[, Chr := chr_map[Chr_old]]
  out <- out[!is.na(Chr)]
  if (nrow(out) == 0) stop("After mapping, no rows remain for ", path, ". Check chr_map keys vs CHROM values.")
  out[, `:=`(Start = pmax(1L, Start), End = pmax(Start, End))]
  setcolorder(out, c("Chr","Start","End","Value"))
  out[, Chr_old := NULL]
  setorder(out, Chr, Start, End)
  as.data.frame(out)
}

# centromere_file expected columns (tab-delimited):
# chr  arm  start  end  color_group
# (like the table you pasted; chr values are B1..C9)
read_centromeres <- function(path){
  stopifnot(file.exists(path))
  cdf <- fread(path)
  setnames(cdf, names(cdf), tolower(names(cdf)))
  must <- c("chr","arm","start","end")
  if (!all(must %in% names(cdf)))
    stop("centromere_file must have columns: chr, arm, start, end (and optional color_group)")
  cdf[, `:=`(chr = as.character(chr),
             arm = tolower(as.character(arm)),
             start = as.integer(start),
             end   = as.integer(end))]
  # Centromere breakpoint = start of q arm
  ce <- cdf[arm == "q", .(Chr = chr, CE_break = start)]
  # Chromosome length from arm 'end'
  lens <- cdf[, .(End_arm = max(end, na.rm = TRUE)), by = chr]
  setnames(lens, "chr", "Chr")
  ce <- merge(ce, lens, by = "Chr", all.x = TRUE)
  # Define a visible centromere region around the breakpoint (~1% of length, clamped)
  ce[, width := pmax(1e5, pmin(2e6, as.numeric(round(0.01 * End_arm))))]
  ce[, `:=`(CE_start = pmax(1L, as.integer(CE_break - width/2)),
            CE_end   = pmin(as.integer(End_arm), as.integer(CE_break + width/2)))]
  ce[, .(Chr, CE_start, CE_end, ChromLen = as.integer(End_arm))]
}

# Build 5-column karyotype: Chr, Start, End, CE_start, CE_end
build_karyotype_with_centromere <- function(over_upper, centromeres){
  # infer End from SNP bins too, just in case
  ends_from_over <- as.data.table(over_upper)[, .(End_over = max(End, na.rm = TRUE)), by = Chr]
  k <- merge(centromeres, ends_from_over, by = "Chr", all = TRUE)
  k[, End := pmax(ChromLen, End_over, na.rm = TRUE)]
  k[!is.finite(End) | is.na(End), End := ChromLen]
  k[, Start := 1L]
  k <- k[, .(Chr, Start, End, CE_start, CE_end)]
  # keep/order requested chromosomes
  ord <- desired_order[desired_order %in% k$Chr]
  if (length(ord) == 0) ord <- sort(unique(k$Chr))
  k <- k[match(ord, k$Chr), , drop = FALSE]
  # safety clamps
  k$CE_start <- pmax(1L, pmin(k$CE_start, k$End-1L))
  k$CE_end   <- pmax(k$CE_start+1L, pmin(k$CE_end,   k$End))
  as.data.frame(k)
}

plot_one_ideogram <- function(karyotype5, overlaid4, out_prefix,
                              palette = c("#f7fbff","#6baed6","#08306b")) {
  # EXACT columns/case for older RIdeogram
  karyo_use <- data.frame(
    Chr      = as.character(karyotype5$Chr),
    Start    = as.integer(karyotype5$Start),
    End      = as.integer(karyotype5$End),
    CE_start = as.integer(karyotype5$CE_start),
    CE_end   = as.integer(karyotype5$CE_end),
    stringsAsFactors = FALSE
  )
  over_use <- data.frame(
    Chr   = as.character(overlaid4$Chr),
    Start = as.integer(overlaid4$Start),
    End   = as.integer(overlaid4$End),
    Value = as.numeric(overlaid4$Value),
    stringsAsFactors = FALSE
  )
  over_use$Value[!is.finite(over_use$Value) | is.na(over_use$Value)] <- 0
  # clamp bins to chrom end
  end_map <- setNames(karyo_use$End, karyo_use$Chr)
  over_use$End <- pmin(over_use$End, end_map[over_use$Chr])
  
  cat("\n[Diag] Karyotype (5 cols) head:\n"); print(utils::head(karyo_use))
  cat("\n[Diag] Overlaid head (", out_prefix, "):\n", sep=""); print(utils::head(over_use))
  cat("\n[Diag] Value summary:\n"); print(summary(over_use$Value)); cat("\n")
  
  svg_file <- paste0(out_prefix, ".svg")
  ideogram(
    karyotype = karyo_use,   # 5 columns → centromere drawn
    overlaid  = over_use,
    colorset1 = palette,
    output    = svg_file
  )
  if (!file.exists(svg_file)) {  # very old versions write a default name
    for (cand in c("chromosome.svg","ideogram.svg","RIdeogram.svg"))
      if (file.exists(cand)) { file.rename(cand, svg_file); break }
  }
  message("Saved SVG: ", svg_file)
}

# ---------- Run: one dataset at a time ----------
stopifnot(file.exists(centromere_file))
centro <- read_centromeres(centromere_file)

make_for_file <- function(infile){
  over <- read_snpden_upper(infile)
  # Build karyotype (with centromeres) from THIS dataset + centromere table
  karyo <- build_karyotype_with_centromere(over, centro)
  prefix <- tools::file_path_sans_ext(basename(infile))
  plot_one_ideogram(karyo, over, out_prefix = prefix)
}

make_for_file(fileA)
make_for_file(fileB)
```

---

**For further reference or troubleshooting, contact: [shimahm](https://github.com/shimahm)**
