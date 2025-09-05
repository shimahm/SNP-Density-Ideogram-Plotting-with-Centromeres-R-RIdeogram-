# --- SNP density ideograms with centromeres + custom bins (0 grey; blue→yellow→red)
#     gap-filling fixed to preserve your original windows ---

options(stringsAsFactors = FALSE)

# Packages
need <- c("data.table", "RIdeogram")
for (p in need) if (!requireNamespace(p, quietly = TRUE))
  install.packages(p, repos = "https://cloud.r-project.org")
library(data.table)
library(RIdeogram)

# --------- SET YOUR FILES ----------
fileA <- "snp_density_all.snpden"
fileB <- "snp_density_20.snpden"
centromere_file <- "centromeres.txt"   # your centromere table
# -----------------------------------

# Chromosome mapping (CM... -> B/C)
chr_map <- c(
  "CM081008.1"="B1","CM081009.1"="B2","CM081010.1"="B3",
  "CM081011.1"="B4","CM081012.1"="B5","CM081013.1"="B6",
  "CM081014.1"="B7","CM081015.1"="B8","CM081016.1"="C1",
  "CM081017.1"="C2","CM081018.1"="C3","CM081019.1"="C4",
  "CM081020.1"="C5","CM081021.1"="C6","CM081022.1"="C7",
  "CM081023.1"="C8","CM081024.1"="C9"
)
desired_order <- c(paste0("B",1:8), paste0("C",1:9))

# ----- Bins & colors (0, then groups of 10; blue→yellow→red) -----
bin_breaks <- c(-Inf, 0, 9, 19, 29, 39, 49, Inf)
bin_labels <- c("0","1–9","10–19","20–29","30–39","40–49","≥50")
bin_colors <- c(
  "#BDBDBD",  # 0
  colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(length(bin_labels)-1)
)

# ---------- Helpers ----------
read_snpden_upper <- function(path){
  stopifnot(file.exists(path))
  dt <- fread(path, showProgress = FALSE)
  if (nrow(dt) == 0) stop("File has 0 rows: ", path)

  orig <- names(dt)
  setnames(dt, names(dt), toupper(names(dt)))

  # Count → Value
  cnt_candidates <- c("N_SNPS","N_SNP","SNP_COUNT")
  cnt_col <- intersect(cnt_candidates, names(dt))
  if (length(cnt_col) == 0)
    stop("No SNP count column in ", path, ". Columns: ", paste(orig, collapse=", "))
  setnames(dt, cnt_col[1], "COUNT")

  # Coordinates → Start/End
  if ("BIN_START" %in% names(dt)) {
    START <- dt$BIN_START
    END   <- if ("BIN_END" %in% names(dt)) dt$BIN_END else {
      bw <- unique(na.omit(diff(sort(unique(dt$BIN_START)))))
      win <- if (length(bw) >= 1) bw[1] else 1e6
      START + win - 1
    }
  } else if ("BIN" %in% names(dt)) {
    START <- dt$BIN
    bw <- unique(na.omit(diff(sort(unique(dt$BIN)))))
    win <- if (length(bw) >= 1) bw[1] else 1e6
    END <- START + win - 1
  } else if ("POS" %in% names(dt)) {
    START <- dt$POS; END <- dt$POS
  } else {
    stop("No BIN_START/BIN_END/BIN/POS in ", path, ".")
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
  if (nrow(out) == 0) stop("After mapping, no rows remain for ", path)
  out[, `:=`(Start = pmax(1L, Start), End = pmax(Start, End))]
  setcolorder(out, c("Chr","Start","End","Value"))
  out[, Chr_old := NULL]
  setorder(out, Chr, Start, End)
  as.data.frame(out)
}

# centromeres.txt: chr arm start end color_group
read_centromeres <- function(path){
  stopifnot(file.exists(path))
  cdf <- fread(path)
  setnames(cdf, names(cdf), tolower(names(cdf)))
  must <- c("chr","arm","start","end")
  if (!all(must %in% names(cdf)))
    stop("centromeres.txt must have columns: chr, arm, start, end")
  cdf[, `:=`(chr=as.character(chr), arm=tolower(as.character(arm)),
             start=as.integer(start), end=as.integer(end))]
  ce <- cdf[arm=="q", .(Chr=chr, CE_break=start)]
  lens <- cdf[, .(ChromLen=max(end, na.rm=TRUE)), by=chr]; setnames(lens,"chr","Chr")
  ce <- merge(ce, lens, by="Chr", all.x=TRUE)
  ce[, width := pmax(1e5, pmin(2e6, as.numeric(round(0.01*ChromLen))))]
  ce[, `:=`(CE_start=pmax(1L, as.integer(CE_break - width/2)),
            CE_end  =pmin(as.integer(ChromLen), as.integer(CE_break + width/2)))]
  ce[, .(Chr, ChromLen, CE_start, CE_end)]
}

# Build karyotype (5 cols, with centromere)
build_karyotype_with_centromere <- function(over_upper, centromeres){
  ends_from_over <- as.data.table(over_upper)[, .(End_over=max(End, na.rm=TRUE)), by=Chr]
  k <- merge(centromeres, ends_from_over, by="Chr", all=TRUE)
  k[, End := pmax(ChromLen, End_over, na.rm=TRUE)]
  k[!is.finite(End) | is.na(End), End := ChromLen]
  k[, Start := 1L]
  k <- k[, .(Chr, Start, End, CE_start, CE_end)]
  ord <- desired_order[desired_order %in% k$Chr]; if (length(ord)==0) ord <- sort(unique(k$Chr))
  k <- k[match(ord, k$Chr), , drop=FALSE]
  k$CE_start <- pmax(1L, pmin(k$CE_start, k$End-1L))
  k$CE_end   <- pmax(k$CE_start+1L, pmin(k$CE_end, k$End))
  as.data.frame(k)
}

# Discretize into fixed bins → integer codes 1..K; palette aligned 1:K
discretize_values <- function(over_upper, breaks, labels, colors){
  codes <- cut(over_upper$Value, breaks=breaks, labels=FALSE, right=TRUE)
  codes[is.na(codes)] <- 1L
  over_disc <- over_upper
  over_disc$Value <- as.integer(codes)
  list(over=over_disc, palette=colors)
}

# NEW: fill only the *gaps between* existing windows (don’t rebuild a grid)
fill_gaps_between_windows <- function(over_disc, karyo5){
  DT <- as.data.table(over_disc)
  res <- vector("list", length(unique(DT$Chr))); i <- 1L
  for (chr in unique(DT$Chr)) {
    sub <- DT[Chr==chr][order(Start, End)]
    if (nrow(sub)==0) next
    chrEnd <- karyo5$End[karyo5$Chr==chr][1]
    out <- list()

    # gap at start
    if (sub$Start[1] > 1L) {
      out[[length(out)+1]] <- data.table(Chr=chr, Start=1L, End=sub$Start[1]-1L, Value=1L)
    }

    # gaps between consecutive bins
    if (nrow(sub) >= 2) {
      for (j in 1:(nrow(sub)-1)) {
        gap_start <- sub$End[j] + 1L
        gap_end   <- sub$Start[j+1] - 1L
        if (gap_end >= gap_start) {
          out[[length(out)+1]] <- data.table(Chr=chr, Start=gap_start, End=gap_end, Value=1L)
        }
      }
    }

    # gap at end
    if (!is.na(chrEnd) && sub$End[nrow(sub)] < chrEnd) {
      out[[length(out)+1]] <- data.table(Chr=chr, Start=sub$End[nrow(sub)]+1L, End=chrEnd, Value=1L)
    }

    res[[i]] <- rbindlist(c(list(sub), out), use.names=TRUE, fill=TRUE)[order(Start, End)]
    i <- i + 1L
  }
  as.data.frame(rbindlist(res, use.names=TRUE, fill=TRUE))
}

# Legend SVG
make_legend_svg <- function(out_prefix, labels, colors, box_w=40, box_h=20, pad=10, font_cex=1){
  svg_file <- paste0(out_prefix, "_legend.svg")
  n <- length(labels); width <- box_w + 200; height <- n*(box_h+pad)+pad
  grDevices::svg(svg_file, width=width/72, height=height/72)
  op <- par(mar=c(0,0,0,0))
  plot(0,0,type="n", xlim=c(0,width), ylim=c(height,0), axes=FALSE, xlab="", ylab="")
  y <- seq(pad, by=box_h+pad, length.out=n)
  for (i in seq_len(n)) {
    rect(pad, y[i], pad+box_w, y[i]+box_h, col=colors[i], border="grey20", lwd=0.6)
    text(pad+box_w+8, y[i]+box_h/2, labels=labels[i], adj=c(0,0.5), cex=font_cex)
  }
  par(op); dev.off()
  message("Saved legend: ", svg_file)
}

# Plot one dataset
plot_one <- function(over_upper, centro_df, out_prefix){
  karyo5 <- build_karyotype_with_centromere(over_upper, centro_df)
  disc   <- discretize_values(over_upper, bin_breaks, bin_labels, bin_colors)
  overD  <- fill_gaps_between_windows(disc$over, karyo5)  # <- safer gap fill
  pal    <- disc$palette

  # Clamp to chrom ends
  end_map <- setNames(karyo5$End, karyo5$Chr)
  overD$End <- pmin(overD$End, end_map[overD$Chr])

  cat("\n[Diag] Value codes present for", out_prefix, ":", paste(sort(unique(overD$Value)), collapse=", "), "\n")

  svg_file <- paste0(out_prefix, ".svg")
  ideogram(
    karyotype = data.frame(Chr=karyo5$Chr, Start=karyo5$Start, End=karyo5$End,
                           CE_start=karyo5$CE_start, CE_end=karyo5$CE_end,
                           stringsAsFactors=FALSE),
    overlaid  = data.frame(Chr=overD$Chr, Start=overD$Start, End=overD$End, Value=overD$Value,
                           stringsAsFactors=FALSE),
    colorset1 = pal,
    output    = svg_file
  )
  if (!file.exists(svg_file)) {
    for (cand in c("chromosome.svg","ideogram.svg","RIdeogram.svg"))
      if (file.exists(cand)) { file.rename(cand, svg_file); break }
  }
  message("Saved SVG: ", svg_file)
  make_legend_svg(out_prefix, bin_labels, bin_colors)
}

# --------- Run (one dataset at a time) ----------
stopifnot(file.exists(centromere_file))
centro <- read_centromeres(centromere_file)

overA <- read_snpden_upper(fileA)
plot_one(overA, centro, tools::file_path_sans_ext(basename(fileA)))

overB <- read_snpden_upper(fileB)
plot_one(overB, centro, tools::file_path_sans_ext(basename(fileB)))

