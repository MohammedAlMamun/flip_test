library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)

process_sample_minimal <- function(bam, origins, window = 500) {
  
  cat("\n=== Minimal Fragment-Orientation Processing ===\n")
  cat("BAM:", bam, "\n")
  
  # ---- Step 1: Basic BAM counts ----
  cat("\n-- Basic counts --\n")
  print(countBam(bam))
  print(countBam(bam, param = ScanBamParam(flag = scanBamFlag(isPaired = TRUE))))
  print(countBam(bam, param = ScanBamParam(flag = scanBamFlag(isProperPair = TRUE))))
  
  # ---- Step 2: Read alignments + make pairs ----
  cat("\n-- Reading alignments --\n")
  
  param <- ScanBamParam(what = scanBamWhat(),
                        flag = scanBamFlag(isUnmappedQuery = FALSE))
  
  gal_chunk <- readGAlignments(bam, param = param, use.names = TRUE)
  cat("Alignments read:", length(gal_chunk), "\n")
  
  # Always pair by name (robust)
  gal_pairs <- makeGAlignmentPairs(gal_chunk)
  cat("Paired fragments:", length(gal_pairs), "\n")
  
  # ---- Step 3: clean pairs (remove chrM, NA) ----
  keep_first  <- !is.na(seqnames(gal_pairs@first))  & seqnames(gal_pairs@first)  != "chrM"
  keep_last   <- !is.na(seqnames(gal_pairs@last))   & seqnames(gal_pairs@last)   != "chrM"
  keep_pairs <- keep_first & keep_last
  
  gal_clean <- gal_pairs[keep_pairs]
  cat("Kept fragments:", length(gal_clean), "\n")
  
  # ---- Step 4: Convert to fragment GRanges ----
  frag_gr <- granges(gal_clean)
  cat("Fragments after cleaning:", length(frag_gr), "\n")
  
  # ---- Step 5: Build origin windows ----
  origins_gr <- GRanges(
    seqnames = origins$chrom,
    ranges   = IRanges(start = origins$oriCenter, width = 1),
    oriName  = origins$oriName
  )
  
  origin_win <- GRanges(
    seqnames = seqnames(origins_gr),
    ranges   = IRanges(start(origins_gr) - window, start(origins_gr) + window)
  )
  
  # ---- Step 6: Count origin overlaps ----
  hits <- findOverlaps(frag_gr, origin_win, ignore.strand = TRUE)
  idx <- unique(queryHits(hits))
  
  cat("\nFragments overlapping origins:", length(idx),
      "out of", length(frag_gr), "\n")
  
  # ---- Step 7: Compute signed distances ----
  frag_mid <- floor((start(frag_gr) + end(frag_gr)) / 2)
  
  frag_origin <- rep(NA_integer_, length(frag_gr))
  frag_origin[queryHits(hits)] <- subjectHits(hits)
  
  frag_mid_in <- frag_mid[idx]
  origin_mid_in <- start(origins_gr)[frag_origin[idx]]
  signed_dist <- frag_mid_in - origin_mid_in
  
  cat("\nSigned distance summary:\n")
  print(summary(signed_dist))
  
  # ---- Step 8: Orientation classification ----
  strand_obs <- as.character(strand(gal_clean@first)[idx])
  pos_sign <- ifelse(signed_dist < 0, "left", "right")
  
  expected_normal <- ifelse(pos_sign == "left", "+", "-")
  expected_flipped <- ifelse(pos_sign == "left", "-", "+")
  
  match_normal  <- mean(strand_obs == expected_normal)
  match_flipped <- mean(strand_obs == expected_flipped)
  
  cat("\nNormal-match rate:", match_normal, "\n")
  cat("Flipped-match rate:", match_flipped, "\n")
  
  # Output final list
  out <- list(
    frag_gr = frag_gr,
    idx = idx,
    strand_obs = strand_obs,
    signed_dist = signed_dist,
    normal_rate = match_normal,
    flipped_rate = match_flipped,
    frag_lens_overlap = width(frag_gr)[idx],
    origin_mid_overlap = origin_mid_in,
    mate1_5p = GRanges( seqnames = seqnames(gal_clean@first[idx]),
                        ranges = IRanges(start = start(gal_clean@first[idx]), width = 1),
                        strand = strand(gal_clean@first[idx]))
  )
  
  # Add attributes carrying the observed results
  attr(out$mate1_5p, "match_normal")  <- match_normal
  attr(out$mate1_5p, "match_flipped") <- match_flipped
  
  return(out)
  
}

# Corrected monte carlo matching your minimal pipeline behaviour
monte_carlo_flip_test <- function(mate1_pos_gr,        # 1bp GRanges at mate1 5' positions
                                  frag_lens_vec,       # numeric vector of fragment lengths (same order)
                                  origin_mid_vec,      # numeric vector of origin centers (same order)
                                  strand_obs_vec,      # character vector ("+" or "-") observed for mate1 (same order)
                                  nrep = 2000,
                                  sample_with_repl = TRUE,
                                  plot_it = TRUE,
                                  verbose = TRUE) {
  
  # basic checks
  n <- length(mate1_pos_gr)
  stopifnot(length(frag_lens_vec)   == n,
            length(origin_mid_vec)  == n,
            length(strand_obs_vec)  == n)
  
  # ensure numeric vectors
  mate1_starts <- as.integer(start(mate1_pos_gr))
  frag_lens_vec <- as.integer(frag_lens_vec)
  origin_mid_vec  <- as.integer(origin_mid_vec)
  strand_obs_vec  <- as.character(strand_obs_vec)
  
  # Pull observed values directly from wrapper
  match_normal_obs  <- attr(mate1_pos_gr, "match_normal")
  match_flipped_obs <- attr(mate1_pos_gr, "match_flipped")
  
  if (verbose) {
    message("Observed normal-match = ", signif(match_normal_obs, 6),
            "; observed flipped-match = ", signif(match_flipped_obs, 6),
            "  (n = ", n, ")")
  }
  
  # Prepare storage for sims
  sim_norm <- numeric(nrep)
  sim_flip <- numeric(nrep)
  
  # run Monte Carlo
  for (r in seq_len(nrep)) {
    # sample fragment lengths
    if (sample_with_repl) {
      sampled_Ls <- sample(frag_lens_vec, size = n, replace = TRUE)
    } else {
      sampled_Ls <- sample(frag_lens_vec, size = n, replace = FALSE)
    }
    
    # simulated midpoints using sampled lengths
    sim_mid <- mate1_starts + floor((sampled_Ls - 1) / 2)
    
    # simulated signed distances relative to the same origins
    sim_signed <- sim_mid - origin_mid_vec
    sim_pos_sign <- ifelse(sim_signed < 0, "left", "right")
    
    sim_expected_normal  <- ifelse(sim_pos_sign == "left", "+", "-")
    sim_expected_flipped <- ifelse(sim_pos_sign == "left", "-", "+")
    
    sim_norm[r] <- mean(strand_obs_vec == sim_expected_normal)
    sim_flip[r] <- mean(strand_obs_vec == sim_expected_flipped)
  }
  
  # empirical p and flip fraction
  p_emp_norm_ge <- mean(sim_norm >= match_normal_obs)
  flip_fraction <- mean(sim_flip > sim_norm)
  
  if (verbose) {
    message("Monte Carlo done (nrep=", nrep, "). flip_fraction = ", signif(flip_fraction, 6))
  }
  
  # # optional plots (histogram of sim_norm and flip fraction bar)
  # if (plot_it) {
  #   op <- par(no.readonly = TRUE)
  #   par(mfrow = c(1, 2))
  #   hist(sim_norm, breaks = 40, main = "Simulated normal-match rates",
  #        xlab = "normal-match (sim)", col = "grey80", border = "white")
  #   abline(v = match_normal_obs, col = "red", lwd = 2)
  #   legend("topright",
  #          legend = c(paste0("obs=", round(match_normal_obs, 3)),
  #                     paste0("p(>=obs)=", signif(p_emp_norm_ge, 3))),
  #          bty = "n")
  #   
  #   barplot(height = c(flip_fraction, 1 - flip_fraction),
  #           names.arg = c("flipped>normal", "not flipped>normal"),
  #           col = c("red", "darkgreen"), main = "Flip fraction", ylab = "Proportion")
  #   par(op)
  # }
  
  list(
    observed = list(normal = match_normal_obs, flipped = match_flipped_obs),
    sim = list(normal = sim_norm, flipped = sim_flip),
    p_values = list(p_norm_ge = p_emp_norm_ge),
    flip_fraction = flip_fraction,
    n_frag = n
  )
}

#-------------------------
# Run for samples
#-------------------------

Sample_Name <- "Rev1-myc-0422-tScript"

bam <- Sys.glob(paste0("~/Desktop/", Sample_Name, "/Bam/*.bam"))[5]
origins <- read.table(paste0("~/Desktop/", Sample_Name, "/Peaks/", Sample_Name, "_Primary_Peaks.bed"), header = TRUE)

res <- process_sample_minimal(bam, origins, window = 500)
mc_res <- monte_carlo_flip_test(mate1_pos_gr   = res$mate1_5p,
                                frag_lens_vec  = res$frag_lens_overlap,
                                origin_mid_vec = res$origin_mid_overlap,
                                strand_obs_vec = res$strand_obs,
                                nrep = 1000,
                                plot_it = TRUE,
                                verbose = TRUE)

