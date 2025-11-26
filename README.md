# flip_test
Minimal R tools for diagnosing strand-orientation bias around replication origins using 5′ ends, fragment midpoints, and Monte Carlo flip testing.

# Overview

Accurate strand assignment near replication origins is essential for interpreting BrdU-seq, eSPAN-seq, and replication-coupled ChIP libraries (e.g., Pol2, Pol3).
However, strand calling can shift depending on how fragments are represented:

Mate1 5′ end (read start)

Fragment midpoint (average position of both mates)

Surprisingly, these two representations can produce opposite strand assignments around origins—leading to apparent flips (Pol2) or shrinkage of strand bias (Pol3).
These effects arise not from biology but from fragment-length geometry.

This repository provides two minimal, transparent R functions to analyze and diagnose these effects:

process_sample_minimal()
Extracts properly paired fragments from a BAM file, finds overlaps with replication origins, and computes normal vs flipped strand-match rates.

monte_carlo_flip_test()
A Monte Carlo simulation that quantifies whether observed strand flips can be explained purely by fragment length distribution.

# Features

Robust pairing using makeGAlignmentPairs()

Fragment midpoint computation

Strand-assignment consistency testing

Origin-centered signed-distance calculation

Monte Carlo evaluation of flip probability

Small, interpretable codebase (easy to audit)

# Installation

Clone this repository:

git clone https://github.com/MohammedAlMamun/flip_test.git


Then in R:

library(GenomicAlignments)
library(GenomicRanges)
source("flip_test.R")

# Biological Rationale
Why Pol2 Appears Flipped (Leading-Strand)

Pol2 tracks continuous leading-strand synthesis.

Leading fragments tend to be long.

Long fragments have midpoints that frequently cross the origin center, relative to their 5′ ends.

Observed in data:

Using 5′ ends → expected Watson/Crick bias

Using midpoints → apparent flip (~0.57 flipped)

Key point:
This is not biological — it’s a geometric effect of long fragments.

Why Pol3 Appears Shrunk (Lagging-Strand)

Pol3 tracks lagging-strand synthesis (Okazaki fragments).

Okazaki fragments are short.

Short fragments have midpoints close to their 5′ ends and rarely cross the origin.

Observed in data:

Midpoints show shrinkage (higher normal, lower flipped; ~0.59 / 0.41).

No flip effect.

Key point:
Short fragments → midpoints stay on the same side → strand bias compresses, not reverses.

# Results Summary
Pol2 Example
Metric	Value
Normal-match (midpoints)	0.43
Flipped-match (midpoints)	0.57

Interpretation:
Midpoints artificially assign a higher fraction of Pol2 fragments to the opposite strand.

Pol3 Example
Metric	Value
Normal-match (midpoints)	0.59
Flipped-match (midpoints)	0.41

Interpretation:
Short fragments → reduced crossing of the origin → shrink, not flip.

# Monte Carlo Flip Test

This simulation answers:

Are the observed strand flips (or shrinks) explained purely by fragment length distribution?

Algorithm Steps (simplified)

Keep mate1 5′ positions fixed.

Randomly resample fragment lengths (with replacement).

Recalculate midpoints.

Recompute expected strand orientation (normal vs flipped).

Compare simulated match rates to observed match rates.

Outputs:

observed$normal

observed$flipped

sim$normal

sim$flipped

flip_fraction

p_values

# Example Usage
Load origins (bed: chrom, oriName, oriCenter)
origins <- read.table("data/origins_example.bed", header=TRUE, sep="\t")

Process BAM
res <- process_sample_minimal("Pol2.bam", origins, window = 500)

Run Monte Carlo test
mc_res <- monte_carlo_flip_test(
    mate1_pos_gr   = res$mate1_5p,
    frag_lens_vec  = res$frag_lens_overlap,
    origin_mid_vec = res$origin_mid_overlap,
    strand_obs_vec = res$strand_obs,
    nrep = 2000
)

# Interpretation
flip_fraction

Fraction of simulations where simulated flipped-match > simulated normal-match.

flip_fraction ≈ 1 → data strongly supports flip

flip_fraction ≈ 0 → data strongly supports shrink

flip_fraction ≈ 0.5 → ambiguous / balanced

This matches biological intuition:

Pol2 → flip_fraction ≈ 1

Pol3 → flip_fraction ≈ 0
