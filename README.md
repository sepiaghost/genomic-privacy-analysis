# Differential Privacy for Genomic Data Sharing

## Overview

This project demonstrates and quantifies the **privacy-utility tradeoff** in
genomic data sharing using differential privacy — the same class of
technique real institutions use (or should use) when releasing gene
expression data for research while protecting individual participants from
re-identification.

**The question this project answers:** How much statistical accuracy do you
have to give up to protect participants from having their presence in a
study detected — and how much protection do you actually get in return?

## Why this project

Genomic data is uniquely re-identifying: an individual's expression profile
can be distinctive enough that even "anonymized" release of aggregate
statistics can leak whether they participated in a study (a **membership
inference attack**). This project implements the standard defense
(the **Laplace mechanism**, the foundational tool of differential privacy)
and empirically measures both sides of the tradeoff on real RNA-seq data,
rather than treating it as a purely theoretical concern.

This sits at the intersection of my cybersecurity background and
computational genomics — most bioinformatics portfolios demonstrate the
biology side; this one demonstrates the security side computationally, on
real data.

## Data source

- **Dataset:** [`airway`](https://bioconductor.org/packages/release/data/experiment/html/airway.html) —
  a built-in Bioconductor RNA-seq dataset: airway smooth muscle cells,
  treated vs. untreated with dexamethasone (n=8, 4 vs. 4)
- Ships directly inside the `airway` R package — **no external download
  required**, unlike TCGA-scale projects
- Chosen deliberately for this project: small enough to re-run differential
  expression dozens of times across privacy levels in seconds rather than
  hours, and it's the same dataset used in DESeq2's own official
  documentation, so it's a recognized benchmark

## Pipeline

| Script | Purpose |
|---|---|
| `scripts/01_load_data.R` | Loads the airway dataset |
| `scripts/02_ground_truth_de.R` | Runs real (non-private) DESeq2 differential expression — the baseline everything else is measured against |
| `scripts/03_dp_mechanism.R` | Implements the Laplace mechanism; generates privatized count matrices at 7 privacy levels |
| `scripts/04_utility_tradeoff.R` | Re-runs DE on each privatized matrix; measures Jaccard/precision/recall/F1 vs. ground truth |
| `scripts/05_membership_inference.R` | Simulates a membership inference attack; measures attack success rate at each privacy level |
| `scripts/06_visualization.R` | Produces the two headline tradeoff figures |

Run in order:

```r
source("scripts/01_load_data.R")
source("scripts/02_ground_truth_de.R")
source("scripts/03_dp_mechanism.R")
source("scripts/04_utility_tradeoff.R")
source("scripts/05_membership_inference.R")
source("scripts/06_visualization.R")
```

Total runtime: a few minutes, thanks to the small dataset — no multi-hour
downloads like a TCGA-scale project.

## Method notes 

- **Sensitivity via clipping:** the Laplace mechanism requires bounding how
  much one individual could possibly change the output. Every count is
  clipped to the 95th percentile before noise is added — otherwise a single
  extreme value could break the privacy guarantee.
- **Scoping honesty:** noise is applied independently per gene/sample as an
  illustration of the mechanism and its tradeoff. A production-grade private
  release of an entire count matrix requires formal **privacy budget
  composition** across every value released — that accounting is out of
  scope here. This project demonstrates the mechanism and measures its
  effect, not a deployment-ready private data release system.
- **Membership inference simplification:** the simulated attacker is
  assumed to already know both candidate "true" means (a common, somewhat
  strong assumption in membership inference research, used here to make the
  privacy/no-privacy contrast clean and interpretable) — this represents a
  strong-attacker scenario, useful for demonstrating the mechanism's
  protective effect at its edges.

## Results

### Privacy Utility Tradeof

Applying the Laplace mechanism at seven privacy levels (ε = 0.1 to 10, plus 
ε = ∞ as the no-privacy baseline) and re-running differential expression on 
each privatized dataset produced the following overlap with the ground-truth 
DE gene set (878 genes at ε = ∞):

| ε   | # Significant Genes | Jaccard | Precision | Recall | F1     |
|-----|---------------------|---------|-----------|--------|--------|
| 0.1 | 4,897               | 0.0114  | 0.0133    | 0.0740 | 0.0225 |
| 0.5 | 4,515               | 0.0097  | 0.0115    | 0.0592 | 0.0193 |
| 1   | 4,399               | 0.0057  | 0.0068    | 0.0342 | 0.0114 |
| 2   | 3,665               | 0.0062  | 0.0076    | 0.0319 | 0.0123 |
| 5   | 2,236               | 0.0104  | 0.0143    | 0.0365 | 0.0206 |
| 10  | 922                 | 0.0262  | 0.0499    | 0.0524 | 0.0511 |
| ∞   | 878                 | 1.0000  | 1.0000    | 1.0000 | 1.0000 |

Recovery of the true DE gene set is poor across nearly the entire privacy 
range tested — F1 never exceeds 0.051, even at the weakest privacy setting 
(ε = 10). Notably, utility does **not** degrade monotonically with ε: F1 is 
actually higher at ε = 0.1 (0.0225) than at ε = 1–2 (~0.011–0.012), only 
recovering meaningfully at ε = 10. This non-monotonic behavior reflects how 
noise interacts with a fixed significance threshold — at very high noise, 
nearly every gene becomes "significant" by chance, artificially inflating 
recall despite the result being scientifically meaningless (note the 4,897 
genes flagged as significant at ε = 0.1, over five times the true count).

### Membership Inference Attack

| ε   | Attack Accuracy |
|-----|------------------|
| 0.1 | 0.506            |
| 0.5 | 0.514            |
| 1   | 0.494            |
| 2   | 0.511            |
| 5   | 0.521            |
| 10  | 0.480            |
| ∞   | 1.000            |

Without privacy protection, an attacker can perfectly determine whether a 
given sample was part of the release (accuracy = 1.0). The moment any 
Laplace noise is introduced — even at the weakest tested privacy budget 
(ε = 10) — attack accuracy collapses to the random-guessing baseline 
(0.48–0.52), and stays there across the entire range down to ε = 0.1.

### Key Finding

Privacy protection against membership inference is cheap: even light noise 
defeats the attack almost completely, and accuracy never meaningfully rises 
above chance at any tested ε. Statistical utility, by contrast, is 
expensive and behaves unpredictably: recovery of true differentially 
expressed genes stays poor (F1 < 0.06) across the entire tested range, and 
does not improve smoothly as ε increases — noise-inflated false positives 
can make weaker privacy settings look artificially "better" by recall alone 
without being scientifically meaningful. This asymmetry — privacy is easy 
to achieve, but preserving genuine scientific utility under it is hard and 
non-trivial to measure correctly — is the central finding of this project.


![Privacy-utility tradeoff](figures/privacy_utility_tradeoff.png)
![Membership inference attack accuracy](figures/membership_inference_attack.png)

## Tools & packages

- R / Bioconductor
- `airway` — built-in benchmark RNA-seq dataset
- `DESeq2` — differential expression (both ground-truth and privatized runs)
- `ggplot2` — visualization

## Project structure

```
genomic-privacy-analysis/
├── README.md
├── .gitignore
├── LICENSE
├── data/processed/       # gitignored intermediate files
├── scripts/               # numbered pipeline
├── results/                # tradeoff and attack CSVs
└── figures/                 # tradeoff figures
```

## License

MIT — see LICENSE file.
