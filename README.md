# Sleep as the brain's "bowhead mode": sleep loss lowers CIRBP

**Hypothesis.** Sleep-dependent DNA repair in neurons is partly temperature-gated through CIRBP. CIRBP is the cold-inducible RNA-binding protein that gives the long-lived bowhead whale its enhanced DNA double-strand-break repair, and PARP1 recruits it to DNA breaks. The brain cools during sleep; waking warms it and lowers *Cirbp*.

This repository re-analyses public gene-expression data to test the first link of that chain, and lists the experiments needed to test the rest.

## Main result

Sleep loss lowered *Cirbp* in all three mouse brain datasets tested (3 labs, 2 brain regions).

| Dataset | Tissue | n (SD vs control) | *Cirbp* log2 change | p | Rank among genes |
|---|---|---|---|---|---|
| GSE166831 | Hippocampus (RNA-seq) | 9 vs 9 | -0.51 | 1.3e-7 | 3 / 23,295 |
| GSE33302 | Hippocampus (microarray) | 8 vs 9 | -0.68 | 4.4e-5 | 75 / 19,077 |
| GSE113754 | Frontal cortex (RNA-seq) | 5 vs 5 | -1.06 | 6.7e-6 | 76 / 16,900 |
| **Pooled** (random effects) | | | **Hedges' g -4.85** (95% CI -7.36 to -2.33) | **0.00016** | |

Each dataset passed a label check: known sleep-loss genes (*Arc*, *Fos*, *Nr4a1*, *Egr1*, *Homer1*) rose as expected.

## What did not hold, or could not be tested

- *Rbm3*, the second cold-inducible gene, was inconsistent: down in hippocampus, up in cortex.
- *Cirbp* expression did not track DNA-repair gene expression within groups (T2).
- In human blood (GSE39445; 22 people, about a week of 6 h vs 10 h sleep opportunity), *CIRBP* did not change (p = 0.61). *PARP1* was lower after short sleep (17/22 people, p = 0.005); this is exploratory.
- **Exercise:** *Cirbp* was absent from the available rat file and from the mouse array, so it was not testable.
- **Recovery sleep:** not testable with GSE140345, because its time points do not pair sleep loss and controls.

## How to reproduce

1. Open `notebooks/sleep_cirbp_analysis.ipynb` in Google Colab.
2. Choose **Runtime > Run all**. Takes about 15-20 minutes and downloads from NCBI GEO.
3. If NCBI is busy and a dataset fails, run again later.

`notebooks/sleep_cirbp_analysis_run_2026-10-31.ipynb` is the executed run with all outputs. `notebooks/sleep_cirbp_analysis_single_cell.py` is the same code for pasting into one cell.

Methods:
- RNA-seq counts are converted to log2 CPM (genes with at least 10 counts in at least 3 samples).
- Microarray probes are collapsed to genes (highest-mean probe).
- Welch t-test and bias-corrected Hedges' g; DerSimonian-Laird random-effects pooling.
- Co-expression is tested against 1,000 random gene sets.
- Sample groups are read from the authors' count-file column names where GEO order differs.

## Status and limits

This is a hypothesis with supporting re-analysis, not a demonstration that CIRBP drives sleep-dependent repair. Groups are small (5-12 per arm), and the cortex effect size is inflated by n = 5. The decisive tests are wet-lab:
1. Neuronal gamma-H2AX after sleep vs sleep loss in *Cirbp* knockout mice.
2. Sleep loss at cool vs normal room temperature.
3. Human cells at 36 vs 37 C (CIRBP protein, NHEJ reporter).

## Data and credit

All data are public. Please cite the original studies:
- Gaine et al. 2021, *Mol Brain* 14:125 (GSE166831)
- GSE33302 (hippocampus microarray, sleep deprivation)
- GSE113754 (frontal cortex, sleep deprivation; wild-type samples only)
- Moller-Levet et al. 2013, *PNAS* 110:E1132 (GSE39445)
- Firsanov et al. 2025, *Nature* 648:717 (bowhead CIRBP)
- Hoekstra et al. 2019, *eLife* 8:e43400
- Chen et al. 2018, *PNAS* 115:E1759
- Zada et al. 2021, *Mol Cell*

## Author

Riddika D , B.Tech Biotechnology, Sathyabama Institute of Science and Technology, Chennai, India.

Code: MIT licence. Text and figures: CC BY 4.0.
