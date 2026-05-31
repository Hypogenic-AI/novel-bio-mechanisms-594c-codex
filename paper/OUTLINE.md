# Outline: Baseline-Controlled PLM Feature Triage

## Title
Baseline Controls Temper Novel-Mechanism Claims from Sparse Protein Language Model Features

## Core Claim
Sparse ESM-2 SAE features often correlate with DMS-sensitive residues, but in this pilot they do not reliably outperform raw hidden dimensions or shuffled sparse controls. The defensible contribution is a reproducible triage workflow and one exploratory causal candidate.

## Evidence Map
- Main association: median best SAE absolute Spearman 0.547 across 16 assays.
- Baselines: raw hidden dimensions median 0.544; shuffled SAE 95th percentile median 0.543.
- Statistics: SAE vs raw Wilcoxon one-sided p = 0.798; SAE vs shuffled p = 0.217; bootstrap median difference CIs cross zero.
- Candidates: three low target-motif and low Swiss-Prot motif-alignment feature-assay pairs.
- Intervention: feature 6419 in YNZC_BACSU_Tsuboyama_2023_2JVD gives high-minus-low activation effect 0.0223 log-prob units, Mann-Whitney p = 0.00495.

## Section Plan
- Introduction: motivate PLM interpretability for biology; state gap between recovering known annotations and claiming novel mechanisms; preview conservative findings.
- Related work: PLMs and ProteinGym; SAE interpretability in PLMs; motif databases as annotation controls; causal intervention need.
- Methodology: datasets, model, SAE, residue sensitivity, correlations, baselines, motif scans, interventions, reproducibility.
- Results: aggregate comparison, candidate table, intervention table, interpretation of controls.
- Discussion: strongest support is workflow; novelty claim remains unproven; limitations and next experiments.
- Conclusion: summarize cautious contribution and future validation path.
