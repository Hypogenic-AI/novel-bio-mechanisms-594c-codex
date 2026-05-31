# Discovering Novel Biological Mechanisms from Protein Language Models: Research Plan

## Motivation & Novelty Assessment

### Why This Research Matters
Protein language models (PLMs) are now used to score variants, design proteins, and prioritize experiments, but most interpretability work asks only whether model features recover already-known annotations. If PLMs encode coherent, functionally predictive features that are not explained by current motif databases, those features can become machine-generated biological hypotheses for follow-up mutagenesis or structural study.

### Gap in Existing Work
The gathered literature shows strong evidence that sparse autoencoder (SAE) features in ESM-2 align with Swiss-Prot, GO, PROSITE, ELM, domains, and structural labels. The unresolved gap is the next step: quantitatively triaging unannotated latent features by functional relevance, not just calling them "unknown" because annotations are incomplete.

### Our Novel Contribution
This project tests a compact feature-to-hypothesis workflow: identify residue-level SAE features from a real PLM, measure whether they predict ProteinGym mutation sensitivity, explicitly screen them against known ELM/PROSITE motifs, and run internal feature ablations to ask whether candidate features causally change ESM-2 variant scores.

### Experiment Justification
- Experiment 1: Known-motif alignment screen. This establishes which SAE features are already explainable by motif databases and gives a control for novelty claims.
- Experiment 2: ProteinGym functional association. This tests whether latent features mark residues where mutations affect measured fitness/function.
- Experiment 3: Internal feature intervention. This tests whether top candidate features influence ESM-2 masked-marginal variant predictions rather than only correlating with labels.
- Experiment 4: Baseline comparison. Raw hidden dimensions and shuffled controls are required to show that sparse features are more useful than arbitrary model coordinates or sequence-position artifacts.

## Research Question
Do ESM-2 sparse latent features identify functionally important protein positions that are not aligned with known ELM/PROSITE motif annotations, and do interventions on those features alter ESM-2 mutation-effect predictions?

## Background and Motivation
InterPLM, InterProt, MotifAE, and related studies show that PLM sparse features recover known motifs, domains, binding sites, localization signals, and family-specific patterns. The stronger scientific claim is not merely that features are interpretable, but that model internals can propose new biological mechanisms. ProteinGym provides an empirical fitness/function benchmark, while ELM and PROSITE provide conservative known-motif controls.

## Hypothesis Decomposition
- H1: SAE features have stronger associations with ProteinGym residue-level mutation sensitivity than raw ESM hidden dimensions and shuffled controls.
- H2: Some strongly function-associated SAE features have low overlap with known ELM/PROSITE motif scans.
- H3: Ablating a top unannotated candidate feature changes ESM-2 masked-marginal variant scores more at high-activation positions than low-activation positions.
- Alternative explanations: candidate features may reflect solvent exposure, conservation, assay-specific noise, motif-database false negatives, or generic sequence statistics rather than a novel mechanism.

## Proposed Methodology

### Approach
Use real ESM-2-8M activations and pretrained InterPLM SAE features. For tractability, analyze short ProteinGym substitution assays with complete single-mutant coverage. Aggregate DMS scores per residue to obtain position sensitivity, compute masked-marginal ESM scores, extract SAE and raw hidden-state activations at wild-type residues, and compare feature/sensitivity correlations. Scan ELM and PROSITE patterns on the target sequences and a Swiss-Prot sample to quantify known-motif alignment.

### Experimental Steps
1. Validate datasets and select ProteinGym assays with target length under 256 and sufficient single substitutions.
2. Load ESM-2-8M and InterPLM layer-4 SAE on GPU; use batch sizes around 128 because RTX A6000 GPUs have 49 GB memory.
3. For each selected assay, compute residue-level DMS sensitivity and ESM masked-marginal sensitivity.
4. Extract layer-4 hidden states and SAE features for target sequences.
5. Compute per-feature Spearman correlations with DMS sensitivity and ESM sensitivity; compute the same for raw hidden dimensions and shuffled-position controls.
6. Scan ELM/PROSITE motifs and classify top features as known-aligned or low-alignment candidates.
7. Ablate decoder directions for the strongest low-alignment candidates during ESM masked-position scoring and compare score changes at high- versus low-feature residues.
8. Save raw outputs, summaries, plots, and candidate tables.

### Baselines
- Raw ESM-2 layer-4 hidden dimensions, tested with the same correlation protocol.
- Shuffled SAE activations within each assay, preserving marginal distributions but breaking positional structure.
- ESM masked-marginal scores as the standard zero-shot functional prediction baseline.
- Known-motif overlap from ELM/PROSITE as an annotation-alignment control.

### Evaluation Metrics
- Spearman correlation between feature activation and residue-level DMS sensitivity.
- Spearman correlation between feature activation and ESM masked-marginal sensitivity.
- Empirical permutation p-values with Benjamini-Hochberg FDR correction across features.
- AUROC/enrichment for known motif-position labels when motif scans yield positives.
- Intervention effect size: difference in absolute masked-marginal score change between high- and low-activation residues.

### Statistical Analysis Plan
Use paired comparisons at the assay level for best SAE feature versus best raw hidden dimension and versus shuffled controls. Use Spearman correlation because DMS scales differ across assays and monotonic association is the target. Use permutation tests for feature significance, Benjamini-Hochberg correction for multiple feature tests, and bootstrap confidence intervals for aggregate median differences. Report effect sizes alongside p-values.

## Expected Outcomes
Support for the hypothesis would be: top SAE features outperform raw/shuffled baselines, at least one top feature has low ELM/PROSITE alignment, and feature ablation causes larger variant-score changes at high-activation residues. Refutation would be: SAE associations do not exceed baselines, all functional features map to known motifs, or ablations do not change mutation scores.

## Timeline and Milestones
- Resource review and planning: complete before implementation.
- Environment setup and dependency installation: 10-20 minutes.
- Pipeline implementation: 60-90 minutes.
- Experiment execution: 60-90 minutes, adjusted to available GPU time.
- Statistical analysis, visualizations, and documentation: 45-75 minutes.

## Potential Challenges
- Pretrained SAE loading may fail due to dependency or Hugging Face issues; fallback is the smaller `plm-sae` ESM2-8M layer-3 SAE or a small locally trained SAE.
- PROSITE pattern syntax is not identical to Python regex; the parser will support common pattern syntax and report unsupported patterns.
- DMS assays differ in score direction and scale; analysis will use within-assay rank and z-score summaries.
- Motif databases are incomplete, so "low alignment" is interpreted as candidate novelty, not proof of unknown biology.

## Success Criteria
The research succeeds if it produces a reproducible pipeline, runs real ESM-2/SAE experiments on local ProteinGym data, reports baseline-controlled statistics, and outputs a ranked list of candidate latent features with explicit known-annotation checks and limitations.
