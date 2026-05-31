# Literature Review: Discovering Novel Biological Mechanisms from Protein Language Models

## Review Scope

### Research Question

Can latent features, neurons, or circuits inside protein language models reveal biologically meaningful mechanisms, especially motifs, interaction patterns, or constraints that are not already captured by standard annotations?

### Inclusion Criteria

- Protein language model interpretability, mechanistic interpretability, SAE/transcoder, neuron, circuit, attention, CAV, or attribution methods.
- Protein function, fitness, motif, domain, structure, or mutation-effect prediction tasks.
- Papers with actionable datasets, code, metrics, or experimental protocols.

### Exclusion Criteria

- General biological ML papers without PLM internals.
- Protein prediction papers with no interpretability or latent-feature analysis.
- Papers where no accessible PDF, abstract, or implementation was available.

### Search Log

| Date | Query / Source | Results | Notes |
|---|---|---:|---|
| 2026-05-31 | paper-finder: "protein language model interpretability latent features circuits protein function motifs" | 151 | 26 high-relevance, 28 medium-relevance; downloaded core set after abstract screening. |
| 2026-05-31 | Web/arXiv: "protein language model sparse autoencoder" | 10+ | Found InterPLM, ProtSAE, InterProt, SAEFold, low-N SAE, circuit tracing. |
| 2026-05-31 | Web/Hugging Face: ProteinGym, TAPE, Swiss-Prot, PROSITE, ELM | 4 dataset families | Selected datasets with direct relevance to feature alignment and function prediction. |

## Research Area Overview

Protein language models (PLMs), especially ESM-2 and related transformer models, learn sequence representations that support structure, function, localization, and mutation-effect prediction. The current interpretability literature has moved from attention visualization and supervised probing toward sparse feature decompositions. Sparse autoencoders (SAEs) and transcoders are now the dominant tools because they address superposition: biological concepts appear to be distributed across dense PLM activations rather than represented by single neurons.

The research hypothesis is well aligned with the most recent literature. InterPLM, Gujral et al., and InterProt all show that SAE features align with known protein families, binding sites, motifs, domains, GO terms, and structural properties. More importantly, several papers report coherent high-activation feature clusters that do not map cleanly to existing annotations, which is the opening for novel mechanism discovery. The main unresolved challenge is turning "unannotated but coherent" latent features into experimentally credible biological hypotheses rather than annotation artifacts.

## Key Papers

### InterPLM: Discovering Interpretable Features in Protein Language Models via Sparse Autoencoders

- Authors: Elana Simon, James Zou
- Year: 2024 preprint; Nature Methods 2025
- Methodology: Trains SAEs over ESM-2 amino-acid embeddings across layers, compares SAE features to raw neurons, associates binary feature activations with Swiss-Prot concepts, uses automated LLM feature descriptions, and demonstrates steering.
- Datasets: UniRef50/Swiss-Prot sequences, Swiss-Prot annotations, AlphaFold structures.
- Results: SAE features recover far more concepts than neurons; reported up to 2,548 interpretable features per layer and associations across binding sites, domains, motifs, active sites, disulfides, phosphorylation, disorder, and targeting signals.
- Code: `code/interPLM/`, pretrained SAEs on Hugging Face.
- Relevance: Primary blueprint for feature discovery and missing-annotation hypothesis generation.

### Sparse autoencoders uncover biologically interpretable features in protein language model representations

- Authors: Onkar Gujral, Mihir Bafna, Eric Alm, Bonnie Berger
- Year: 2025
- Methodology: Applies SAEs and transcoders to protein-level and amino-acid-level ESM2 representations; evaluates GO enrichment and LLM-assisted interpretation.
- Datasets: ESM2 representations, Swiss-Prot/UniProt metadata, GO terms, protein families.
- Results: Sparse features are more interpretable than raw PLM neurons; features align with GO terms and specific families/functions including metabolic and sensory pathways.
- Code: `code/Rep_SAEs_PLMs/`
- Relevance: Strong support for unsupervised feature extraction and automated interpretation, but requires careful control for GO annotation biases.

### From Mechanistic Interpretability to Mechanistic Biology

- Authors: Etowah Adams, Li Bai, Minji Lee, Yiyang Yu, Mohammed AlQuraishi
- Year: 2025
- Methodology: Trains SAEs on ESM-2 residual streams, classifies features as generic/family-specific, uses linear probes to connect features to thermostability and localization, and visualizes uninterpreted candidate features.
- Datasets: PLM activations, protein families, property labels.
- Results: Demonstrates a practical workflow from latent feature to biological hypothesis; emphasizes that unknown features should be treated as testable hypotheses.
- Code: `code/interprot/`
- Relevance: Most directly matches this project goal: "study the model to study biology."

### MotifAE Reveals Functional Motifs from Protein Language Model

- Authors: Chao Hou, Di Liu, Yufeng Shen
- Year: 2025
- Methodology: Adds a local smoothness/coherence loss to SAE training over ESM2 embeddings, extracts PSSMs from latent features, evaluates against ELM motifs and CATH domains, and aligns selected features with ProteinGym/DMS stability data.
- Datasets: AlphaFold cluster representative sequences, ELM, CATH, ProteinGym DMS.
- Results: MotifAE improves known motif recovery over standard SAE, with reported ELM median AUROC 0.88 versus 0.80 for SAE; features can support stability landscape prediction and in silico design.
- Code: No official repo found during search.
- Relevance: Best motif-specific method and a strong baseline for novel motif discovery.

### Protein Circuit Tracing via Cross-layer Transcoders

- Authors: Darin Tsui, Kunal Talreja, Daniel Saeedi, Amirali Aghazadeh
- Year: 2026
- Methodology: Uses cross-layer transcoders to approximate ESM2 computation across layers and extract compact circuits for family classification and function prediction.
- Datasets: Swiss-Prot-style family/function labels and PLM activations.
- Results: Reports 82-89% recovery of original model performance and compact circuits using less than 1% of latent space while retaining substantial task accuracy.
- Code/models: model release found on Hugging Face; no GitHub repo found in this run.
- Relevance: Important extension beyond independent SAE features; useful if experiments prioritize causal circuits over static features.

### ProteinGym

- Authors: Pascal Notin et al.
- Year: 2023
- Methodology: Large benchmark of DMS and clinical mutation-effect datasets with zero-shot and supervised evaluation.
- Datasets: 250+ DMS assays and clinical variants; local Hugging Face snapshot downloaded.
- Results: Defines standard metrics and baselines across mutation types, taxa, MSA depth, and assay classes.
- Code: `code/ProteinGym/`
- Relevance: Primary benchmark for asking whether discovered features influence function/fitness predictions.

## Common Methodologies

- Sparse decomposition: Standard ReLU/Top-K/BatchTopK/Matryoshka SAEs expand dense PLM activations into sparse features. Evaluation usually includes reconstruction error, fraction variance explained, L0 sparsity, dead features, and downstream fidelity.
- Transcoders and circuits: Transcoders learn sparse approximations of transformations between layers. Cross-layer transcoders are more appropriate for "circuits" because they preserve computation across layers.
- Annotation association: Feature activations are binarized or scored and compared to Swiss-Prot, GO, ELM, PROSITE, CATH, active-site, binding-site, transmembrane, localization, and domain annotations using F1, AUROC, enrichment tests, and manual inspection.
- Automated interpretation: Several papers prompt LLMs with top-activating and inactive proteins/residues. This scales feature labeling, but must be validated with held-out proteins and known annotations.
- Causal steering: Feature or neuron activations are modified during generation/prediction to test whether they causally affect protein properties. This is stronger evidence than correlation but can introduce off-manifold artifacts.
- Attribution baselines: Attention analysis, integrated gradients, and concept activation vectors remain useful comparators, especially for residue-level localization.

## Standard Baselines

- Raw ESM neurons or dimensions versus SAE features.
- Attention-head relevance and contact/motif attention from ProVis/BERTology.
- Integrated gradients over input and internal transformer heads from `xai-proteins`.
- Linear probes over PLM embeddings and SAE latents for known labels.
- Concept activation vectors for known motif localization.
- ESM masked-marginal or pseudo-log-likelihood scoring on ProteinGym.
- ProteinGym baselines: ESM-1v/ESM-2, MSA Transformer, EVE, Tranception, ProteinMPNN, ESM-IF1, and newer leaderboard models.

## Evaluation Metrics

- Feature quality: reconstruction MSE, fraction variance explained, L0 sparsity, dead feature rate, activation frequency.
- Concept alignment: AUROC, AUPRC, F1, precision/recall, enrichment p-values/FDR, held-out validation of feature labels.
- Residue localization: residue-level AUROC/F1 against ELM/PROSITE/Swiss-Prot feature coordinates, distance to active/binding sites, structural clustering.
- Functional influence: Spearman correlation with DMS scores, MSE for supervised fitness, AUC/MCC for binary pathogenic/benign or fit/unfit labels.
- Circuit fidelity: percentage of original PLM performance recovered by replacement model/circuit; performance retained by sparse circuit subsets.
- Novelty triage: no/low overlap with existing annotations, coherence across homologous but evolutionarily distinct proteins, structural colocalization, and predictive value on held-out DMS/clinical data.

## Datasets in the Literature

- Swiss-Prot/UniProtKB: curated protein sequences and rich annotations; used by InterPLM and Rep_SAEs_PLMs for feature association.
- UniRef50/UniRef90: large-scale pretraining or SAE training sequences; useful but large.
- GO/GOA: function labels and enrichment analysis for protein-level feature interpretation.
- ELM and PROSITE: known motif/domain/function-site libraries for residue-level feature validation.
- ProteinGym: mutation-effect and fitness benchmark; best local target for measuring whether unknown features influence function predictions.
- CATH/AlphaFold/PDB: structural domains and residue geometry for testing whether features correspond to structure-mediated mechanisms.
- TAPE: older transfer-learning benchmark; useful for historical comparison but less central than ProteinGym for this hypothesis.

## Gaps and Opportunities

- Existing annotations are incomplete. A feature failing to align with Swiss-Prot/GO/ELM may be novel, but may also reflect missing, noisy, or coarse annotation.
- Most studies evaluate feature-concept correlation more than experimental mechanism. Stronger evidence needs held-out prediction, intervention, and structural/phylogenetic controls.
- SAE features can be dataset- and hyperparameter-sensitive. Results should be checked across seeds, sparsity settings, layers, model sizes, and sequence subsets.
- LLM-generated feature labels are useful for triage but should not be treated as evidence without quantitative validation.
- Circuit-level work is emerging. Cross-layer transcoders are promising for causal mechanisms, while ordinary per-layer SAEs are better for feature discovery.

## Recommendations for Our Experiment

- Primary dataset: use `datasets/proteingym_v1/DMS_substitutions` for mutation-effect association and `datasets/swissprot/uniprot_sprot.fasta.gz` for activation collection. Use PROSITE and ELM as known-motif controls.
- Primary model: start with ESM-2 8M or 35M for tractable activation extraction; scale to ESM-2 650M only after the pipeline works.
- Primary methods: begin with InterPLM pretrained SAEs and/or train a small Top-K SAE on Swiss-Prot shards. Add MotifAE-style local smoothness if residue-level motif coherence is weak.
- Unknown-feature workflow: identify high-activation sparse features that show low overlap with GO/Swiss-Prot/ELM/PROSITE, cluster their top proteins/residues, test enrichment for ProteinGym fitness sensitivity, then inspect structure/localization.
- Baselines: compare against raw neurons, PCA/linear probes, CAV motif detectors, integrated gradients, and ESM masked-marginal scores.
- Metrics: report concept AUROC/F1 for known annotations, Spearman correlation with ProteinGym residue importance or DMS scores, and causal steering/intervention effects where feasible.
- Validation standard: a candidate "novel mechanism" should pass at least three filters: coherent sequence or structural pattern, predictive association with function/fitness, and absence from existing annotation databases after explicit checking.
