# Cloned Repositories

Total repositories cloned: 12.

## Primary Tooling

### interPLM

- URL: https://github.com/ElanaPearl/interPLM
- Location: `code/interPLM/`
- Purpose: SAE training, evaluation, feature activation collection, annotation association, and dashboard generation for PLM features.
- Key files: `scripts/extract_embeddings.py`, `scripts/evaluate_sae.py`, `scripts/embed_annotations.py`, `examples/train_basic_sae.py`
- Notes: Provides pretrained ESM-2 SAEs on Hugging Face and a walkthrough using Swiss-Prot FASTA plus UniProtKB annotation TSVs.

### InterProt

- URL: https://github.com/etowahadams/interprot
- Location: `code/interprot/`
- Purpose: SAE inference, visualization-file generation, feature analysis, steering, and linear probes over SAE latents.
- Key files: `notebooks/sae_inference.ipynb`, `interprot/make_viz_files/`, `interprot/scripts/`
- Notes: Most useful for hypothesis generation and visual inspection of candidate latent features.

### Rep_SAEs_PLMs

- URL: https://github.com/onkarsg10/Rep_SAEs_PLMs
- Location: `code/Rep_SAEs_PLMs/`
- Purpose: Training and interpreting SAEs/transcoders on ESM2 protein-level and amino-acid-level representations.
- Key files: `Folder_Random_Seed_Regular_Pooling_Scripts/`, `Folder_Random_Seed_TC_Pooling_Scripts/`, `goe_analysis/`, `Folder_Pooling_Autointerp/`
- Notes: Requires Swiss-Prot TSV, UniRef50 FASTA, GO OBO, GOA human GAF, checkpoints, GPU, and optional Claude API key for automated interpretation.

### ProteinGym

- URL: https://github.com/OATML-Markslab/ProteinGym
- Location: `code/ProteinGym/`
- Purpose: Benchmark data handling, scoring scripts, and performance aggregation for mutation-effect prediction.
- Key files: `scripts/scoring_DMS_zero_shot/`, `benchmarks/`, `proteingym/`
- Notes: Local data snapshot is in `datasets/proteingym_v1/`; full upstream score archives are much larger and were not downloaded.

### ESM

- URL: https://github.com/facebookresearch/esm
- Location: `code/esm/`
- Purpose: Official ESM/ESM-2/ESMFold implementation, embedding extraction, zero-shot variant prediction, inverse folding.
- Key files: `esm/`, `scripts/`, `examples/variant-prediction/`
- Notes: Useful backbone for generating PLM activations. ESMFold dependencies may require older Python/PyTorch for full structure prediction.

## Comparators and Baselines

### plm-sae

- URL: https://github.com/edithvillegas/plm-sae
- Location: `code/plm-sae/`
- Purpose: Lightweight SAE inference code for ESM2-8M layer 3 and sequence steering experiments.
- Notes: README points to Hugging Face weights `evillegasgarcia/sae_esm2_6_l3`.

### reticular-sae

- URL: https://github.com/johnyang101/reticular-sae
- Location: `code/reticular-sae/`
- Purpose: SAEFold/Matryoshka SAE training and evaluation for ESM2-3B/ESMFold interpretability.
- Notes: Cloned with `GIT_LFS_SKIP_SMUDGE=1` because pretrained checkpoint download hit the repo's Git LFS quota. Code and pointer files are present; checkpoint blobs are not local.

### provis

- URL: https://github.com/salesforce/provis
- Location: `code/provis/`
- Purpose: Attention visualization and probing from "BERTology Meets Biology".
- Key files: `notebooks/provis.ipynb`, `protein_attention/attention_analysis/`, `protein_attention/probing/`
- Notes: Useful baseline for attention-based interpretability, though dependencies are older.

### xai-proteins

- URL: https://github.com/markuswenzel/xai-proteins
- Location: `code/xai-proteins/`
- Purpose: Integrated gradients extended to transformer internals for GO/EC protein-function predictors.
- Key files: `go/`, `ec/`
- Notes: Useful attribution baseline for comparing SAE feature maps to supervised explanation methods.

### PLMNeuron

- URL: https://github.com/arjun-banerjee/PLMNeuron
- Location: `code/PLMNeuron/`
- Purpose: Automated neuron labeling and neuron activation-guided steering for PLMs.
- Key files: `datageneration/`, `analysis/`, `benchmarks/`
- Notes: Includes generated CSVs/images and notebooks; not deeply documented beyond a short README.

### PLM-eXplain

- URL: https://github.com/AIT4LIFE-UU/PLM-eXplain
- Location: `code/PLM-eXplain/`
- Purpose: Pointer repository for the PLM-eXplain paper.
- Notes: Only contains a README linking to the implementation repository.

### PLM-eXplain-full

- URL: https://github.com/janvaneck1994/PLM-eXplain
- Location: `code/PLM-eXplain-full/`
- Purpose: Implementation for partitioning PLM embeddings into interpretable biochemical feature subspace plus residual predictive subspace.
- Key files: `train_plmx/`, `train_downstream_tasks/`
- Notes: README points to Google Drive data; code is available but documentation is sparse.
