# Resources Catalog

## Summary

Resources gathered for "Discovering Novel Biological Mechanisms from Protein Language Models":

- Papers downloaded: 24 main PDFs
- Dataset families downloaded: 4
- Code repositories cloned: 12
- Local environment: `.venv` managed by `uv`, dependencies in `pyproject.toml`

## Papers

| Title | Authors | Year | File | Key Info |
|---|---|---:|---|---|
| InterPLM | Simon, Zou | 2024/2025 | `papers/2412.12101_interplm_sparse_autoencoders.pdf` | Main SAE feature-discovery framework. |
| Sparse autoencoders uncover biologically interpretable features | Gujral et al. | 2025 | `papers/gujral_2025_sparse_autoencoders_plm_representations.pdf` | SAE/transcoder analysis with GO and automated interpretation. |
| From Mechanistic Interpretability to Mechanistic Biology | Adams et al. | 2025 | `papers/adams_2025_mechanistic_interpretability_to_mechanistic_biology.pdf` | Feature-to-hypothesis workflow and InterProt. |
| MotifAE | Hou et al. | 2025 | `papers/hou_2025_motifae.pdf` | Motif-specific SAE with smoothness loss. |
| Protein Circuit Tracing via Cross-layer Transcoders | Tsui et al. | 2026 | `papers/2602.12026_protein_circuit_tracing_cross_layer_transcoders.pdf` | Circuit-level method for PLM internals. |
| ProteinGym | Notin et al. | 2023 | `papers/notin_2023_proteingym.pdf` | Fitness/mutation-effect benchmark. |
| ESMFold / ESM-2 | Lin et al. | 2023 | `papers/lin_2023_esmfold_evolutionary_scale_structure_prediction.pdf` | PLM backbone and structure prediction. |
| ESM scaling | Rives et al. | 2021 | `papers/rives_2021_biological_structure_function_scaling_protein_lms.pdf` | Foundation for PLM biological representations. |
| BERTology Meets Biology | Vig et al. | 2020 | `papers/2006.15222_bertology_meets_biology.pdf` | Attention/probing baseline. |
| Inner workings of transformer models for protein function prediction | Wenzel et al. | 2024 | `papers/2309.03631_inner_workings_transformers_protein_function.pdf` | Integrated-gradient attribution baseline. |

See `papers/README.md` for the full list of 24 PDFs.

## Datasets

| Name | Source | Size / Count | Task | Location | Notes |
|---|---|---:|---|---|---|
| ProteinGym v1 | https://huggingface.co/datasets/OATML-Markslab/ProteinGym_v1 | 2,818,579 rows | Fitness and mutation-effect prediction | `datasets/proteingym_v1/` | Downloaded as parquet snapshot. |
| UniProtKB/Swiss-Prot FASTA | https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | 574,627 records | SAE training/evaluation sequences | `datasets/swissprot/` | Reviewed sequences, compressed FASTA. |
| PROSITE | https://ftp.expasy.org/databases/prosite/ | 2,731 PROSITE entries, 1,442 ProRule entries | Known motif/domain/site alignment | `datasets/prosite/` | Includes `prosite.dat` and `prorule.dat`. |
| ELM | http://elm.eu.org/ | 354 motif class rows; 101 instance rows in current export | Short linear motif controls | `datasets/elm/` | TSV exports downloaded. |

See `datasets/README.md` and `datasets/dataset_summary.json` for schemas, counts, and loading examples.

## Code Repositories

| Name | URL | Purpose | Location | Notes |
|---|---|---|---|---|
| interPLM | https://github.com/ElanaPearl/interPLM | SAE training, evaluation, annotation association, dashboarding | `code/interPLM/` | Most immediately reusable. |
| interprot | https://github.com/etowahadams/interprot | SAE inference, visualization, linear probes, steering | `code/interprot/` | Good hypothesis triage tooling. |
| Rep_SAEs_PLMs | https://github.com/onkarsg10/Rep_SAEs_PLMs | SAE/transcoder training and GO/LLM interpretation | `code/Rep_SAEs_PLMs/` | Requires large inputs and GPU. |
| ProteinGym | https://github.com/OATML-Markslab/ProteinGym | Benchmark scoring and model-comparison utilities | `code/ProteinGym/` | Pair with local HF snapshot. |
| esm | https://github.com/facebookresearch/esm | ESM/ESM-2/ESMFold implementation | `code/esm/` | Activation extraction backbone. |
| plm-sae | https://github.com/edithvillegas/plm-sae | Lightweight ESM2-8M SAE inference | `code/plm-sae/` | Hugging Face SAE weights. |
| reticular-sae | https://github.com/johnyang101/reticular-sae | ESM2-3B/ESMFold SAEFold code | `code/reticular-sae/` | LFS checkpoint blobs skipped due quota. |
| provis | https://github.com/salesforce/provis | Attention visualization/probing | `code/provis/` | Older dependencies. |
| xai-proteins | https://github.com/markuswenzel/xai-proteins | Integrated-gradient transformer attribution | `code/xai-proteins/` | GO/EC explainability baseline. |
| PLMNeuron | https://github.com/arjun-banerjee/PLMNeuron | Neuron labeling and steering | `code/PLMNeuron/` | Contains notebooks, CSVs, images. |
| PLM-eXplain | https://github.com/AIT4LIFE-UU/PLM-eXplain | Pointer repo | `code/PLM-eXplain/` | Links to implementation. |
| PLM-eXplain-full | https://github.com/janvaneck1994/PLM-eXplain | Embedding partitioning implementation | `code/PLM-eXplain-full/` | Sparse README, code present. |

See `code/README.md` for entry points and requirements.

## Search Strategy

The search began with the local paper-finder service using the query "protein language model interpretability latent features circuits protein function motifs". It returned 151 papers, including 26 high-relevance and 28 medium-relevance entries. I then performed targeted arXiv, PubMed/PMC, Nature, Hugging Face, GitHub, UniProt, Expasy, and ELM searches for:

- protein language model sparse autoencoders
- mechanistic interpretability and circuits in PLMs
- motif discovery from ESM2/PLM embeddings
- protein function and mutation-effect benchmarks
- Swiss-Prot, PROSITE, ELM, and ProteinGym data sources

Selection prioritized direct fit to the hypothesis, accessible PDFs, available code/data, and benchmark usefulness.

## Challenges Encountered

- PNAS and PMC direct PDF endpoints used bot-protection or returned placeholder HTML. The final PDFs were downloaded through the PMC Open Access S3 bucket where available.
- `datasets.load_dataset()` could not infer ProteinGym's repository layout, so the Hugging Face repository was downloaded with `huggingface_hub.snapshot_download()`.
- `reticular-sae` attempted to download 700 MB Git LFS checkpoint files and hit an LFS quota error. It was recloned with `GIT_LFS_SKIP_SMUDGE=1`, so code is present but checkpoint blobs are not.
- Some relevant papers, especially MotifAE and circuit-tracing work, have limited or no official code available.

## Recommendations for Experiment Design

1. Primary implementation path: use `code/interPLM/` with `datasets/swissprot/uniprot_sprot.fasta.gz` to extract ESM-2 activations and run SAE feature analysis.
2. Primary evaluation data: use `datasets/proteingym_v1/DMS_substitutions` for mutation-effect association and residue/function sensitivity.
3. Known-annotation controls: use PROSITE and ELM for motif/domain/site overlap; add UniProtKB TSV annotation export if larger supervised concept alignment is needed.
4. Baselines: raw neurons, InterPLM SAE features, CAV motif localization, integrated gradients from `xai-proteins`, and ESM zero-shot ProteinGym scoring.
5. Novel-feature criterion: prioritize features with coherent top-activating proteins/residues, low overlap with known annotations, and measurable association with ProteinGym fitness or clinical labels.

## Research Execution Outputs

The automated research run used the resources above to execute a real ESM-2/InterPLM SAE pilot experiment. Outputs are:

- `planning.md`: motivation, novelty assessment, and preregistered experimental plan.
- `src/run_experiment.py`: reproducible pipeline for ProteinGym assay selection, ESM masked-marginal scoring, SAE feature extraction, motif screening, shuffled/raw baselines, and feature ablation.
- `results/aggregate_stats.json`: aggregate statistics for 16 ProteinGym assays.
- `results/candidate_features.csv`: ranked SAE feature candidates and motif-alignment metrics.
- `results/intervention_results.csv`: ablation effects for low-motif-alignment candidates.
- `figures/`: baseline comparison, top candidate, and intervention plots.
- `REPORT.md`: final research report with actual results and limitations.
