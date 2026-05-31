# Discovering Novel Biological Mechanisms from Protein Language Models

This workspace runs a pilot mechanistic-biology experiment using real ESM-2 internals. It tests whether InterPLM sparse autoencoder features mark ProteinGym mutation-sensitive residues that are not explained by ELM/PROSITE motif scans.

## Key Findings

- On 16 short ProteinGym assays, best SAE features correlated with DMS residue sensitivity (median abs Spearman 0.547), but did not significantly beat raw hidden dimensions or shuffled SAE controls.
- Three low-motif-alignment candidate feature-assay pairs were found.
- One candidate, feature 6419 in `YNZC_BACSU_Tsuboyama_2023_2JVD`, showed an exploratory ablation effect on ESM masked-marginal scores.
- The result is hypothesis-generating, not evidence of a confirmed novel biological mechanism.

See [REPORT.md](REPORT.md) for the full methodology, results, and limitations.

## Reproduce

```bash
source .venv/bin/activate
uv sync
python src/run_experiment.py
```

The script writes outputs to `results/` and figures to `figures/`. A quick validation run is available:

```bash
python src/run_experiment.py --quick
```

## File Structure

- `planning.md`: preregistered motivation, novelty, and experiment plan.
- `src/run_experiment.py`: full experiment pipeline.
- `results/`: CSV/JSON outputs, including aggregate stats and candidate tables.
- `figures/`: generated plots.
- `datasets/`: local ProteinGym, Swiss-Prot, ELM, and PROSITE resources.
- `code/`: cloned baseline/tooling repositories.
- `literature_review.md` and `resources.md`: gathered research context.
