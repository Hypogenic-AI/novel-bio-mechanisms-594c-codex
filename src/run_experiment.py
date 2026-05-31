"""Run PLM latent-feature experiments for mechanism discovery.

The pipeline uses real ESM-2-8M activations and InterPLM SAE weights. It
compares SAE features to raw hidden dimensions on ProteinGym DMS assays, screens
top candidates against ELM/PROSITE motif scans, and ablates selected features
inside ESM-2 masked-marginal scoring.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import platform
import random
import re
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from Bio import SeqIO
from huggingface_hub import hf_hub_download
from scipy import stats
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from transformers import AutoTokenizer, EsmForMaskedLM


WORKSPACE = Path(__file__).resolve().parents[1]
INTERPLM_PATH = WORKSPACE / "code" / "interPLM"
if str(INTERPLM_PATH) not in sys.path:
    sys.path.insert(0, str(INTERPLM_PATH))

from interplm.sae.dictionary import ReLUSAE  # noqa: E402


AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")


@dataclass
class Config:
    seed: int = 42
    model_name: str = "facebook/esm2_t6_8M_UR50D"
    sae_repo: str = "Elana/InterPLM-esm2-8m"
    sae_layer: int = 4
    max_assays: int = 16
    max_target_length: int = 80
    min_single_mutants: int = 500
    min_position_coverage: float = 0.55
    mask_batch_size: int = 128
    swissprot_sample_size: int = 250
    max_swissprot_length: int = 256
    min_swissprot_length: int = 20
    shuffled_repeats: int = 100
    permutation_repeats: int = 200
    top_candidate_pairs: int = 24
    intervention_candidates: int = 6
    intervention_strength: float = 1.0
    device: str = "cuda" if torch.cuda.is_available() else "cpu"


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.benchmark = True


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)


def load_proteingym_substitutions() -> pd.DataFrame:
    files = sorted((WORKSPACE / "datasets/proteingym_v1/DMS_substitutions").glob("*.parquet"))
    cols = ["DMS_id", "DMS_score", "DMS_score_bin", "mutated_sequence", "target_seq", "mutant"]
    frames = [pd.read_parquet(path, columns=cols) for path in files]
    return pd.concat(frames, ignore_index=True)


_MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def add_single_mutation_columns(df: pd.DataFrame) -> pd.DataFrame:
    parsed = df["mutant"].astype(str).str.extract(_MUT_RE)
    out = df.copy()
    out["wt_aa"] = parsed[0]
    out["position"] = pd.to_numeric(parsed[1], errors="coerce")
    out["mut_aa"] = parsed[2]
    out = out.dropna(subset=["wt_aa", "position", "mut_aa", "DMS_score"])
    out["position"] = out["position"].astype(int)
    seqs = out["target_seq"].astype(str)
    valid_pos = (out["position"] >= 1) & (out["position"] <= seqs.str.len())
    out = out.loc[valid_pos].copy()
    target_chars = [seq[pos - 1] for seq, pos in zip(out["target_seq"], out["position"], strict=False)]
    out["target_wt_aa"] = target_chars
    out = out[out["target_wt_aa"] == out["wt_aa"]].copy()
    out = out[out["mut_aa"].isin(AA_ALPHABET) & out["wt_aa"].isin(AA_ALPHABET)]
    return out


def select_assays(df: pd.DataFrame, cfg: Config) -> list[str]:
    g = (
        df.groupby("DMS_id")
        .agg(
            n=("mutant", "size"),
            length=("target_seq", lambda x: len(str(x.iloc[0]))),
            unique_targets=("target_seq", "nunique"),
            score_std=("DMS_score", "std"),
        )
        .reset_index()
    )
    eligible = g[
        (g["unique_targets"] == 1)
        & (g["length"] <= cfg.max_target_length)
        & (g["n"] >= cfg.min_single_mutants)
        & (g["score_std"] > 0)
    ].copy()
    eligible = eligible.sort_values(["length", "n", "DMS_id"])
    selected: list[str] = []
    for assay_id in eligible["DMS_id"]:
        assay = add_single_mutation_columns(df[df["DMS_id"] == assay_id])
        if assay.empty:
            continue
        length = len(str(assay["target_seq"].iloc[0]))
        covered = assay["position"].nunique() / length
        if covered < cfg.min_position_coverage:
            continue
        selected.append(assay_id)
        if len(selected) >= cfg.max_assays:
            break
    return selected


def residue_sensitivity(assay_df: pd.DataFrame) -> pd.DataFrame:
    assay = add_single_mutation_columns(assay_df)
    assay = assay.copy()
    score_std = assay["DMS_score"].std(ddof=0)
    if score_std == 0 or math.isnan(score_std):
        score_std = 1.0
    assay["DMS_score_z"] = (assay["DMS_score"] - assay["DMS_score"].mean()) / score_std
    agg = (
        assay.groupby("position")
        .agg(
            wt_aa=("wt_aa", "first"),
            n_mutants=("mutant", "size"),
            mean_score=("DMS_score", "mean"),
            mean_score_z=("DMS_score_z", "mean"),
            mean_bin=("DMS_score_bin", "mean"),
        )
        .reset_index()
    )
    agg["dms_sensitivity"] = -agg["mean_score_z"]
    if agg["mean_bin"].notna().any():
        agg["binary_sensitivity"] = 1.0 - agg["mean_bin"]
    else:
        agg["binary_sensitivity"] = np.nan
    return agg


def load_models(cfg: Config):
    tokenizer = AutoTokenizer.from_pretrained(cfg.model_name)
    model = EsmForMaskedLM.from_pretrained(cfg.model_name)
    model = model.to(cfg.device)
    model.eval()

    sae_filename = f"layer_{cfg.sae_layer}/ae_unnormalized.pt"
    sae_path = hf_hub_download(repo_id=cfg.sae_repo, filename=sae_filename)
    sae = ReLUSAE.from_pretrained(sae_path, device=cfg.device)
    sae.eval()
    return tokenizer, model, sae, sae_path


def get_token_ids(tokenizer) -> dict[str, int]:
    ids = {}
    for aa in AA_ALPHABET:
        tok_id = tokenizer.convert_tokens_to_ids(aa)
        if tok_id is None or tok_id == tokenizer.unk_token_id:
            raise ValueError(f"Could not map amino acid token {aa}")
        ids[aa] = tok_id
    return ids


@torch.no_grad()
def hidden_and_features(
    sequence: str,
    tokenizer,
    model,
    sae: ReLUSAE,
    cfg: Config,
) -> tuple[np.ndarray, np.ndarray]:
    inputs = tokenizer([sequence], return_tensors="pt").to(cfg.device)
    outputs = model(**inputs, output_hidden_states=True)
    hidden = outputs.hidden_states[cfg.sae_layer][0, 1 : len(sequence) + 1].detach()
    features = sae.encode(hidden).detach()
    return hidden.float().cpu().numpy(), features.float().cpu().numpy()


@torch.no_grad()
def masked_marginal_scores(
    sequence: str,
    mutations: pd.DataFrame,
    tokenizer,
    model,
    cfg: Config,
    hook_handle=None,
) -> pd.DataFrame:
    token_ids = get_token_ids(tokenizer)
    base = tokenizer([sequence], return_tensors="pt").to(cfg.device)
    input_ids = base["input_ids"][0]
    attention_mask = base["attention_mask"][0]
    mask_id = tokenizer.mask_token_id
    if mask_id is None:
        raise ValueError("Tokenizer has no mask token")

    positions = sorted(mutations["position"].unique().tolist())
    log_probs_by_pos: dict[int, torch.Tensor] = {}
    for start in range(0, len(positions), cfg.mask_batch_size):
        pos_batch = positions[start : start + cfg.mask_batch_size]
        batch_ids = input_ids.repeat(len(pos_batch), 1)
        batch_mask = attention_mask.repeat(len(pos_batch), 1)
        for i, pos in enumerate(pos_batch):
            batch_ids[i, pos] = mask_id  # position is 1-indexed; token 0 is CLS
        outputs = model(input_ids=batch_ids, attention_mask=batch_mask)
        log_probs = torch.log_softmax(outputs.logits, dim=-1)
        for i, pos in enumerate(pos_batch):
            log_probs_by_pos[pos] = log_probs[i, pos].detach().float().cpu()

    rows = []
    for row in mutations.itertuples(index=False):
        wt_id = token_ids[row.wt_aa]
        mut_id = token_ids[row.mut_aa]
        lp = log_probs_by_pos[int(row.position)]
        rows.append(
            {
                "mutant": row.mutant,
                "position": int(row.position),
                "wt_aa": row.wt_aa,
                "mut_aa": row.mut_aa,
                "esm_mutation_score": float(lp[mut_id] - lp[wt_id]),
            }
        )
    return pd.DataFrame(rows)


def prepare_ranked_features(features: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rank and center feature columns once for repeated Spearman tests."""
    xr = stats.rankdata(features, axis=0, method="average").astype(np.float64)
    xr = xr - xr.mean(axis=0, keepdims=True)
    xden = np.sqrt(np.sum(xr**2, axis=0))
    return xr, xden


def spearman_from_ranked(xr: np.ndarray, xden: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Spearman correlations from pre-ranked feature columns."""
    if len(y) < 4:
        return np.full(xr.shape[1], np.nan)
    yr = stats.rankdata(y, method="average").astype(np.float64)
    yr = yr - yr.mean()
    yden = np.sqrt(np.sum(yr**2))
    den = xden * yden
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = (xr.T @ yr) / den
    corr[~np.isfinite(corr)] = 0.0
    return corr


def vectorized_spearman(features: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Spearman correlations between each column in features and y."""
    valid = np.isfinite(y)
    features = features[valid]
    y = y[valid]
    xr, xden = prepare_ranked_features(features)
    return spearman_from_ranked(xr, xden, y)


def permutation_p_values_from_ranked(
    xr: np.ndarray,
    xden: np.ndarray,
    y: np.ndarray,
    observed_corr: np.ndarray,
    repeats: int,
    rng: np.random.Generator,
) -> np.ndarray:
    counts = np.zeros_like(observed_corr, dtype=np.int32)
    abs_obs = np.abs(observed_corr)
    for _ in range(repeats):
        y_perm = rng.permutation(y)
        perm_corr = np.abs(spearman_from_ranked(xr, xden, y_perm))
        counts += perm_corr >= abs_obs
    return (counts + 1) / (repeats + 1)


def shuffled_max_abs_corr_from_ranked(
    xr: np.ndarray,
    xden: np.ndarray,
    y: np.ndarray,
    repeats: int,
    rng: np.random.Generator,
) -> tuple[float, float, float]:
    vals = []
    for _ in range(repeats):
        shuffled_xr = xr[rng.permutation(xr.shape[0])]
        vals.append(float(np.nanmax(np.abs(spearman_from_ranked(shuffled_xr, xden, y)))))
    arr = np.array(vals)
    return float(arr.mean()), float(np.percentile(arr, 95)), float((arr >= arr[0]).mean())


def parse_elm_patterns() -> list[tuple[str, re.Pattern]]:
    path = WORKSPACE / "datasets/elm/elms_index.tsv"
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t", comment="#")
    patterns = []
    for row in df.itertuples(index=False):
        ident = getattr(row, "ELMIdentifier")
        regex = getattr(row, "Regex")
        if not isinstance(regex, str) or not regex:
            continue
        try:
            patterns.append((f"ELM:{ident}", re.compile(f"(?=({regex}))")))
        except re.error:
            continue
    return patterns


def prosite_to_regex(pattern: str) -> str | None:
    pat = pattern.strip().rstrip(".").replace(" ", "")
    if not pat or "/" in pat:
        return None
    pat = pat.replace("-", "")
    pat = pat.replace("<", "^").replace(">", "$")
    pat = re.sub(r"\{([A-Z]+)\}", lambda m: f"[^{m.group(1)}]", pat)
    pat = re.sub(r"[xX]\((\d+),(\d+)\)", r"[A-Z]{\1,\2}", pat)
    pat = re.sub(r"[xX]\((\d+)\)", r"[A-Z]{\1}", pat)
    pat = re.sub(r"[xX]", r"[A-Z]", pat)
    allowed = set("ABCDEFGHIJKLMNOPQRSTUVWXYZ[]^${},0123456789")
    if any(ch not in allowed for ch in pat):
        return None
    return pat


def parse_prosite_patterns() -> list[tuple[str, re.Pattern]]:
    path = WORKSPACE / "datasets/prosite/prosite.dat"
    if not path.exists():
        return []
    patterns: list[tuple[str, re.Pattern]] = []
    accession = None
    pa_lines: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("AC"):
            accession = line.split()[1].rstrip(";")
        elif line.startswith("PA"):
            pa_lines.append(line[5:].strip())
        elif line.startswith("//"):
            if accession and pa_lines:
                raw = "".join(pa_lines)
                converted = prosite_to_regex(raw)
                if converted:
                    try:
                        patterns.append((f"PROSITE:{accession}", re.compile(f"(?=({converted}))")))
                    except re.error:
                        pass
            accession = None
            pa_lines = []
    return patterns


def motif_mask(sequence: str, patterns: list[tuple[str, re.Pattern]]) -> tuple[np.ndarray, list[str]]:
    mask = np.zeros(len(sequence), dtype=bool)
    hits: list[str] = []
    for name, pattern in patterns:
        try:
            for match in pattern.finditer(sequence):
                matched = match.group(1)
                if matched is None:
                    continue
                start = match.start(1)
                end = start + len(matched)
                if end > start:
                    mask[start:end] = True
                    hits.append(name)
        except Exception:
            continue
    return mask, sorted(set(hits))


def feature_motif_overlap(values: np.ndarray, mask: np.ndarray) -> dict:
    active = values > np.quantile(values, 0.85)
    if active.sum() == 0:
        active[np.argmax(values)] = True
    overlap = active & mask
    table = np.array(
        [
            [int(overlap.sum()), int((active & ~mask).sum())],
            [int((~active & mask).sum()), int((~active & ~mask).sum())],
        ]
    )
    if table.min() >= 0 and table.sum() > 0 and active.any() and mask.any():
        _, p = stats.fisher_exact(table, alternative="greater")
    else:
        p = np.nan
    return {
        "active_positions": int(active.sum()),
        "motif_positions": int(mask.sum()),
        "overlap_positions": int(overlap.sum()),
        "overlap_fraction_active": float(overlap.sum() / max(active.sum(), 1)),
        "fisher_p": float(p) if np.isfinite(p) else np.nan,
    }


def load_swissprot_sample(cfg: Config, rng: np.random.Generator) -> list[tuple[str, str]]:
    fasta = WORKSPACE / "datasets/swissprot/uniprot_sprot.fasta.gz"
    candidates: list[tuple[str, str]] = []
    with gzip.open(fasta, "rt") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            seq = str(rec.seq).replace("U", "X").replace("O", "X")
            if (
                cfg.min_swissprot_length <= len(seq) <= cfg.max_swissprot_length
                and set(seq).issubset(set(AA_ALPHABET))
            ):
                candidates.append((rec.id, seq))
            if len(candidates) >= cfg.swissprot_sample_size * 20:
                break
    idx = rng.choice(len(candidates), size=min(cfg.swissprot_sample_size, len(candidates)), replace=False)
    return [candidates[i] for i in idx]


@torch.no_grad()
def candidate_swissprot_alignment(
    candidate_features: list[int],
    tokenizer,
    model,
    sae: ReLUSAE,
    patterns: list[tuple[str, re.Pattern]],
    cfg: Config,
    rng: np.random.Generator,
) -> pd.DataFrame:
    if not candidate_features:
        return pd.DataFrame()
    records = load_swissprot_sample(cfg, rng)
    feature_values = {feat: [] for feat in candidate_features}
    motif_labels: list[bool] = []
    for _, seq in tqdm(records, desc="Swiss-Prot candidate motif screen"):
        mask, _ = motif_mask(seq, patterns)
        inputs = tokenizer([seq], return_tensors="pt").to(cfg.device)
        outputs = model(**inputs, output_hidden_states=True)
        hidden = outputs.hidden_states[cfg.sae_layer][0, 1 : len(seq) + 1]
        feats = sae.encode_feat_subset(hidden, candidate_features).detach().float().cpu().numpy()
        for j, feat in enumerate(candidate_features):
            feature_values[feat].extend(feats[:, j].tolist())
        motif_labels.extend(mask.tolist())
    y = np.array(motif_labels, dtype=bool)
    rows = []
    for feat in candidate_features:
        vals = np.array(feature_values[feat], dtype=float)
        if y.any() and (~y).any() and np.unique(vals).size > 1:
            auc = roc_auc_score(y.astype(int), vals)
            rho, p = stats.spearmanr(vals, y.astype(float))
        else:
            auc = np.nan
            rho, p = np.nan, np.nan
        rows.append(
            {
                "feature_id": feat,
                "swissprot_tokens": int(len(y)),
                "swissprot_motif_positive_fraction": float(y.mean()) if len(y) else np.nan,
                "motif_auc": float(auc) if np.isfinite(auc) else np.nan,
                "motif_spearman": float(rho) if np.isfinite(rho) else np.nan,
                "motif_spearman_p": float(p) if np.isfinite(p) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def make_ablation_hook(sae: ReLUSAE, feature_id: int, strength: float):
    decoder_vec = sae.decoder.weight[:, feature_id].detach()

    def hook(_module, _inputs, output):
        hidden = output[0] if isinstance(output, tuple) else output
        flat = hidden.reshape(-1, hidden.shape[-1])
        feat = sae.encode_feat_subset(flat, [feature_id]).reshape(*hidden.shape[:-1])
        patched = hidden - strength * feat.unsqueeze(-1) * decoder_vec
        if isinstance(output, tuple):
            return (patched,) + output[1:]
        return patched

    return hook


def run_interventions(
    candidates: pd.DataFrame,
    assay_frames: dict[str, pd.DataFrame],
    tokenizer,
    model,
    sae: ReLUSAE,
    cfg: Config,
) -> pd.DataFrame:
    rows = []
    layer_module = model.esm.encoder.layer[cfg.sae_layer - 1].output
    selected = candidates.head(cfg.intervention_candidates)
    for cand in tqdm(selected.itertuples(index=False), total=len(selected), desc="Feature ablations"):
        assay_df = add_single_mutation_columns(assay_frames[cand.DMS_id])
        sequence = str(assay_df["target_seq"].iloc[0])
        _, features = hidden_and_features(sequence, tokenizer, model, sae, cfg)
        feat_vals = features[:, int(cand.feature_id)]
        pos_df = residue_sensitivity(assay_frames[cand.DMS_id])
        mutations = assay_df[["mutant", "position", "wt_aa", "mut_aa"]].drop_duplicates()
        baseline = masked_marginal_scores(sequence, mutations, tokenizer, model, cfg)
        handle = layer_module.register_forward_hook(
            make_ablation_hook(sae, int(cand.feature_id), cfg.intervention_strength)
        )
        try:
            patched = masked_marginal_scores(sequence, mutations, tokenizer, model, cfg)
        finally:
            handle.remove()
        merged = baseline.merge(
            patched,
            on=["mutant", "position", "wt_aa", "mut_aa"],
            suffixes=("_baseline", "_ablated"),
        )
        merged["abs_score_change"] = (
            merged["esm_mutation_score_ablated"] - merged["esm_mutation_score_baseline"]
        ).abs()
        pos_change = (
            merged.groupby("position")["abs_score_change"].mean().reset_index()
            .merge(pos_df[["position", "dms_sensitivity"]], on="position", how="left")
        )
        q75 = np.quantile(feat_vals, 0.75)
        q25 = np.quantile(feat_vals, 0.25)
        high_pos = {i + 1 for i, v in enumerate(feat_vals) if v >= q75}
        low_pos = {i + 1 for i, v in enumerate(feat_vals) if v <= q25}
        high = pos_change[pos_change["position"].isin(high_pos)]["abs_score_change"].values
        low = pos_change[pos_change["position"].isin(low_pos)]["abs_score_change"].values
        if len(high) > 0 and len(low) > 0:
            u_p = stats.mannwhitneyu(high, low, alternative="greater").pvalue
            effect = float(np.mean(high) - np.mean(low))
        else:
            u_p = np.nan
            effect = np.nan
        rows.append(
            {
                "DMS_id": cand.DMS_id,
                "feature_id": int(cand.feature_id),
                "n_mutations": int(len(merged)),
                "n_positions": int(pos_change["position"].nunique()),
                "mean_abs_change_high_activation": float(np.mean(high)) if len(high) else np.nan,
                "mean_abs_change_low_activation": float(np.mean(low)) if len(low) else np.nan,
                "intervention_effect_high_minus_low": effect,
                "mannwhitney_p": float(u_p) if np.isfinite(u_p) else np.nan,
                "strength": cfg.intervention_strength,
            }
        )
    return pd.DataFrame(rows)


def analyze_assay(
    assay_id: str,
    assay_df: pd.DataFrame,
    tokenizer,
    model,
    sae: ReLUSAE,
    patterns: list[tuple[str, re.Pattern]],
    cfg: Config,
    rng: np.random.Generator,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    assay = add_single_mutation_columns(assay_df)
    sequence = str(assay["target_seq"].iloc[0])
    pos_df = residue_sensitivity(assay_df)
    mutations = assay[["mutant", "position", "wt_aa", "mut_aa"]].drop_duplicates()
    esm_scores = masked_marginal_scores(sequence, mutations, tokenizer, model, cfg)
    esm_pos = esm_scores.groupby("position")["esm_mutation_score"].mean().reset_index()
    esm_pos["esm_sensitivity"] = -esm_pos["esm_mutation_score"]
    pos_df = pos_df.merge(esm_pos[["position", "esm_sensitivity"]], on="position", how="left")
    hidden, features = hidden_and_features(sequence, tokenizer, model, sae, cfg)

    positions = pos_df["position"].to_numpy(dtype=int) - 1
    y = pos_df["dms_sensitivity"].to_numpy(dtype=float)
    y_esm = pos_df["esm_sensitivity"].to_numpy(dtype=float)
    feat_pos = features[positions]
    hidden_pos = hidden[positions]

    sae_xr, sae_xden = prepare_ranked_features(feat_pos)
    raw_xr, raw_xden = prepare_ranked_features(hidden_pos)
    sae_corr = spearman_from_ranked(sae_xr, sae_xden, y)
    raw_corr = spearman_from_ranked(raw_xr, raw_xden, y)
    sae_esm_corr = spearman_from_ranked(sae_xr, sae_xden, y_esm)
    raw_esm_corr = spearman_from_ranked(raw_xr, raw_xden, y_esm)
    pvals = permutation_p_values_from_ranked(sae_xr, sae_xden, y, sae_corr, cfg.permutation_repeats, rng)
    qvals = multipletests(pvals, method="fdr_bh")[1]
    shuffle_mean, shuffle_p95, _ = shuffled_max_abs_corr_from_ranked(
        sae_xr, sae_xden, y, cfg.shuffled_repeats, rng
    )

    motif_positions, motif_hits = motif_mask(sequence, patterns)
    top_idx = np.argsort(-np.abs(sae_corr))[: min(20, len(sae_corr))]
    candidate_rows = []
    for rank, feat_id in enumerate(top_idx, start=1):
        overlap = feature_motif_overlap(features[:, feat_id], motif_positions)
        candidate_rows.append(
            {
                "DMS_id": assay_id,
                "feature_id": int(feat_id),
                "rank_within_assay": rank,
                "spearman_dms": float(sae_corr[feat_id]),
                "abs_spearman_dms": float(abs(sae_corr[feat_id])),
                "spearman_esm": float(sae_esm_corr[feat_id]),
                "permutation_p": float(pvals[feat_id]),
                "bh_q": float(qvals[feat_id]),
                "target_length": len(sequence),
                "assay_positions": int(len(pos_df)),
                "motif_hit_count": int(len(motif_hits)),
                **overlap,
            }
        )

    top_raw = int(np.nanargmax(np.abs(raw_corr)))
    top_sae = int(np.nanargmax(np.abs(sae_corr)))
    assay_summary = {
        "DMS_id": assay_id,
        "target_length": len(sequence),
        "n_single_mutations": int(len(assay)),
        "n_positions": int(pos_df["position"].nunique()),
        "dms_esm_spearman": float(stats.spearmanr(y, y_esm, nan_policy="omit").statistic),
        "best_sae_feature": top_sae,
        "best_sae_abs_spearman": float(abs(sae_corr[top_sae])),
        "best_sae_spearman": float(sae_corr[top_sae]),
        "best_sae_esm_spearman": float(sae_esm_corr[top_sae]),
        "best_sae_perm_p": float(pvals[top_sae]),
        "best_sae_bh_q": float(qvals[top_sae]),
        "best_raw_dim": top_raw,
        "best_raw_abs_spearman": float(abs(raw_corr[top_raw])),
        "best_raw_spearman": float(raw_corr[top_raw]),
        "best_raw_esm_abs_spearman": float(np.nanmax(np.abs(raw_esm_corr))),
        "shuffled_sae_best_abs_spearman_mean": shuffle_mean,
        "shuffled_sae_best_abs_spearman_p95": shuffle_p95,
        "motif_positions": int(motif_positions.sum()),
        "motif_hit_count": int(len(motif_hits)),
        "motif_hits_preview": ";".join(motif_hits[:8]),
    }

    pos_out = pos_df.copy()
    pos_out["target_aa"] = [sequence[p - 1] for p in pos_out["position"]]
    pos_out["motif_position"] = [bool(motif_positions[p - 1]) for p in pos_out["position"]]
    pos_out["best_sae_feature_activation"] = features[positions, top_sae]
    pos_out["best_raw_dim_activation"] = hidden[positions, top_raw]
    pos_out["DMS_id"] = assay_id
    return assay_summary, pd.DataFrame(candidate_rows), pos_out


def make_figures(assay_summary: pd.DataFrame, candidates: pd.DataFrame, interventions: pd.DataFrame) -> None:
    sns.set_theme(style="whitegrid")
    fig_dir = WORKSPACE / "figures"
    fig_dir.mkdir(exist_ok=True)

    plt.figure(figsize=(9, 5))
    plot_df = assay_summary.melt(
        id_vars=["DMS_id"],
        value_vars=[
            "best_sae_abs_spearman",
            "best_raw_abs_spearman",
            "shuffled_sae_best_abs_spearman_p95",
        ],
        var_name="condition",
        value_name="abs_spearman",
    )
    sns.boxplot(data=plot_df, x="condition", y="abs_spearman", color="#d9e2ec")
    sns.stripplot(data=plot_df, x="condition", y="abs_spearman", color="#334e68", alpha=0.75)
    plt.xticks(rotation=15, ha="right")
    plt.xlabel("")
    plt.ylabel("Best absolute Spearman with DMS sensitivity")
    plt.tight_layout()
    plt.savefig(fig_dir / "baseline_comparison.png", dpi=180)
    plt.close()

    top = candidates.sort_values("abs_spearman_dms", ascending=False).head(20).copy()
    top["feature_label"] = top["DMS_id"].str.split("_").str[:2].str.join("_") + " / F" + top["feature_id"].astype(str)
    plt.figure(figsize=(9, 7))
    sns.barplot(data=top, y="feature_label", x="abs_spearman_dms", hue="overlap_fraction_active", dodge=False)
    plt.xlabel("Absolute Spearman with DMS sensitivity")
    plt.ylabel("Assay / SAE feature")
    plt.legend(title="Motif overlap\nfraction", loc="lower right")
    plt.tight_layout()
    plt.savefig(fig_dir / "top_candidates.png", dpi=180)
    plt.close()

    if not interventions.empty:
        plt.figure(figsize=(8, 4))
        inter = interventions.copy()
        inter["label"] = inter["DMS_id"].str.split("_").str[:2].str.join("_") + " / F" + inter["feature_id"].astype(str)
        sns.barplot(data=inter, x="intervention_effect_high_minus_low", y="label", color="#2f855a")
        plt.axvline(0, color="black", linewidth=1)
        plt.xlabel("Mean |masked-score change|: high activation minus low activation")
        plt.ylabel("Ablated candidate")
        plt.tight_layout()
        plt.savefig(fig_dir / "intervention_effects.png", dpi=180)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-assays", type=int, default=None)
    parser.add_argument("--quick", action="store_true", help="Smaller run for validation.")
    args = parser.parse_args()

    cfg = Config()
    if args.max_assays is not None:
        cfg.max_assays = args.max_assays
    if args.quick:
        cfg.max_assays = min(cfg.max_assays, 4)
        cfg.swissprot_sample_size = 60
        cfg.shuffled_repeats = 30
        cfg.permutation_repeats = 50
        cfg.intervention_candidates = 2
    set_seed(cfg.seed)
    rng = np.random.default_rng(cfg.seed)
    start_time = time.time()

    results_dir = WORKSPACE / "results"
    results_dir.mkdir(exist_ok=True)
    write_json(results_dir / "config.json", asdict(cfg))
    hardware = {
        "python": sys.version,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "cuda_available": torch.cuda.is_available(),
        "cuda_version": torch.version.cuda,
        "gpu_count": torch.cuda.device_count(),
        "gpus": [
            {
                "name": torch.cuda.get_device_name(i),
                "memory_total_mib": torch.cuda.get_device_properties(i).total_memory // 2**20,
            }
            for i in range(torch.cuda.device_count())
        ],
    }
    write_json(results_dir / "environment.json", hardware)

    df = load_proteingym_substitutions()
    selected = select_assays(df, cfg)
    selected_df = pd.DataFrame({"DMS_id": selected})
    selected_df.to_csv(results_dir / "selected_assays.csv", index=False)
    assay_frames = {assay_id: df[df["DMS_id"] == assay_id].copy() for assay_id in selected}

    patterns = parse_elm_patterns() + parse_prosite_patterns()
    write_json(
        results_dir / "motif_pattern_summary.json",
        {"compiled_patterns": len(patterns), "elm_patterns": len(parse_elm_patterns()), "prosite_patterns": len(parse_prosite_patterns())},
    )

    tokenizer, model, sae, sae_path = load_models(cfg)
    write_json(results_dir / "model_paths.json", {"sae_path": str(sae_path), "model_name": cfg.model_name})

    summaries = []
    candidate_frames = []
    pos_frames = []
    for assay_id in tqdm(selected, desc="Assays"):
        summary, cand, pos = analyze_assay(
            assay_id, assay_frames[assay_id], tokenizer, model, sae, patterns, cfg, rng
        )
        summaries.append(summary)
        candidate_frames.append(cand)
        pos_frames.append(pos)

    assay_summary = pd.DataFrame(summaries)
    candidates = pd.concat(candidate_frames, ignore_index=True)
    positions = pd.concat(pos_frames, ignore_index=True)

    # Candidate novelty heuristic: strong DMS association and little direct motif overlap.
    candidates = candidates.sort_values("abs_spearman_dms", ascending=False).reset_index(drop=True)
    candidates["low_target_motif_alignment"] = (
        (candidates["overlap_fraction_active"] <= 0.10)
        & ((candidates["fisher_p"].isna()) | (candidates["fisher_p"] > 0.05))
    )
    low_align_features = (
        candidates[candidates["low_target_motif_alignment"]]
        .head(cfg.top_candidate_pairs)["feature_id"]
        .drop_duplicates()
        .astype(int)
        .tolist()
    )
    swiss_align = candidate_swissprot_alignment(
        low_align_features, tokenizer, model, sae, patterns, cfg, rng
    )
    if not swiss_align.empty:
        candidates = candidates.merge(swiss_align, on="feature_id", how="left")
        candidates["low_swissprot_motif_alignment"] = candidates["motif_auc"].isna() | (candidates["motif_auc"] <= 0.60)
    else:
        candidates["low_swissprot_motif_alignment"] = np.nan

    top_for_intervention = candidates[
        candidates["low_target_motif_alignment"]
        & (candidates["low_swissprot_motif_alignment"].fillna(True))
    ].sort_values("abs_spearman_dms", ascending=False)
    interventions = run_interventions(top_for_intervention, assay_frames, tokenizer, model, sae, cfg)

    assay_summary.to_csv(results_dir / "assay_summary.csv", index=False)
    candidates.to_csv(results_dir / "candidate_features.csv", index=False)
    positions.to_csv(results_dir / "position_level_results.csv", index=False)
    swiss_align.to_csv(results_dir / "swissprot_candidate_motif_alignment.csv", index=False)
    interventions.to_csv(results_dir / "intervention_results.csv", index=False)

    # Aggregate statistical tests.
    paired = assay_summary.dropna(subset=["best_sae_abs_spearman", "best_raw_abs_spearman"])
    if len(paired) >= 3:
        sae_vs_raw = stats.wilcoxon(
            paired["best_sae_abs_spearman"], paired["best_raw_abs_spearman"], alternative="greater"
        )
        sae_vs_shuffle = stats.wilcoxon(
            paired["best_sae_abs_spearman"],
            paired["shuffled_sae_best_abs_spearman_p95"],
            alternative="greater",
        )
    else:
        sae_vs_raw = sae_vs_shuffle = None
    aggregate = {
        "n_assays": int(len(assay_summary)),
        "median_best_sae_abs_spearman": float(assay_summary["best_sae_abs_spearman"].median()),
        "median_best_raw_abs_spearman": float(assay_summary["best_raw_abs_spearman"].median()),
        "median_shuffled_p95_abs_spearman": float(assay_summary["shuffled_sae_best_abs_spearman_p95"].median()),
        "median_dms_esm_spearman": float(assay_summary["dms_esm_spearman"].median()),
        "n_low_target_alignment_candidates": int(candidates["low_target_motif_alignment"].sum()),
        "n_interventions": int(len(interventions)),
        "sae_vs_raw_wilcoxon_stat": float(sae_vs_raw.statistic) if sae_vs_raw else None,
        "sae_vs_raw_wilcoxon_p": float(sae_vs_raw.pvalue) if sae_vs_raw else None,
        "sae_vs_shuffle95_wilcoxon_stat": float(sae_vs_shuffle.statistic) if sae_vs_shuffle else None,
        "sae_vs_shuffle95_wilcoxon_p": float(sae_vs_shuffle.pvalue) if sae_vs_shuffle else None,
        "runtime_seconds": float(time.time() - start_time),
    }
    write_json(results_dir / "aggregate_stats.json", aggregate)
    make_figures(assay_summary, candidates, interventions)

    print(json.dumps(aggregate, indent=2))


if __name__ == "__main__":
    main()
