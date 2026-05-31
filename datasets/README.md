# Downloaded Datasets

Data files are local for the experiment runner but excluded from git by `datasets/.gitignore`.

## ProteinGym v1

- Source: https://huggingface.co/datasets/OATML-Markslab/ProteinGym_v1
- Location: `datasets/proteingym_v1/`
- Size: 142 MB repository download; 2,818,579 total local rows across four parquet configurations
- Task: protein fitness prediction, mutation-effect prediction, clinical variant classification
- Splits/configs:
  - `DMS_substitutions`: 2,465,767 rows
  - `DMS_indels`: 287,207 rows
  - `clinical_substitutions`: 62,727 rows
  - `clinical_indels`: 2,878 rows
- Columns:
  - DMS: `DMS_score`, `DMS_score_bin`, `mutated_sequence`, `target_seq`, `mutant`, `DMS_id`
  - Clinical: `mutated_sequence`, `target_seq`, `mutant`, `protein_id`, `annotation`
- Samples: `datasets/proteingym_v1/*/samples/sample_records.json`

Download:

```python
from huggingface_hub import snapshot_download

snapshot_download(
    repo_id="OATML-Markslab/ProteinGym_v1",
    repo_type="dataset",
    local_dir="datasets/proteingym_v1",
)
```

Load:

```python
import pandas as pd

df = pd.read_parquet("datasets/proteingym_v1/DMS_substitutions/train-00000-of-00005.parquet")
```

## UniProtKB/Swiss-Prot FASTA

- Source: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
- Location: `datasets/swissprot/uniprot_sprot.fasta.gz`
- Size: 93,457,057 bytes compressed
- Records: 574,627 reviewed protein sequences
- Task: activation collection and SAE training/evaluation sequence source; matches InterPLM walkthrough guidance
- Sample: `datasets/swissprot/samples/sample_sequences.fasta`

Download:

```bash
curl -L -o datasets/swissprot/uniprot_sprot.fasta.gz \
  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

Load:

```python
import gzip
from Bio import SeqIO

with gzip.open("datasets/swissprot/uniprot_sprot.fasta.gz", "rt") as handle:
    records = list(SeqIO.parse(handle, "fasta"))
```

## PROSITE

- Source: https://ftp.expasy.org/databases/prosite/
- Location: `datasets/prosite/`
- Files:
  - `prosite.dat`: 2,731 entries
  - `prorule.dat`: 1,442 entries
  - `prosite_readme.txt`
- Task: known motifs, domains, profiles, and functional sites for feature-alignment baselines
- Sample: `datasets/prosite/samples/prosite_first_entries.txt`

Download:

```bash
curl -L -o datasets/prosite/prosite.dat https://ftp.expasy.org/databases/prosite/prosite.dat
curl -L -o datasets/prosite/prorule.dat https://ftp.expasy.org/databases/prosite/prorule.dat
```

## ELM Linear Motifs

- Source: http://elm.eu.org/
- Location: `datasets/elm/`
- Files:
  - `elms_index.tsv`: 354 motif class rows excluding comments
  - `instances.tsv`: 101 instance rows excluding comments in current export
- Task: eukaryotic short linear motif comparison and CAV/motif-localization baselines
- Samples: `datasets/elm/samples/`

Download:

```bash
curl -L -o datasets/elm/elms_index.tsv http://elm.eu.org/elms/elms_index.tsv
curl -L -o datasets/elm/instances.tsv http://elm.eu.org/instances.tsv
```

## Validation

Validation summary is saved in `datasets/dataset_summary.json`. The data were loaded or counted after download, and small samples were written for each dataset family.
