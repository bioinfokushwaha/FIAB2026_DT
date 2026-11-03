# In-Silico Target Identification in Bovine Mastitis-Causing Pathogens

## Overview & Workflow

This tutorial describes an in-silico subtractive genomics and reverse vaccinology / drug-target identification pipeline for four major bovine mastitis-causing pathogens:
1. *Staphylococcus aureus* (Gram-positive)
2. *Streptococcus uberis* (Gram-positive)
3. *Escherichia coli* (Gram-negative)
4. *Klebsiella pneumoniae* (Gram-negative)

### Pipeline Flowchart

```
NCBI Proteomes (.faa)
       │
       ▼
  OrthoFinder ────► Core Orthogroups (present in all 4 species)
       │
       ▼
    CD-HIT ───────► Non-Redundant Core Proteome (100% identity cutoff)
       │
       ▼
 Bovine Host BLASTP ──► Remove Bovine Homologs (>40% identity or >40% coverage)
       │
       ▼
 Human Host BLASTP ───► Remove Human Homologs (>40% identity or >40% coverage)
       │
       ▼
   DEG BLASTP ────► Retain Essential Genes (≥40% identity and ≥40% coverage)
       │
       ▼
   TTD BLASTP ────► Retain Druggable Targets (≥40% identity and ≥40% coverage)
       │
       ▼
 ProtParam (BioPython) ─► Retain Stable Proteins (Instability Index < 40)
       │
       ▼
 Annotation Recovery ──► Map Protein IDs to Original FASTA Descriptions
       │
       ▼
 Gram Classification ──► Separate Gram-Positive vs Gram-Negative FASTA
       │
       ▼
     PSORTb ──────► Subcellular Localization (Cytoplasmic, Membrane, etc.)
```

---

## Prerequisites & Required Files

To execute this entire workflow from start to finish, the following input files, reference databases, and software tools are required:

### 1. Bacterial Proteome FASTA Files (Required in `Proteome/`)
You need the translated protein FASTA files (`.faa`) for the target bacterial strains:
* *Staphylococcus aureus* proteome (`.faa`)
* *Escherichia coli* proteome (`.faa`)
* *Streptococcus uberis* proteome (`.faa`)
* *Klebsiella pneumoniae* proteome (`.faa`)

> **Download source:** NCBI GenBank / RefSeq genomes or NCBI Datasets:
> `https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/`

### 2. Reference BLAST Databases (Required in `Database/`)
You need formatted NCBI BLAST protein databases (`makeblastdb -in <file.fa> -dbtype prot -out Database/<name>`):
* `Database/BT` — *Bos taurus* (bovine) host proteome database
* `Database/HS` — *Homo sapiens* (human) host proteome database
* `Database/DEGaa` — Database of Essential Genes (DEG) protein database
* `Database/TTD` — Therapeutic Target Database (TTD) protein database

### 3. Software Dependencies
* **OrthoFinder** (`orthofinder` or `/path/to/orthofinder.py`)
* **CD-HIT** (`cd-hit`)
* **NCBI BLAST+** (`blastp`, `makeblastdb`)
* **Python 3** with packages:
  * `pandas`
  * `matplotlib`
  * `biopython`

---

## Step 1: Directory Setup & Input Verification

### Step 1.1: Create Project Directory Structure

```bash
# Create organized workshop directory tree
mkdir -p 02_Annotation 03_TargetHunt_Core 04_Druggability Database Proteome results

# Verify directory structure
ls -d 02_Annotation 03_TargetHunt_Core 04_Druggability Database Proteome results
```

### Step 1.2: Place Bacterial Proteome Files
Upload your `.faa` files into the `Proteome/` directory.

### Step 1.3: Verify Uploaded FASTA Files
```bash
# Count total .faa files in Proteome/ (should be 4)
find Proteome/ -maxdepth 1 -type f -name "*.faa" | wc -l
```

### Step 1.4: Run OrthoFinder
Run OrthoFinder from the parent directory containing the `Proteome/` folder:

```bash
# Run OrthoFinder with 15 threads (adjust threads according to your CPU cores)
orthofinder -f Proteome -t 15 -a 15
```

---

## Step 2: Verification of OrthoFinder Results

### Step 2.1: Locate the Results Directory
```bash
# List OrthoFinder result folders
ls -d Proteome/OrthoFinder/Results_*/
```

### Step 2.2: Verify Output Files
```bash
# List contents of the results folder
ls Proteome/OrthoFinder/Results_*/
```

### Step 2.3: Verify the Orthogroups Table
```bash
# Locate the Orthogroups.tsv file
find Proteome/OrthoFinder/Results_*/ -name Orthogroups.tsv

# Display the first 3 rows of the orthogroups table
head -n 3 Proteome/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv

# Count the number of columns (expected: 5 -> 1 Orthogroup ID column + 4 species columns)
awk -F '\t' 'NR==1 {print NF}' Proteome/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv
```

### Step 2.4: Count Total Orthogroups
```bash
# Count orthogroups excluding header
awk 'END {print NR-1}' Proteome/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv
```

### Step 2.5: Verify Orthogroup Sequence Files
```bash
# Check the generated individual FASTA sequence files
ls Proteome/OrthoFinder/Results_*/Orthogroup_Sequences/ | head
```

### Step 2.6: Confirm Input File for Pangenome Analysis
```bash
# Ensure Orthogroups.tsv exists and is not empty
test -s Proteome/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv && echo "Orthogroups.tsv OK"
```

---

## Step 3: Pangenome Analysis Using Python

### Step 3.1: Check Python Environment & Required Libraries
```bash
python3 -c "import pandas as pd; import matplotlib; import Bio; print('pandas:', pd.__version__); print('matplotlib:', matplotlib.__version__); print('biopython:', Bio.__version__)"
```

### Step 3.2: Classify Core, Accessory, and Unique Orthogroups
This script dynamically finds `Orthogroups.tsv`, converts gene presence/absence into a binary matrix, and classifies orthogroups:
* **Core:** Present in all 4 species (`presence_count == 4`)
* **Accessory:** Present in 2 or 3 species (`2 <= presence_count < 4`)
* **Unique:** Present in only 1 species (`presence_count == 1`)

Run this Python script:

```python
import glob
import os
import pandas as pd

# Automatically locate the Orthogroups.tsv file
tsv_files = glob.glob("Proteome/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv")
if not tsv_files:
    raise FileNotFoundError("Could not find Orthogroups.tsv in Proteome/OrthoFinder/Results_*/Orthogroups/")

orthogroups_file = sorted(tsv_files)[-1]
print(f"Reading: {orthogroups_file}")

# Load Orthogroups table
df = pd.read_csv(orthogroups_file, sep="\t")

print("\n===== INPUT FILE CHECK =====")
print("Number of orthogroups:", len(df))
print("Number of species:", len(df.columns) - 1)
print("\nSpecies names:")
for species in df.columns[1:]:
    print(" -", species)

# Species columns (excluding 'Orthogroup')
species_cols = df.columns[1:]

# Convert gene presence/absence into 1 and 0
binary_df = df.copy()
binary_df[species_cols] = df[species_cols].notna().astype(int)

# Save presence/absence matrix
binary_df.to_csv("Ortho_pangenome_matrix.csv", index=False)

# Count species containing each orthogroup
presence_count = binary_df[species_cols].sum(axis=1)
num_species = len(species_cols)

# Classify orthogroups
core = binary_df[presence_count == num_species]
accessory = binary_df[(presence_count >= 2) & (presence_count < num_species)]
unique = binary_df[presence_count == 1]

# Count unique orthogroups for each species
unique_counts = unique[species_cols].sum()

# Save categorized tables
core.to_csv("Ortho_core_genes.csv", index=False)
accessory.to_csv("Ortho_accessory_genes.csv", index=False)
unique.to_csv("Ortho_unique_genes.csv", index=False)

# Display summary
print("\n===== PANGENOME SUMMARY =====")
print("Core orthogroups:     ", len(core))
print("Accessory orthogroups:", len(accessory))
print("Unique orthogroups:   ", len(unique))

print("\n===== UNIQUE ORTHOGROUPS PER SPECIES =====")
print(unique_counts.to_string())

# Validation
print("\n===== VALIDATION =====")
print("All orthogroups classified:", len(core) + len(accessory) + len(unique) == len(df))
```

### Step 3.3: Visualize Pangenome Distribution
Generate publication-quality plots (300 DPI) for pangenome composition and species-specific unique orthogroups:

```python
import matplotlib.pyplot as plt

# 1. Pangenome composition bar chart
counts = {"Core": len(core), "Accessory": len(accessory), "Unique": len(unique)}

plt.figure(figsize=(6, 5))
plt.bar(counts.keys(), counts.values(), color=["green", "blue", "red"])
plt.ylabel("Number of Orthogroups")
plt.title("Pangenome Composition")
plt.tight_layout()
plt.savefig("pangenome_composition.png", dpi=300)
plt.close()

# 2. Unique orthogroups per species bar chart
plt.figure(figsize=(8, 5))
unique_counts.plot(kind="bar", color="purple")
plt.ylabel("Number of Unique Orthogroups")
plt.title("Unique Orthogroups per Species")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("unique_orthogroups_per_species.png", dpi=300)
plt.close()

print("Both pangenome plots saved successfully!")
```

---

## Step 4: Extract Core Sequences & Remove Redundancy

### Step 4.1: Extract and Merge Core Orthogroup Sequences
Save the following bash script as `01_Extract_Core_Sequences.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

# Locate Orthogroup_Sequences directory dynamically
ORTHO=$(ls -d ./Proteome/OrthoFinder/Results_*/Orthogroup_Sequences | head -n 1)
OUT="02_Annotation"

mkdir -p "$OUT/Core_OrthoGP_Seq"

# Extract core orthogroup IDs (skip header)
tail -n +2 Ortho_core_genes.csv | cut -d',' -f1 > "$OUT/Core_Ortho_Ids"

# Empty or initialize destination FASTA file
: > "$OUT/Core_OrthoGp_proteome.fa"

# Copy and merge sequences
while read -r OG; do
    if [ -f "$ORTHO/$OG.fa" ]; then
        cp "$ORTHO/$OG.fa" "$OUT/Core_OrthoGP_Seq/"
        cat "$ORTHO/$OG.fa" >> "$OUT/Core_OrthoGp_proteome.fa"
    fi
done < "$OUT/Core_Ortho_Ids"

# Display summary counts
echo "Core orthogroups: $(wc -l < "$OUT/Core_Ortho_Ids")"
echo "Total proteins:   $(grep -c '^>' "$OUT/Core_OrthoGp_proteome.fa")"
echo "Extraction completed!"
```

Run the script:
```bash
bash 01_Extract_Core_Sequences.sh
```

### Step 4.2: Remove Redundant Sequences with CD-HIT (100% Identity)
Remove 100% identical duplicates among core proteins:

```bash
#!/usr/bin/env bash
set -euo pipefail

INPUT="02_Annotation/Core_OrthoGp_proteome.fa"
OUTPUT="02_Annotation/Core_OrthoGp_proteome_NR1.fa"

# Cluster at 100% identity (-c 1.0 -n 5)
cd-hit -i "$INPUT" -o "$OUTPUT" -c 1.0 -n 5

# Compare counts
echo "Proteins before CD-HIT: $(grep -c '^>' "$INPUT")"
echo "Proteins after CD-HIT:  $(grep -c '^>' "$OUTPUT")"
echo "CD-HIT completed!"
```

---

## Step 5: Host Non-Homology Screening (Bovine & Human)

Filter out bacterial proteins homologous to host proteins to prevent cross-reactivity and adverse autoimmune/toxic side effects.

### Thresholds:
* **BLASTP E-value:** $\le 10^{-20}$
* **Sequence Identity:** $> 40\%$
* **Query Coverage:** $> 40\%$
* **Action:** Exclude any protein matching host with Identity $> 40\%$ OR Query Coverage $> 40\%$. Retain non-homologs.

---

### Step 5.1 & 5.2: Bovine Host Screening (`Database/BT`)

```bash
#!/usr/bin/env bash
set -euo pipefail

INPUT="02_Annotation/Core_OrthoGp_proteome_NR1.fa"
OUT="03_TargetHunt_Core"
mkdir -p "$OUT"

# 1. BLASTP against bovine proteome
blastp -query "$INPUT" \
       -db Database/BT \
       -out "$OUT/Core_BT.tsv" \
       -outfmt "6 qseqid pident qcovhsp evalue" \
       -evalue 1e-20 \
       -num_threads 25

echo "Bovine BLASTP completed!"

# 2. Identify proteins with identity > 40% OR coverage > 40%
awk '$2>40 || $3>40 {print $1}' "$OUT/Core_BT.tsv" | sort -u > "$OUT/Core_BT_Exclude.ids"

# 3. Extract retained (non-homologous) proteins
python3 - <<'PY'
from Bio import SeqIO

out = "03_TargetHunt_Core"
with open(f"{out}/Core_BT_Exclude.ids") as f:
    excluded = set(f.read().splitlines())

records = [
    r for r in SeqIO.parse("02_Annotation/Core_OrthoGp_proteome_NR1.fa", "fasta")
    if r.id not in excluded
]

n = SeqIO.write(records, f"{out}/Core_BT_NonHomo.fasta", "fasta")
print("Bovine non-homologous proteins retained:", n)
PY
```

---

### Step 5.3: Human Host Screening (`Database/HS`)

Screen the bovine non-homologous proteins against the human proteome:

```bash
#!/usr/bin/env bash
set -euo pipefail

OUT="03_TargetHunt_Core"

# 1. BLASTP against human proteome
blastp -query "$OUT/Core_BT_NonHomo.fasta" \
       -db Database/HS \
       -out "$OUT/Core_HS.tsv" \
       -outfmt "6 qseqid pident qcovhsp evalue" \
       -evalue 1e-20 \
       -num_threads 25

echo "Human BLASTP completed!"

# 2. Identify proteins with identity > 40% OR coverage > 40%
awk '$2>40 || $3>40 {print $1}' "$OUT/Core_HS.tsv" | sort -u > "$OUT/Core_HS_Exclude.ids"

# 3. Extract final host-filtered non-homologous proteins
python3 - <<'PY'
from Bio import SeqIO

out = "03_TargetHunt_Core"
with open(f"{out}/Core_HS_Exclude.ids") as f:
    excluded = set(f.read().splitlines())

records = [
    r for r in SeqIO.parse(f"{out}/Core_BT_NonHomo.fasta", "fasta")
    if r.id not in excluded
]

n = SeqIO.write(records, f"{out}/Core_HS_NonHomo.fasta", "fasta")
print("Final host-filtered proteins retained:", n)
PY

# Verify protein counts at each stage
echo "=== Summary Counts ==="
echo -n "Initial non-redundant proteins: "
grep -c '^>' 02_Annotation/Core_OrthoGp_proteome_NR1.fa
echo -n "After bovine filtering:         "
grep -c '^>' 03_TargetHunt_Core/Core_BT_NonHomo.fasta
echo -n "After human filtering:          "
grep -c '^>' 03_TargetHunt_Core/Core_HS_NonHomo.fasta
```

---

## Step 6: Screening Against Database of Essential Genes (DEG)

Identify essential bacterial genes (survival-critical proteins) by searching against DEG.

### Thresholds:
* **BLASTP E-value:** $\le 10^{-20}$
* **Sequence Identity:** $\ge 40\%$
* **Query Coverage:** $\ge 40\%$
* **Action:** Retain proteins matching essential genes.

```bash
#!/usr/bin/env bash
set -euo pipefail

OUT="04_Druggability"
mkdir -p "$OUT"

# 1. BLASTP against DEG database
blastp -query 03_TargetHunt_Core/Core_HS_NonHomo.fasta \
       -db Database/DEGaa \
       -out "$OUT/Core_HS_DEG.tsv" \
       -outfmt "6 qseqid pident qcovhsp evalue" \
       -evalue 1e-20 \
       -num_threads 25

echo "DEG BLASTP completed!"

# 2. Select proteins meeting identity >= 40% AND coverage >= 40%
awk '$2>=40 && $3>=40 {print $1}' "$OUT/Core_HS_DEG.tsv" | sort -u > "$OUT/Core_HS_DEG.ids"
echo "Candidate essential protein IDs: $(wc -l < "$OUT/Core_HS_DEG.ids")"

# 3. Extract essential protein sequences
python3 - <<'PY'
from Bio import SeqIO

out = "04_Druggability"
with open(f"{out}/Core_HS_DEG.ids") as f:
    ids = set(f.read().splitlines())

records = [
    r for r in SeqIO.parse("03_TargetHunt_Core/Core_HS_NonHomo.fasta", "fasta")
    if r.id in ids
]

n = SeqIO.write(records, f"{out}/Core_HS_DEG.fasta", "fasta")
print("Candidate essential proteins extracted:", n)
PY
```

---

## Step 7: Druggability Screening Against Therapeutic Target Database (TTD)

Identify druggable targets among essential proteins by searching against TTD.

### Thresholds:
* **BLASTP E-value:** $\le 10^{-20}$
* **Sequence Identity:** $\ge 40\%$
* **Query Coverage:** $\ge 40\%$
* **Action:** Retain matching therapeutic targets.

```bash
#!/usr/bin/env bash
set -euo pipefail

OUT="04_Druggability"

# 1. BLASTP against TTD database
blastp -query "$OUT/Core_HS_DEG.fasta" \
       -db Database/TTD \
       -out "$OUT/Core_HS_DEG_TTD.tsv" \
       -outfmt "6 qseqid pident qcovhsp evalue" \
       -evalue 1e-20 \
       -num_threads 25

echo "TTD BLASTP completed!"

# 2. Select matching targets
awk '$2>=40 && $3>=40 {print $1}' "$OUT/Core_HS_DEG_TTD.tsv" | sort -u > "$OUT/Core_HS_DEG_TTD.ids"
echo "Candidate drug target IDs: $(wc -l < "$OUT/Core_HS_DEG_TTD.ids")"

# 3. Extract drug-target protein sequences
python3 - <<'PY'
from Bio import SeqIO

out = "04_Druggability"
with open(f"{out}/Core_HS_DEG_TTD.ids") as f:
    ids = set(f.read().splitlines())

records = [
    r for r in SeqIO.parse(f"{out}/Core_HS_DEG.fasta", "fasta")
    if r.id in ids
]

n = SeqIO.write(records, f"{out}/Core_HS_DEG_TTD.fasta", "fasta")
print("Candidate drug-target proteins extracted:", n)
PY
```

---

## Step 8: Physicochemical Characterization & Stability Screening (ProtParam)

Calculate physicochemical properties using Biopython's ProtParam and filter for stable proteins:
* **Stable:** Instability Index $< 40$
* **Unstable:** Instability Index $\ge 40$

```bash
#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

out = "04_Druggability"
input_fasta = f"{out}/Core_HS_DEG_TTD.fasta"

rows = []
for r in SeqIO.parse(input_fasta, "fasta"):
    try:
        # Strip internal stops/unknown chars if present
        cleaned_seq = str(r.seq).replace("*", "")
        p = ProteinAnalysis(cleaned_seq)
        rows.append([
            r.id,
            len(r.seq),
            round(p.molecular_weight(), 2),
            round(p.isoelectric_point(), 2),
            round(p.instability_index(), 2),
            round(p.gravy(), 3)
        ])
    except ValueError as e:
        print(f"Skipped {r.id}: {e}")

df = pd.DataFrame(rows, columns=[
    "Protein_ID", "Length", "Molecular_Weight",
    "Isoelectric_Point", "Instability_Index", "GRAVY"
])

# Filter stable proteins (Instability Index < 40)
stable = df[df["Instability_Index"] < 40]
stable.to_csv(f"{out}/Core_HS_DEG_TTD_ProtParam.csv", index=False)

print("Total analyzed: ", len(df))
print("Stable proteins:", len(stable))

# Extract stable FASTA sequences
stable_ids = set(stable["Protein_ID"])
stable_seqs = [r for r in SeqIO.parse(input_fasta, "fasta") if r.id in stable_ids]
n = SeqIO.write(stable_seqs, f"{out}/Core_HS_DEG_TTD_Stable.fasta", "fasta")
print(f"Stable FASTA written: {n} sequences -> {out}/Core_HS_DEG_TTD_Stable.fasta")
PY
```

---

## Step 9: Functional Annotation Recovery

Retrieve original protein descriptions from the source proteome `.faa` files for all retained stable proteins:

```bash
#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
import glob
import pandas as pd
from Bio import SeqIO

out = "04_Druggability"

# 1. Build dictionary of all original protein annotations
ann = {}
for f in glob.glob("Proteome/*.faa"):
    for r in SeqIO.parse(f, "fasta"):
        ann[r.id] = r.description

# 2. Map stable protein IDs to descriptions
stable_fasta = f"{out}/Core_HS_DEG_TTD_Stable.fasta"
rows = [
    (r.id, ann.get(r.id, "NOT_FOUND"))
    for r in SeqIO.parse(stable_fasta, "fasta")
]

# 3. Save annotation table
tsv_out = f"{out}/Core_HS_DEG_TTD_Stable_Annotation.tsv"
df_ann = pd.DataFrame(rows, columns=["Protein_ID", "Annotation"])
df_ann.to_csv(tsv_out, sep="\t", index=False)

print("Total stable proteins processed:", len(rows))
print("Annotated:                      ", sum(a != "NOT_FOUND" for _, a in rows))
print("Not found:                      ", sum(a == "NOT_FOUND" for _, a in rows))
print(f"Saved: {tsv_out}")
PY
```

---

## Step 10: Separation into Gram-Positive and Gram-Negative Targets

Separate target candidates according to bacterial envelope architecture:
* **Gram-Positive:** *Staphylococcus*, *Streptococcus*
* **Gram-Negative:** *Escherichia*, *Klebsiella*

```bash
#!/usr/bin/env bash
set -euo pipefail

python3 - <<'PY'
import pandas as pd
from Bio import SeqIO

out = "04_Druggability"
df = pd.read_csv(f"{out}/Core_HS_DEG_TTD_Stable_Annotation.tsv", sep="\t")
seqs = list(SeqIO.parse(f"{out}/Core_HS_DEG_TTD_Stable.fasta", "fasta"))

categories = [
    ("GramPositive", r"Staphylococcus|Streptococcus"),
    ("GramNegative", r"Escherichia|Klebsiella")
]

for name, pattern in categories:
    ids = set(df.loc[df["Annotation"].str.contains(pattern, case=False, na=False), "Protein_ID"])
    matched_seqs = [r for r in seqs if r.id in ids]
    out_file = f"{out}/Core_HS_DEG_TTD_Stable_{name}.fasta"
    n = SeqIO.write(matched_seqs, out_file, "fasta")
    print(f"{name} sequences extracted: {n} -> {out_file}")
PY

# Verify generated FASTA files
ls -lh 04_Druggability/Core_HS_DEG_TTD_Stable_GramPositive.fasta \
       04_Druggability/Core_HS_DEG_TTD_Stable_GramNegative.fasta
```

---

## Step 11: Subcellular Localization (PSORTb Web Server)

Subcellular localization determines whether a target is surface-exposed/membrane-associated (ideal vaccine candidates or easily accessible drug targets) or cytoplasmic (requires membrane-permeable small-molecule inhibitors).

1. Open the **PSORTb v3.0** web server:
   [http://www.psort.org/psortb/](http://www.psort.org/psortb/)
2. For Gram-positive candidates:
   * Select **Bacteria** $\rightarrow$ **Gram-Positive**
   * Upload `04_Druggability/Core_HS_DEG_TTD_Stable_GramPositive.fasta`
   * Submit analysis.
3. For Gram-negative candidates:
   * Select **Bacteria** $\rightarrow$ **Gram-Negative**
   * Upload `04_Druggability/Core_HS_DEG_TTD_Stable_GramNegative.fasta`
   * Submit analysis.
4. Record predicted localizations:
   * Cytoplasmic
   * Cytoplasmic Membrane
   * Periplasmic (Gram-negative)
   * Outer Membrane (Gram-negative)
   * Extracellular
