# In-Silico Target Identification in Bovine Mastitis-Causing Pathogens

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1d0-paRZanuewG6yL09DlaJTx0vYxR0H3)

Computational subtractive genomics, reverse vaccinology, and druggability pipeline for identifying therapeutic targets against four major mastitis-causing bacterial pathogens:
- *Staphylococcus aureus* (Gram-positive)
- *Streptococcus uberis* (Gram-positive)
- *Escherichia coli* (Gram-negative)
- *Klebsiella pneumoniae* (Gram-negative)

---

## Repository Structure

```
FIAB2026_DT/
├── Proteome/                                       # Bacterial proteome FASTA files (.faa)
│   ├── Escherichia_coli_K-12.faa
│   ├── Klebsiella_pneumoniae.faa
│   ├── Staphylococcus_aureus.faa
│   └── Streptococcus_uberis.faa
│   └── OrthoFinder/Results_Oct27/                  # Precomputed OrthoFinder orthogroups
├── 02_Annotation/                                  # Core orthogroup extraction & CD-HIT NR1 proteome
├── 03_TargetHunt_Core/                             # Bovine & Human host non-homology BLASTP results
├── 04_Druggability/                                # DEG, TTD, ProtParam, & Gram-specific target sets
├── In_Silico_Target_Identification_Pipeline.ipynb  # Interactive Google Colab / Jupyter notebook
├── In_Silico_Target_Identification_Workflow.md    # Complete protocol and command-line reference
├── scripts/
│   └── generate_demo_data.py                       # Self-contained demo dataset generator
├── gram_positive_psortb.txt                        # Precomputed Gram-positive PSORTb predictions
├── gram_negative_psortb.txt                        # Precomputed Gram-negative PSORTb predictions
└── README.md
```

---

## Quick Start on Google Colab

Click the badge above or use the link below to open the pipeline directly in Google Colab:

👉 **[Open Workshop Pipeline in Google Colab](https://colab.research.google.com/drive/1d0-paRZanuewG6yL09DlaJTx0vYxR0H3)**

The notebook is configured to automatically pull this repository (`https://github.com/bioinfokushwaha/FIAB2026_DT.git`) and run every step with instant visualization and analysis in **Fast Workshop Demo Mode**.

---

## Pipeline Overview

1. **Pan-Pathogen Core Proteome:** OrthoFinder identifies proteins conserved across all 4 bacterial pathogens to discover broad-spectrum targets.
2. **Redundancy Removal:** CD-HIT clusters sequences at 100% identity to eliminate duplicate entries.
3. **Host Non-Homology:** Eliminating Bovine (*Bos taurus*) and Human (*Homo sapiens*) homologs avoids on-target toxicity and autoimmune cross-reactivity.
4. **Essentiality (DEG):** Ensures disrupting the target halts bacterial growth or kills the cell.
5. **Druggability (TTD / BindingDB):** Confirms target tractability and binding pocket feasibility.
6. **In-Vitro Stability (ProtParam):** Instability index $< 40$ prioritizes stable proteins suitable for recombinant expression and structural validation.
7. **Cellular Localization (PSORTb):** Identifies surface/membrane proteins (vaccines / antibodies) and cytoplasmic enzymes (small-molecule antibiotics).
