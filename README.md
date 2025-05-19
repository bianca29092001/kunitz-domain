# A Bioinformatics Workflow for Modeling and Detecting the Kunitz Domain Using Profile HMMs 

This repository implements a **computational pipeline** to build and evaluate a Profile Hidden Markov Model (HMM) for the **Kunitz-type protease inhibitor domain** (Pfam: PF00014).

The workflow integrates **sequence extraction**, **redundancy reduction**, **multiple sequence alignment**, **HMM construction**, and **model evaluation** using curated positive and negative datasets.

Developed as part of the *Laboratory of Bioinformatics 1* course (MSc in Bioinformatics, University of Bologna).

---

##  Table of Contents

- [Objectives](#objectives)
- [Requirements](#requirements)
- [Pipeline Overview](#pipeline-overview)
  - [1. Data Extraction](#1-data-extraction)
  - [2. Redundancy Filtering](#2-redundancy-filtering)
  - [3. Multiple Sequence Alignment](#3-multiple-sequence-alignment)
  - [4. HMM Construction](#4-hmm-construction)
  - [5. Dataset Preparation](#5-dataset-preparation)
  - [6. Model Evaluation](#6-model-evaluation)
- [Output](#output)
- [Author](#author)


---

## Objectives

The goal of this project is to build and evaluate a profile Hidden Markov Model (HMM) for the detection of Kunitz-type protease inhibitor domains (Pfam: PF00014). The main steps include:

- Extracting Kunitz-domain sequences from PDB and UniProt
- Removing redundancy to ensure representative and diverse data
- Building a structure-based multiple sequence alignment
- Creating a profile HMM from curated sequences
- Testing the model on filtered positive and negative datasets
- Assessing classification performance using hmmsearch outputs

---

##  Requirements

Create and activate a conda environment:

```bash
conda create -n kunitz_env python=3.10
conda activate kunitz_env
conda install -c bioconda cd-hit blast hmmer biopython
```

---
## Pipeline Overview
### 1. Data Extraction
Initial data was collected from the RCSB PDB and UniProt databases.

- **`rcsb_pdb_custom_report_*.csv`**:  
  A report downloaded from the RCSB PDB website using the advanced search interface.  
  It contains metadata such as resolution, chain ID, sequence length, and Pfam domain annotation (PF00014), used to filter for high-quality Kunitz-domain entries.

- **`pdb_kunitz_customreported.fasta`**:  
  The corresponding amino acid sequences of the PDB chains listed in the CSV file.  
  These are the raw input sequences for clustering and alignment.

Filtering criteria applied:
- Resolution ≤ 3.5 Å
- Sequence length between 45–80 amino acids
- Presence of the Kunitz-type domain (Pfam ID: PF00014)

### 2. Redundancy Filtering

To ensure the dataset used to train the HMM is not biased by overrepresented or highly similar sequences, we applied a redundancy reduction step using **CD-HIT**:

- **Input file**: `pdb_kunitz_customreported.fasta`
- **Tool used**: `cd-hit` with a 90% sequence identity threshold

```bash
cd-hit -i pdb_kunitz_customreported.fasta -o pdb_kunitz_customreported_filtered.fasta -c 0.9
```
CD-HIT clusters sequences based on sequence identity and outputs:

A filtered FASTA file `pdb_kunitz_customreported_filtered.fasta` with only representative sequences.

A cluster file `pdb_kunitz_customreported_filtered.clstr` listing all sequences and their cluster assignments.

We then extracted the representative sequences (the first entry in each cluster) using:

```bash
clstr2txt.pl pdb_kunitz_customreported_filtered.clstr > pdb_kunitz.clusters.txt
awk '$5 == 1 {print $1}' pdb_kunitz.clusters.txt > pdb_kunitz_rp.ids
grep -A 1 -Ff pdb_kunitz_rp.ids pdb_kunitz_customreported.fasta > pdb_kunitz_rp.fasta
```
`pdb_kunitz_rp.ids`: IDs of non-redundant representative sequences

`pdb_kunitz_rp.fasta`: FASTA file containing only the representative sequences to be used for alignment and model building


### 3. Multiple Sequence Alignment
The representative sequences were aligned using a structure-aware tool (e.g., PDBeFold) to ensure functional and structural consistency.

Output alignment file: `pdb_kunitz_rp.ali`

Converted to a proper FASTA-like alignment format:
```bash
awk '{if (substr($1,1,1)==">") {print "\n" toupper($1)} else {printf "%s", toupper($1)}}' pdb_kunitz_rp.ali > pdb_kunitz_rp_formatted.ali
```
This formatted file is compatible with HMMER for model building.

### 4. HMM Construction
The profile Hidden Markov Model was built using HMMER:
```bash
hmmbuild structural_model.hmm pdb_kunitz_rp_formatted.ali
```
- `structural_model.hmm` contains the trained model.

This model captures conserved sequence features of the Kunitz domain and can be used to classify new proteins.


### 5. Dataset Preparation

To evaluate the performance of the HMM, we prepared **four test sets** from UniProt sequences:

- **Two positive sets** containing proteins annotated with the Kunitz domain (Pfam: PF00014)
- **Two negative sets** containing proteins with no Kunitz annotation

#### Positive Set Preparation

5.1. **Remove redundancy with training set using BLAST**  
   The `pdb_kunitz_rp.fasta` file (used to train the model) was searched against the full UniProt test set:

   ```bash
   makeblastdb -in all_kunitz.fasta -dbtype prot
   blastp -query pdb_kunitz_rp.fasta -db all_kunitz.fasta -out pdb_kunitz_nr_23.blast -outfmt 7
   ```
5.2. **Extract IDs to remove**
   Based on identity ≥95% and alignment length ≥50%:
   ```bash
   grep -v "^#" pdb_kunitz_nr_23.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > to_remove.ids
   ```
5.3. **Select remaining (non-redundant) IDs:**
  ```bash
    grep ">" all_kunitz.fasta | cut -d "|" -f 2 > all_kunitz.id
    comm -23 <(sort all_kunitz.id) <(sort to_remove.ids) > to_keep.ids
  ```
5.4. **Extract sequences:**
   ```bash
   python scripts/get_seq.py to_keep.ids data/raw/uniprot_sprot.fasta > data/processed/ok_kunitz.fasta
   ```
5.5. **Randomize and split into 2 positive sets:**
  ```bash
  sort -R to_keep.ids > random_ok_kunitz.ids
  head -n 183 random_ok_kunitz.ids > pos_1.ids
  tail -n 183 random_ok_kunitz.ids > pos_2.ids

  python scripts/get_seq.py pos_1.ids data/raw/uniprot_sprot.fasta > pos_1.fasta
  python scripts/get_seq.py pos_2.ids data/raw/uniprot_sprot.fasta > pos_2.fasta
  ```
#### Negative Set Preparation
5.6. **Extract all UniProt IDs:**
   ```bash
   grep ">" data/raw/uniprot_sprot.fasta | cut -d "|" -f 2 > sp.id
   ```
5.7. **Remove all Kunitz positives from the full set:**
   ```bash
   comm -23 <(sort sp.id) <(sort all_kunitz.id) > sp_negs.ids
   ```   
5.8. **Randomize and split into 2 negative sets:**
  ```bash
  sort -R sp_negs.ids > random_sp_negs.ids
  head -n 286288 random_sp_negs.ids > neg_1.ids
  tail -n 286288 random_sp_negs.ids > neg_2.ids

  python scripts/get_seq.py neg_1.ids data/raw/uniprot_sprot.fasta > neg_1.fasta
  python scripts/get_seq.py neg_2.ids data/raw/uniprot_sprot.fasta > neg_2.fasta

  ``` 



### 6. Model Evaluation
The HMM was tested using hmmsearch on each set:
```bash
hmmsearch -Z 1000 --max --tblout pos_1.out structural_model.hmm pos_1.fasta
hmmsearch -Z 1000 --max --tblout pos_2.out structural_model.hmm pos_2.fasta
hmmsearch -Z 1000 --max --tblout neg_1.out structural_model.hmm neg_1.fasta
hmmsearch -Z 1000 --max --tblout neg_2.out structural_model.hmm neg_2.fasta
```
Each output file `.out` contains per-sequence scores and e-values. These were parsed into a tabular format for downstream classification:

```bash
grep -v "^#" pos_1.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_1.class
grep -v "^#" pos_2.out | awk '{split($1,a,"|"); print a[2]"\t1\t"$5"\t"$8}' > pos_2.class
grep -v "^#" neg_1.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_1.class
grep -v "^#" neg_2.out | awk '{split($1,a,"|"); print a[2]"\t0\t"$5"\t"$8}' > neg_2.class
```
To evaluate the classification performance of the HMM, the `.class` files were merged into two global test sets:

```bash
cat pos_1.class neg_1.class > set_1.class
cat pos_2.class neg_2.class > set_2.class
```
A custom Python script (performance.py) was used to calculate:

- Accuracy

- Precision

- Recall (Sensitivity)

- Specificity

- F1 score

- ROC AUC

Optimal threshold based on bit score
```bash
python scripts/performance.mcc.py set_1.class set_2.class
```
The script reads the two sets, plots score distributions, and computes metrics at multiple thresholds, helping visualize the model’s ability to discriminate between positives and negatives.
To further evaluate the HMM and determine the optimal classification threshold, several Python scripts were used.

#### Fixed-threshold evaluation

The performance was tested at a specific e-value threshold (e.g., 1e-5):

```bash
python3 performance.mcc.py set_1.class 1e-5
python3 performance.mcc.py set_2.class 1e-5
```

Confusion matrix and ROC

Additional scripts were used to visualize classification performance:
```bash
python3 confusion_matrix.py
python3 roc.py
```
Threshold sweep (1e-1 to 1e-10)
To evaluate the effect of varying the e-value threshold on performance, a loop was used:
```bash
for i in $(seq 1 10); do
    python3 performance.mcc.py set_1.class 1e-$i
done | sort -nrk 6 > performance_set1_thresholds.txt

for i in $(seq 1 10); do
    python3 performance.mcc.py set_2.class 1e-$i
done | sort -nrk 6 > performance_set2_thresholds.txt
```
These loops:

Test thresholds from 1e-1 to 1e-10

Sort results by the Matthews Correlation Coefficient (MCC)

Output the best-performing thresholds to:

- `performance_set1_thresholds.txt`

- `performance_set2_thresholds.txt`

This allows identification of the optimal cutoff for domain detection, balancing sensitivity and specificity.



## Output

The following files are generated throughout the pipeline:

| File | Description |
|------|-------------|
| `pdb_kunitz_rp.fasta` | Non-redundant PDB sequences used for MSA |
| `pdb_kunitz_rp_formatted.ali` | Formatted alignment for HMM building |
| `structural_model.hmm` | Final HMM profile for the Kunitz domain |
| `ok_kunitz.fasta` | Filtered UniProt sequences (non-redundant positives) |
| `pos_1.fasta`, `pos_2.fasta` | Positive datasets for evaluation |
| `neg_1.fasta`, `neg_2.fasta` | Negative datasets for evaluation |
| `*.out` | Raw HMMER search results |
| `*.class` | Parsed results: sequence ID, label, score, e-value |
| `set_1.class`, `set_2.class` | Merged class files for evaluation |
| `performance_set1_thresholds.txt` | MCC results across thresholds for set 1 |
| `performance_set2_thresholds.txt` | MCC results across thresholds for set 2 |
| `*.png`| Visual outputs (ROC curves, confusion matrices, MCC) |

## Author

**Bianca Mastroddi**  
MSc Student in Bioinformatics  
University of Bologna  
Laboratory of Bioinformatics 1 — A.Y. 2024/2025






