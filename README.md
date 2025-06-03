# Protein-Reciprocal-Hits

_A toolkit for computing pairwise protein alignment scores using BLOSUM62 and identifying best reciprocal hits (BRHs) between human and chicken gene sets._

---

## Overview

Protein sequences from human and chicken are compared via a memoized global‐alignment scoring function (using the BLOSUM62 substitution matrix and affine gap penalties). All‐versus‐all pairwise alignment scores are computed, and best reciprocal hits (BRHs) are reported to suggest putative orthologous relationships between sex-determining genes of the two species.

---

## Repository Structure
`Protein-Reciprocal-Hits/`

├── `blosum62.py`

├── `humanChickenProteins.py.`

├── `sex_determination.py`

└── `README.md`

- `blosum62.py`
  - Defines the standard BLOSUM62 substitution matrix as a Python dictionary for scoring amino acid substitutions.

- `humanChickenProteins.py`
  - Supplies two gene lists (`humanGeneList`, `chickenGeneList`) and a dictionary (`geneD`) containing chromosome coordinates and protein sequences for each gene. Also includes smaller “sample” lists (`sampleHumanGeneList`, `sampleChickenGeneList`) for quick testing.
 
- `sex_determination.py`
  - Implements:
    1. A memoized global‐alignment scorer (`memoAlignScore`) that uses BLOSUM62 and a fixed gap penalty.
    2. Functions to compute all pairwise alignment scores between two gene lists.
    3. Functions to find the “closest match” for any given gene and to print best reciprocal hits.
    4. Example code showing how to run BRH analysis on sample data and, optionally, on the full gene lists.
   
---
## Prerequisites
- `Python 3.6` or higher
- No additional packages are required; all functionality relies on the standard library.

---
## Usage
**1. Clone or download this repository**

    git clone https://github.com/<achandrasek6>/Protein-Reciprocal-Hits.git
    cd HumanChickenOrthology

**2. Run BRH analysis on sample genes**

Execute:

    python3 sex_determination.py

- By default, the script computes alignment scores between the sample human and chicken genes defined in `humanChickenProteins.py`.
- It prints best reciprocal hits in this form:

      human --- chicken
      <chromosome> <start> <geneID> --- <chromosome> <start> <geneID>
      …
      
      chicken --- human
      <chromosome> <start> <geneID> --- <chromosome> <start> <geneID>
      …

**3. (Optional) Full BRH analysis**

- To repeat the analysis on the complete gene sets (`humanGeneList` vs. `chickenGeneList`):
- `Open sex_determination.py`
- Uncomment the call to `runBRH()` near the bottom of the file (and comment out or remove the `runBRHSample()` call).
- Save and re‐run:
  
      python3 sex_determination.py

---

## File Descriptions
### `blosum62.py`
  - Contains a single dictionary named `blosum62` mapping each ordered pair of amino acids `('A','R')` → score, etc. This is used by the alignment routine to score matches and mismatches.
### `humanChickenProteins.py`
  - Defines: `geneD` - a mapping from gene identifiers (e.g. `"h4"`, `"c8"`, etc.) to a tuple:

        (
          chromosome_string,
          start_position_int,
          end_position_int,
          protein_sequence_string
        )

- Each `protein_sequence_string` is the amino acid sequence of that gene’s protein.
  - `sampleHumanGeneList` & `sampleChickenGeneList`
    - Short lists of gene IDs for quick testing (e.g. [`"h4"`,`"h6"`,`"h9"`,`"h17"`] and [`"c8"`,`"c17"`,`"c19"`,`"c22"`]).
  - `humanGeneList` & `chickenGeneList`
    - Full sets of gene IDs to be used in the comprehensive BRH analysis.
### `sex_determination.py`
  - Implements the core logic:
    1. `memoAlignScore(S1, S2, gap, substitutionMatrix, memo)`
        - Computes a global alignment score between protein sequences `S1` and `S2` using recursion with memoization.
        - `gap` is the (negative) penalty for inserting a gap; here it is set to –9.
        - `substitutionMatrix` is passed in as the BLOSUM62 dictionary.
        - `memo` is a Python dict used to cache intermediate results and avoid redundant recursion.
    2. `allScores(geneList1, geneList2)`
        - Iterates over every pair `(gene1, gene2)` where gene1 comes from `geneList1` and `gene2` from `geneList2`.
        - Calls `memoAlignScore` on their protein sequences (`geneD[gene][3]`) with the fixed gap penalty of –9 and `blosum62`.
        - Returns a dict mapping each `(gene1, gene2)` tuple → alignment score (integer).
    3. closestMatch(geneName, allScoresD)
        - Looks through all keys in `allScoresD` to find the gene in the _other_ species with the highest score against `geneName`.
        - Returns that partner gene’s ID.
    4. `printBRH(geneName, allScoresD)`
        - Calls `closestMatch(geneName, allScoresD)` to find the top hit in the other species.
        - Then calls `closestMatch` again to see if that top hit’s best partner is `geneName` itself. If so, it prints both genes and their chromosomal coordinates in the format:
          
              <chrX> <startX> <geneName> --- <chrY> <startY> <bestPartner>
        - If no reciprocal hit is found, it prints nothing for that gene.
    5. `runBRHSample()`
        - Computes `allScoresD = allScores(sampleHumanGeneList, sampleChickenGeneList)`.
        - Prints:

              human --- chicken
              … (BRH lines for each human sample gene) …

              chicken --- human
              … (BRH lines for each chicken sample gene) …
    6. `runBRH()` _(commented out by default)_
       - Identical to `runBRHSample()`, but uses the full lists `humanGeneList` & `chickenGeneList`.
         
**At the bottom, the script calls `runBRHSample()`. To run the full dataset, swap that out for `runBRH()`.**
- _Note that this may take a few minutes to finish running._

--- 

## Biological Relevance
- Identifying best reciprocal hits between human and chicken proteins provides:
    - **Orthology Inference**
      - Detecting one‐to‐one orthologs helps transfer functional annotations (e.g., gene essentiality, pathway membership) from one species to the other.
    - **Comparative Genomics**
      - Chicken serves as a classic vertebrate model. Mapping orthologs highlights conserved developmental and physiological pathways.
    - **Sex Chromosome Evolution** _(Demonstrated in `sex_determination.py` script)_
      - By observing whether human X‐chromosome genes have reciprocal matches on chicken Z (and vice versa), insights into sex‐linked gene conservation or divergence can be gained.
 
---
## Troubleshooting

- **Recursion Limit Errors**
  - Very long protein sequences (over ~1,000 amino acids) may exceed Python’s default recursion limit. In `sex_determination.py`, near the top:
  
        import sys
        sys.setrecursionlimit(100000)
    ensures that deeply nested calls to `memoAlignScore` do not crash with a `RecursionError`.
  
- **Slow Runtime on Full Lists**
  - Computing all‐pairs alignment for large gene sets is inherently O(N₁ × N₂ × len(seq)²) in the worst case. For the full human vs. chicken list (~30–40 genes each), the script may take several minutes. For faster performance:
    - Test first on the short `sample…` lists.
    - Consider replacing memoized recursion with an iterative dynamic‐programming implementation in C or Cython for very large numbers of comparisons.

- **Missing or Malformed Data in `humanChickenProteins.py`**
  - Ensure that:
    - Every gene ID in `sampleHumanGeneList`, `humanGeneList`, `sampleChickenGeneList`, and `chickenGeneList` appears as a key in `geneD`.
    - Each entry in `geneD[gene]` is a 4‐tuple: `(chromosome, start, end, proteinSeq)`, where `proteinSeq` is a string of uppercase letters in {A, C, D, …, Y}.

- **Incorrect BLOSUM62 Keys**
  - The `blosum62` dictionary must contain entries for every ordered pair of amino acids, e.g. `('A','R')`, `('R','A')`, etc. If a substitution key is missing, scoring any alignment will raise a `KeyError`.

