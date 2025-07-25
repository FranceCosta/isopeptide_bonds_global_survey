# Analyze the chemical environment surrounding isopeptide bonds

## Workflow for PDB structures
- In this case it is not possible to compare between positive and negative hits (the dataset only contains positive hits). Therefore, no need of alignment. Perform analysis on one go with: `python bin/proximity_analysis_pdb.py`

## Workflow for AFDB structures
1. Prepare structures for FoldMason alignment with `python bin/prepare_structures.py` considering +/- 30 aa from domain bounderies'
2. Align structures with `python bin/align_foldmason.py`
3. Find Isopeptide bond-like positions in negative entries with `python bin/map_to_msa.py` (consider +/- 30 aa from pfam domain start and end);
4. re-assign the domain boundaries using the rule: +/- 20 aa from first/last positions and cluster sequences at 30% (cov 80%). Consider only sequences with all 3 isopeptide bonds relative positions `python bin/cluster.py`
5. Perform structural analysis of surrounding amino acids `python bin/proximity_analysis_af_scan.py`;