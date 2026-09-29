# Manual database/web steps

Several analyses in the manuscript use web databases or GUI applications and are
not fully represented by a single R script.

## STRING
The SASP-related candidate-gene list was submitted to STRING to obtain
protein-protein association information using the confidence threshold reported
in the manuscript.

## miRNet and Cytoscape
Candidate genes were queried in miRNet for predicted miRNA-gene and TF-gene
relationships. Exported interaction tables were visualized in Cytoscape.

## GeneMANIA
Candidate genes were submitted to GeneMANIA to retrieve predicted functional
associations and related genes.

## Compound enrichment
The original workflow used the R package enrichR and the DSigDB library.

## Molecular docking
Protein structures and compound structures were prepared as described in the
manuscript. Docking was performed with CB-Dock2. These web-server steps are
not automated in this repository.

## Provenance
The authors' original code archive contained machine-specific absolute paths and
an earlier GSE134697 validation workflow. The scripts in ../scripts document the
revised GSE190451 validation and clinical-association analyses in portable form.
