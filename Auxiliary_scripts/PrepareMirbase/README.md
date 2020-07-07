# Prepare miRBase data for MirCure


This Shiny App allows to transform miRBase miRNA annotations to the correct format to be used as input for MirCure.

- Given the [MirBase gff](ftp://mirbase.org/pub/mirbase/CURRENT/genomes) file for any species, this App will generate the 3 necessary gff files for MirCure.
- The 3 gff files that MirCure requires correspond to; precursors, the 5' arms, and the 3' arms coordinates.
- The 3 gff files contain same miRNAs in the same order, with the same ID indicating, and a suffix indicating
 "_3p" or "5_" arm.
- When one arm annotation is missing from  miRBase, this App will predict it based on the length of the annotated arm and position witin precursor.
