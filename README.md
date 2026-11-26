# Toroidal Obstructions with No K5 Minor

This repository accompanies the paper "Toroidal Obstructions with No K5 Minor" by Austin Ulrigg. It holds the supplemental data of the paper, the genus computations of its Appendix A, the programs that make the data, and a checker for the data.

The [guide to the paper's supplemental references](PAPER-REFERENCES.md) gives the place in the paper where each item is used and links to its data. The folder `supplemental_data` contains O1–O4, W1–W16, R1–R12, E1–E30, TO1, TO2 and N1. N1 includes the 120 topological obstructions in the website release of November 11, 2026 that are absent from the Myrvold–Woodcock collection of 250,815 topological obstructions. Of these 120 graphs, 22 are minor obstructions. The [data guide](supplemental_data/README.md) describes the files and witnesses, and [SHA256SUMS.txt](supplemental_data/SHA256SUMS.txt) gives their checksums.

The folder `genus` holds the genus computations. Both PAGE and multi_genus compute the genus of the 74 obstructions T1, ..., T74. multi_genus computes the genus of the 311,739 nontoroidal graphs, up to isomorphism, that the enumerations use without another certificate. The results, with the running times, are in `genus/results/TABLES.md`, and `genus/README.md` says how to run the computation again.

The folder `checker` holds a program that checks every witness of the supplemental data, using only the files of `supplemental_data`. Its `README.md` says what it checks and how to run it. The folder `tools` holds the source of drat-trim (https://github.com/marijnheule/drat-trim, MIT license), which the checker uses for the DRAT proofs of E3, E4 and E12.

The folder `programs` holds the programs that make the supplemental data: the rotation tables, the enumerations and the topological search. Its `README.md` says how to run them and how to compare what they write with `supplemental_data`.

The Python packages used by the checker and the genus computations are listed in `requirements.txt`; the programs need the exact versions listed in `programs/requirements.txt`.
