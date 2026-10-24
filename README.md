# Toroidal Obstructions with No K5 Minor

This repository accompanies the paper "Toroidal Obstructions with No K5 Minor" by Austin Ulrigg. It holds the supplemental data of the paper, the genus computations of its Appendix A, the programs that make the data, and a checker for the data.

The folder `supplemental_data` holds every item the paper cites: O1–O4, W1–W21, R1–R12, E1–E34, TO1 and TO2, each under its own name in the subfolder given by its letters. Its `README.md` describes the files and the kinds of witness they contain, and `SHA256SUMS.txt` gives the SHA-256 checksum of every file.

The folder `genus` holds the genus computations. Both PAGE and multi_genus compute the genus of the 74 obstructions T1, ..., T74 (Table 1). multi_genus computes the genus of the 312,410 nontoroidal graphs, up to isomorphism, that the enumerations use without another certificate (Table 2). The results, with the running times, are in `genus/results/TABLES.md`, and `genus/README.md` says how to run the computation again.

The folder `checker` holds a program that checks every witness of the supplemental data, using only the files of `supplemental_data`. Its `README.md` says what it checks and how to run it. The folder `tools` holds the source of drat-trim (https://github.com/marijnheule/drat-trim, MIT license), which the checker uses for the DRAT proofs of E4, E5 and E13.

The folder `programs` holds the programs that make the supplemental data: the rotation tables, the enumerations and the topological search. Its `README.md` says how to run them and how to compare what they write with `supplemental_data`.

The Python packages used by the checker and the genus computations are listed in `requirements.txt`; the programs need the exact versions listed in `programs/requirements.txt`.
