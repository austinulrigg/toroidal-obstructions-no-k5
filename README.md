# Toroidal Obstructions with No K5 Minor

This repository accompanies the paper "Toroidal Obstructions with No K5 Minor" by Austin Ulrigg. It holds the supplemental data of the paper and the genus computations of its Appendix A.

The folder `supplemental_data` holds every item the paper cites: O1–O4, W1–W21, R1–R11, E1–E34, TO1 and TO2, each under its own name in the subfolder given by its letters. Its `README.md` describes the files and the kinds of witness they contain, and `SHA256SUMS.txt` gives the SHA-256 checksum of every file.

The folder `genus` holds the genus computations. Both PAGE and multi_genus compute the genus of the 74 obstructions T1, ..., T74 (Table 1). multi_genus computes the genus of the 312,410 nontoroidal graphs, up to isomorphism, that the enumerations use without another certificate (Table 2). The results, with the running times, are in `genus/results/TABLES.md`, and `genus/README.md` says how to run the computation again.
