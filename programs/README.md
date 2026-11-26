# Generating the supplemental data

The programs below generate the rotation tables R1–R12, the enumerations E1–E30, the split search TO1–TO2 and the catalogue N1. Run the commands from the root of the repository.

## Setup

Use Python 3.11 or later. The package versions in [requirements.txt](requirements.txt) are those used for the published data. The version of NetworkX can affect which Kuratowski subdivision is returned.

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r programs/requirements.txt
```

The package `graph-genus` includes PAGE and `multi_genus`. N1 also requires nauty's `labelg` (or `nauty-labelg`) on the search path, or its path in `LABELG`. E4 and E12 use `drat-trim` to check their DRAT proofs. The programs look for it in `DRAT_TRIM` or on the search path, then try to compile [drat-trim.c](../tools/drat-trim.c) with `cc`. If it is unavailable, they still write the proofs and report that they were not checked.

## Rotation tables and enumerations

Use an empty output folder for a new run. The following commands regenerate all tables and enumerations and compare them with the published files.

```sh
out="$HOME/regenerated"
python programs/rotation_tables/make_rotation_tables.py "$out"
python programs/compare_files.py supplemental_data "$out" rotations

python programs/enumerations/make_enumerations.py "$out" --E4-covering-of-data
python programs/compare_files.py supplemental_data "$out" enumerations
```

To run individual items, give their labels after the output folder, for example `R1 R5` or `E3 E11`. For each rotation system, the generator computes the genus and records the attachments on each facial walk in boundary order. It then checks the claims listed in the table. The paper passages used by the generators are in [paper_text.json](paper_text.json). [Supplemental references in the paper](../PAPER-REFERENCES.md) identifies the items used in earlier proofs.

Each enumeration has a `generate_E<k>.py` file describing its construction and tests. E1–E29 keep their working files in `$out/work/E<k>` and write their results to `$out/enumerations/E<k>`. Rerunning the command skips completed items. Most interrupted items resume from their saved records. To repeat a completed item, use a fresh output folder or remove its working folder. Parallel enumerations use one fewer worker than the number of cores by default. Set `ENUM_WORKERS` to choose the number. E30 runs PAGE on the 24 attachment-edge completions of T(z), S(x;y) and Q(x;y,z) and writes its results directly to `$out/enumerations/E30`.

E4 can find different coverings depending on which parallel search finishes first. The option `--E4-covering-of-data` uses the published covering and recomputes its DRAT proof, torus embeddings and graph classifications. Omit the option in a fresh run to find a covering independently. The [genus computations](../genus/README.md) refer to the published covering.

`compare_files.py` compares files byte for byte, then compares the contents of JSON files that differ, ignoring formatting and running times. PAGE or `multi_genus` can return a different embedding of the same graph. The program reports differences confined to rotation systems separately. To check regenerated witnesses, replace the corresponding files in a copy of `supplemental_data`, remove `SHA256SUMS.txt` from that copy and run the [checker](../checker/README.md) on it with `--items`.

## Topological obstructions

```sh
python programs/topological_search/topological_search.py supplemental_data "$out/no_K5_topological_obstructions"
python programs/topological_search/compare.py supplemental_data "$out/no_K5_topological_obstructions"
```

Starting from T1, …, T74 in O1, the search splits vertices as in Theorem 8.1. The comparison matches the 381 obstructions in TO1 and the 7,006 graphs in TO2 up to isomorphism and compares the verdict for each graph in TO2. Different runs can give different branch sets or embeddings for the same verdict.

On the Mac used for the reported runs, the rotation tables took about 20 minutes, the split search a few minutes and the enumerations three to four hours.

## The additional obstructions in N1

```sh
python programs/new_obstructions/make_N1.py "$out"
cmp supplemental_data/new_obstructions/N1.json "$out/new_obstructions/N1.json"
```

The program builds N1 from the published website records and retains the original certificates for N1.1–N1.13. It computes canonical forms and lists the entries of the Myrvold–Woodcock collection with matching graph invariants.

With `--recompute`, `multi_genus` computes the genus of all 120 graphs, every edge deletion, every edge contraction of the 22 minor obstructions, and one nontoroidal contraction for each of the other 98 graphs. With `--page`, it performs the same checks and also runs PAGE on the 120 graphs. Both options write the original witnesses.
