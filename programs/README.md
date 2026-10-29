# Programs

These are the programs that make the supplemental data of the paper. `rotation_tables` makes the rotation tables R1–R11, `enumerations` makes the enumerations E1–E34, `new_obstructions` makes N1, and `topological_search` makes TO1 and TO2, the search for the topological obstructions in the proof of Theorem 8.1. Each program writes its files in the same layout as the folder `supplemental_data`, so that what it writes can be compared with that folder file by file.

## Setting up

The programs need Python 3.11 or later and three packages: graph-genus 0.1.3 (which brings the genus programs PAGE and multi_genus), networkx 3.6.1 and python-sat 1.9.dev17. Some outputs depend on the exact version of a package (for example, which Kuratowski subgraph networkx returns), so `requirements.txt` in this folder fixes all three. In a terminal, in the top folder of the repository, run

```
python3 -m venv .venv
source .venv/bin/activate
pip install -r programs/requirements.txt
```

The first line makes a virtual environment, a private folder `.venv` for the packages; it is needed only once. The second line switches the terminal to it, and is needed in every new terminal window.

E5 and E13 also check their DRAT proofs with drat-trim, found in the same way as by the checker: through the environment variable `DRAT_TRIM`, or on the PATH, or by compiling `tools/drat-trim.c` the first time it is needed (this needs a C compiler, which on a Mac comes with the Xcode command line tools; without drat-trim the proofs are written but not checked).

## The folders

`rotation_tables/make_rotation_tables.py` runs through every rotation system of the graphs of each table and writes R1–R11 to `OUTPUT_FOLDER/rotations`. A statement of the paper about a table that fails stops the program.

`topological_search/topological_search.py` starts from T1, ..., T74, splits vertices in every way and writes TO1.json and TO2.json. `topological_search/compare.py` compares them with the data folder up to isomorphism: the same 381 topological obstructions and the same 7,006 graphs with the same verdicts. The witnesses (branch sets and embeddings) may differ from those of the data folder, since the search may find other valid ones.

`enumerations/make_enumerations.py` runs the program `generate_E<k>.py` of each enumeration, which writes its own records to `OUTPUT_FOLDER/work/E<k>`, and then writes the files of the item to `OUTPUT_FOLDER/enumerations/E<k>` (with `package.py`). The docstring at the top of each `generate_E<k>.py` says what it enumerates, where the paper uses it and how it classifies each graph. The other modules of the folder are shared: `enumeration.py` builds the bridges of Construction 3.17 and classifies whole graphs, `joint_test.py` is the test of Lemma 4.3, `covering.py` is the covering used by E5 and E13, and so on, each with its own docstring.

`new_obstructions/make_N1.py` takes the thirteen graphs of N1 and writes `OUTPUT_FOLDER/new_obstructions/N1.json`: the genus of each graph by multi_genus (and also by PAGE with the option `--page`), a minimum genus embedding of the graph and of each of its edge deletions and contractions by multi_genus, K5 branch sets, the connectivity, the canonical form (it needs nauty's `labelg`) and the comparison with the collection in O4. Run it from the top folder of the repository with `python programs/new_obstructions/make_N1.py OUTPUT_FOLDER`.

`compare_files.py` compares the files of an output folder with the data folder: first byte for byte, then, for JSON files whose bytes differ, by content, so that a difference in formatting alone is reported as such. For a file whose content differs, it counts the places that differ and says whether they all lie in rotation systems (see the last section).

`paper_text.json` holds the passages of the paper that the files of the items quote.

## Running everything

From the top folder of the repository, with the virtual environment active, and with `~/regenerated` as the output folder (any empty folder will do):

```
cd programs/rotation_tables
python make_rotation_tables.py ~/regenerated
python ../compare_files.py ../../supplemental_data ~/regenerated rotations

cd ../topological_search
python topological_search.py ../../supplemental_data ~/regenerated/topological
python compare.py ../../supplemental_data ~/regenerated/topological

cd ../enumerations
python make_enumerations.py ~/regenerated --E5-covering-of-data
python ../compare_files.py ../../supplemental_data ~/regenerated enumerations
```

Some enumerations can be run alone by naming them, as in `python make_enumerations.py ~/regenerated E4 E12`. The enumerations use all cores but one; the environment variable `ENUM_WORKERS` sets another number of processes. If a run is interrupted, running the same command again continues it: most of the programs keep their progress in their folder `~/regenerated/work/E<k>` and continue where they stopped, and a few start their item again. A program that has finished is not run again; to run it again, delete its folder `~/regenerated/work/E<k>`.

On a recent Mac, the rotation tables take about twenty minutes (R8 about a quarter of an hour), the topological search a few minutes, and the enumerations three to four hours in all; E22, E30, E9, E5 and E23 take the longest.

## What the comparison shows

The programs fix every choice that could vary from one run to another (the seeds of the random searches, the order of sets, and which of several isomorphisms or minor models is written), so that a run writes the files of the data folder again, and `compare_files.py` names every file that differs. There are two exceptions.

E5 covers its family of graphs with several processes at once, and which graphs it records depends on the order in which the processes finish. A run without the option `--E5-covering-of-data` therefore records a different, equally valid, family of graphs, and the genus computations in the folder `genus` (Table 2 of the paper) do not cover the nontoroidal graphs of that family. With the option, `generate_E5.py` reads the graphs recorded in the data folder instead of covering anew, and does everything after the covering again: the DRAT proof, the torus embeddings and the classification of every recorded graph.

The torus embeddings (rotation systems of genus at most one) of the toroidal graphs can differ from those of the data folder. PAGE, which finds most of them, searches with several threads at once, and for a few graphs the embedding found first changes from run to run; multi_genus, which finds those of E5, can return different embeddings of the same graph on different computers. Every such embedding is valid, and a comparison that reports differences only in the rotation systems of toroidal records has found nothing else.

The checker in the folder `checker` checks the files of the data folder. It can also check a regenerated item: copy the item's folder into a copy of `supplemental_data`, delete that copy's `SHA256SUMS.txt` (the checker then skips the checksums) and run the checker on the copy with `--items`.
