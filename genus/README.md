# Genus computations

This folder holds the genus computations cited in Appendix A of the paper. They concern two lists of graphs.

Table 1 is the 74 graphs T1, ..., T74 of the paper. The file `table1_graphs.txt` gives each one as its graph6 string, copied from item O1 of the supplemental data. Each of them has genus two. Item O2 gives an embedding of genus two for each, and the computations here show that none has genus one or zero.

Table 2 is every nontoroidal graph that the enumerations of the supplemental data use without another certificate of nontoroidality. These are the one-edge minors given by a verdict "proper nontoroidal minor" whose witness names only the deleted or contracted edge (including three in E3 whose recorded reason is this computation, and those of the whole graphs recorded in E8), the recorded nontoroidal graphs of the coverings E5 and E13, and the combinations of E20, and its auxiliary combinations, whose nontoroidal status excludes larger combinations. There are 815,417 such graphs counted with repetition and 312,410 up to isomorphism. `table2_graphs.txt` lists the graphs up to isomorphism, by index, with the canonical graph6 string of each (computed with nauty's `labelg`), its numbers of vertices and edges, and its number of uses. `table2_uses.json.gz` gives every use: the enumeration, the candidate and how the graph is obtained from it. `inputs_summary.json` gives the counts by enumeration. `build_inputs.py` makes these files from the supplemental data folder.

## The programs

PAGE (Metzger and Ulrigg, Discrete Mathematics 349, 2026) and multi_genus (Brinkmann, Ars Mathematica Contemporanea 22, 2022) both compute the genus of a graph exactly, and each also returns an embedding of that genus. They are used here as bundled with the Python package graph-genus 0.1.3 (https://github.com/SanderGi/Genus), whose sources for them are `PAGE/page.c` and `MultiGenus/multi_genus.c` of that repository.

`run_genus.py` runs the two programs at the same time. multi_genus is single-threaded and runs on one core: it computes the genus of T1, ..., T74 one at a time, and then the genus of every graph of Table 2, in one run of multi_genus per enumeration. PAGE runs on the other cores with its own threads (the environment variable THREADS, set to the number of cores minus one; PAGE uses at most eight): it computes the genus of T1, ..., T74 one at a time. PAGE is not run on the graphs of Table 2.

Each genus a program reports is checked against the embedding it returns: the faces of that rotation system are traced and Euler's formula gives its genus. T25 is two disjoint copies of K3,3; both programs need connected graphs, so it is computed component by component, using that the genus of a graph is the sum of the genera of its components (Battle, Harary, Kodama and Youngs, 1962).

## Running it

With Python 3.8 or later and the package graph-genus 0.1.3 installed (`pip install graph-genus==0.1.3`), run, in this folder,

```
python run_genus.py
```

It takes a few minutes. At the end it checks that both programs give genus two for every T_i and that multi_genus gives genus at least two for every graph of Table 2, and it writes the folder `results`: `TABLES.md` (both tables), `table1.csv`, `table2.csv`, `table2_genus.txt.gz` (the genus of each graph of Table 2, by index) and `run_info.json` (the machine and the versions used).
