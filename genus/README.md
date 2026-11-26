# Genus computations

PAGE and `multi_genus` both give genus two for each of the 74 obstructions in [O1](../supplemental_data/no_K5_obstructions/O1.json). The graphs are listed in [obstructions.txt](obstructions.txt), and [O2](../supplemental_data/no_K5_obstructions/O2.json) gives a genus-two embedding of each. The [tables](results/TABLES.md) give the results and running times.

## Graphs used in the enumerations

The enumerations use `multi_genus` to establish nontoroidality in three cases.

- An edge deletion or contraction gives a proper nontoroidal minor, and the witness gives only the operation or refers to the genus computation. This includes the three genus checks in E2 and the minors of the whole graphs in E7.
- E4 and E12 use recorded nontoroidal graphs to exclude combinations containing them.
- E19 uses nontoroidal combinations, including auxiliary combinations, to exclude larger ones.

These cases give 814,514 uses of 311,739 graphs up to isomorphism. The [graph list](enumeration_graphs.txt) gives each graph's index, canonical graph6 string, numbers of vertices and edges, and number of uses. The [list of uses](enumeration_graph_uses.json.gz) identifies the enumeration, candidate and graph operation where applicable. The [input counts](inputs_summary.json) are grouped by enumeration. The [genus of each graph](results/enumeration_graphs_genus.txt.gz) is listed by index.

## Repeating the computations

From this folder, with Python 3.8 or later, run

```sh
python -m pip install graph-genus==0.1.3
python run_genus.py
```

The [graph-genus package](https://github.com/SanderGi/Genus) supplies PAGE (Metzger and Ulrigg, 2026) and `multi_genus` (Brinkmann, 2022). Both compute the minimum orientable genus. The script checks each returned rotation against the input graph and computes its genus by Euler's formula. For $T_{25}$, the disjoint union of two $K_{3,3}$ graphs, it adds the genera of the two components.

The programs run concurrently. After the 74 obstructions, `multi_genus` processes one batch per enumeration. The [results folder](results/) contains both tables in CSV format and the [machine and software versions](results/run_info.json).

The input lists are included. To rebuild them from O1 and E1–E29, run `python build_inputs.py ../supplemental_data` from this folder with nauty's `labelg` on the path, or set `LABELG` to its executable.
