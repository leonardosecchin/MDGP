# MDGP

This is an implementation of the multistart strategy to solve the Molecular Distance Geometry
Problem (MDGP) with interval data, as described in

[Secchin, da Rocha, da Rosa, Liberti, Lavor. A fast heuristic for the molecular distance geometry problem. 2025]()

## Installation

`]add https://github.com/leonardosecchin/MDGP.git`

## Use

### Input file format



### Changing parameters

You can change the algorithm parameters described in the reference paper. For details,
run `?mdgp_multistart`.


### Generating instances from the Protein Data Bank (PDB)

You can generate instances from the [PDB](https://www.rcsb.org/) using the parser developed by Wagner
da Rocha, written in Python. This parser returns the distance and predecessor matrices, as well as the
reference solution file. Please see [Wagner's Github page](https://github.com/wdarocha) for details.

The parser will generate three files, starting with `I_`, `T_` and `X_`. The first is the distance
file (`Dfile`), the second contains the predecessors (`Pfile`) of each atom and the last the reference
solution from PDB (`Xfile`). You can read them into Julia using

`X, Dij, D, P, residues, atoms, torsions = mdgp_read("path to I_ file", "path to T_ file"; Xfile = "path to X_ file")`

For more details, run `?mdgp_read`.


## Scripts

Scripts for reproducing the tests reported in the references can be found in the `scripts` folder.
The required packages for these scripts are not necessarily listed in the MDGP package dependencies.
Therefore, you must install them manually. Also, you must install `python3` and any required package
for the PDB parser, such as `MDAnalysis`.


## Funding

This research was supported by the São Paulo Research Foundation (FAPESP) (grant 2024/12967-8) and the
National Council for Scientific and Technological Development (CNPq) (grant 302520/2025-2), Brazil.


## How to cite

If you use this code in your publications, please cite us. For now, you can cite the preprint:

[Secchin, da Rocha, da Rosa, Liberti, Lavor. A fast heuristic for the molecular distance geometry problem. 2025]()
