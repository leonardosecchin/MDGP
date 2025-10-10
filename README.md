# MDGP

This is an implementation of the multistart strategy to solve the Molecular
Distance Geometry Problem (MDGP) with interval data, as described in

[Secchin, da Rocha, da Rosa, Liberti, Lavor. A fast heuristic for the molecular distance geometry problem. 2025]()

## Installation

`]add https://github.com/leonardosecchin/MDGP.git`

## Usage

The basic usage is

`mdgp_multistart(Dij, D, P, atoms, angles)`

- `Dij` is the `nd x 2` matrix of the indices `i,j` of the distances, where `nd`
is the number of distances.
- `D` is the `nd x 2` matrix of corresponding distances intervals
`[lower d_ij, upper d_ij]`. Entries of `Dij` are `Int64`, while those of `D`
are `Float64`.
- `P` is a `nv x 4` matrix with `Int64` entries, where `nv` is the number of
atoms. `P[i,1:3]` are the indices of predecessors of atom `i` in descending
order, while `P[i,4]` contains `1`, `-1` or `0` indicating the
"side" that atom `i` is located with respect to the plane of its predecessors
(quirality); when `P[i,4] = 0`, both sides are accepted.
- `atoms` is the `nv` vector of `String` containing the name of each atom (for
example, "H1", "N", "CA", "HA", "C" and so on).
- `torsions` is a `nv x 2` matrix of `Float64`'s, whose i-th row `[w, delta]`
represents the torsion angle interval in the format `[w-delta, w+delta]` for
atom `i`, in degrees. `w` must be non-negative (it sign is given by `P`).

For more details, run `?mdgp_multistart`.

### Changing parameters

You can change the algorithm parameters described in the reference paper. For
details, run `?mdgp_multistart`.

### Generating instances from the Protein Data Bank (PDB)

You can generate instances from the [PDB](https://www.rcsb.org/) using the
parser developed by Wagner da Rocha, written in Python. This parser returns the
distance and predecessor matrices, as well as the reference solution file. It is
included in the folder `scripts/MDGP_multistart`. For more details and updates,
please see [Wagner's Github page](https://github.com/wdarocha), in particular
[this link](https://github.com/wdarocha/BP_Algorithms_for_iDDGP).

The parser will generate three files, starting with `I_`, `T_` and `X_`. The
first is the distance file (`Dfile`), the second contains the predecessors
(`Pfile`) of each atom and the last the reference solution from PDB (`Xfile`).
You can read them into Julia using

`X, Dij, D, P, residues, atoms, torsions = mdgp_read("path to I_ file", "path to T_ file"; Xfile = "path to X_ file")`

For more details, run `?mdgp_read`.

## Scripts

Scripts for reproducing the tests reported in the references can be found in the
`scripts` folder. The required packages for these scripts are not necessarily
listed in the MDGP package dependencies. Therefore, you must install them
manually. Also, you must install `python3` and any required package for the PDB
parser, such as `MDAnalysis`.

## Funding

This research was supported by the São Paulo Research Foundation (FAPESP) (grant
2024/12967-8) and the National Council for Scientific and Technological
Development (CNPq) (grant 302520/2025-2), Brazil.

## How to cite

If you use this code in your publications, please cite us. For now, you can cite
the preprint:

[Secchin, da Rocha, da Rosa, Liberti, Lavor. A fast heuristic for the molecular distance geometry problem. 2025]()
