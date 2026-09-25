# Cluster recentering fixes for LAMMPS

Two custom LAMMPS fixes that keep condensates (droplets, slabs, aggregates) at a chosen position in the box. They work like the built-in [`fix recenter`](https://docs.lammps.org/fix_recenter.html), but the reference is the **largest cluster** found by [`compute cluster/atom`](https://docs.lammps.org/compute_cluster_atom.html), not the whole group. Dilute atoms that leave a droplet therefore do not drag its reference point around.

| Fix style | Files | Purpose |
|---|---|---|
| `recenter/cluster` | `fix_recenter_cluster.cpp/.h` | Recenter the single largest cluster |
| `recenter/two_clusters/z` | `fix_recenter_two_clusters_z.cpp/.h` | Recenter two clusters, one below and one above a z plane, independently |

Both fixes shift coordinates only. They change no velocities or forces and add no energy.

## Installation

The fixes use only core LAMMPS classes. Copy them into `src/` and rebuild:

```bash
cp fix_recenter_*.cpp fix_recenter_*.h /path/to/lammps/src/

cd /path/to/lammps
cmake -S cmake -B build     # re-run so CMake picks up the new files
cmake --build build -j
```

`compute cluster/atom` is part of core LAMMPS.

## `fix recenter/cluster`

```
fix ID group-ID recenter/cluster x y z clusterComputeID keyword value ...
```

- `x`, `y`, `z`: target position of the cluster's centre of mass. Each is a number, `NULL` (leave this direction alone) or `INIT` (the cluster's centre of mass when the run is set up).
- `clusterComputeID`: ID of a `compute cluster/atom` (any compute that gives a per-atom cluster ID works).

| Keyword | Values | Default | Meaning |
|---|---|---|---|
| `units` | `box`, `lattice`, `fraction` | `box` | Units of `x y z`. `fraction` means fractions of the box, re-evaluated every step |
| `shift` | `cluster`, `all`, or a group ID | `cluster` | Which atoms are moved: only the selected cluster, every atom, or the atoms of that group |
| `minatoms` | N ≥ 1 | 1 | Ignore clusters with fewer than N atoms. If no cluster qualifies, nothing is moved that step |
| `nevery` | N ≥ 1 | 1 | Recenter only on steps that are multiples of N |

The selected cluster is the one with the most atoms among the atoms in `group-ID`, with ties broken by mass. Its centre of mass is mass-weighted. Along every recentred periodic direction, each atom is taken at the periodic image closest to the target, so a droplet or slab that straddles a boundary is handled correctly. The shift uses the minimum-image convention, and shifted atoms are wrapped back into the box. The fix forces a neighbor-list rebuild on steps where it moves atoms.

With `shift cluster`, only the cluster moves relative to its surroundings, which can create overlaps with nearby atoms. Use `shift all` to translate the whole system rigidly, like `fix recenter`.

**Output** (global, intensive):

| Quantity | Meaning |
|---|---|
| `f_ID` | Length of the applied shift |
| `f_ID[1..3]` | Shift vector |
| `f_ID[4..6]` | Cluster centre of mass before shifting |
| `f_ID[7]` | Number of atoms in the selected cluster |

All values are 0 on steps without recentering (`nevery`, or no valid cluster).

**Example**: keep a droplet's centre at z = −10 by translating the whole system:

```lammps
compute clus all cluster/atom 1.5
fix keep all recenter/cluster NULL NULL -10.0 clus units box shift all minatoms 300
thermo_style custom step temp pe f_keep f_keep[3] f_keep[7]
```

## `fix recenter/two_clusters/z`

```
fix ID group-ID recenter/two_clusters/z xneg yneg zneg xpos ypos zpos clusterComputeID keyword value ...
```

- `xneg yneg zneg`: target centre of mass for the cluster below the split plane.
- `xpos ypos zpos`: target for the cluster above it.
- The x and y targets may be `NULL`. **Both z targets must be numbers.** `INIT` is not supported.

| Keyword | Values | Default | Meaning |
|---|---|---|---|
| `units` | `box`, `lattice`, `fraction` | `box` | Units of the six targets |
| `zsplit` | z | 0.0 | Plane that separates the two sides |
| `minatoms` | N ≥ 1 | 1 | Minimum cluster size |
| `nevery` | N ≥ 1 | 1 | Recenter every N steps |

The "negative" cluster is the largest cluster (≥ `minatoms`) with more of its atoms at `z < zsplit` than above, judged from wrapped coordinates. The "positive" cluster is the largest with more atoms at `z ≥ zsplit`, excluding the one already chosen. Each cluster is shifted to its own target. Only the atoms of the two selected clusters (in `group-ID`) are moved.

**Safety rule:** if either cluster is missing, or if both sides belong to the same merged cluster, the fix moves nothing on that step.

`zsplit` is scaled for `units lattice` but **not** for `units fraction`, where it stays in box units.

**Output** (global, intensive):

| Quantity | Meaning |
|---|---|
| `f_ID` | Combined length of both shift vectors |
| `f_ID[1..3]`, `f_ID[4..6]` | Shift of the negative and of the positive cluster |
| `f_ID[7..9]`, `f_ID[10..12]` | Centre of mass of each cluster before shifting |
| `f_ID[13]`, `f_ID[14]` | Atom count of each cluster |
| `f_ID[15]`, `f_ID[16]` | Cluster ID of each cluster (0 if none) |

**Example**: two condensates held at z = −30 and z = +30 on either side of z = 0:

```lammps
compute clus all cluster/atom 1.5
fix hold all recenter/two_clusters/z NULL NULL -30.0 NULL NULL 30.0 clus zsplit 0.0 minatoms 200
```

## Usage notes

- Define these fixes **after** the time-integration fixes, as with `fix recenter`. Otherwise the fix prints a warning.
- `compute cluster/atom` needs a neighbor list whose cutoff covers the clustering distance. A pair cutoff of at least that distance is enough.
- Recentering is done at the start of each step (`initial_integrate`), and rRESPA runs use the outermost level.
- Instantaneous shifts make unwrapped trajectories discontinuous. Account for this, for example by subtracting the reported shifts, when computing diffusion or other displacement-based quantities.

## License

GNU General Public License v2 or later, the same as LAMMPS.
