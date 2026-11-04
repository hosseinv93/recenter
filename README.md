# fix recenter/cluster for LAMMPS

Two extra [LAMMPS](https://www.lammps.org) fixes that keep a droplet, slab or aggregate at a fixed position in the simulation box.

They work like the built-in [`fix recenter`](https://docs.lammps.org/fix_recenter.html), with one difference. They follow only the **largest cluster** found by [`compute cluster/atom`](https://docs.lammps.org/compute_cluster_atom.html), so loose atoms in the dilute phase do not pull the reference point around.

| Fix | What it does |
|---|---|
| `recenter/cluster` | Keeps the largest cluster at a target position |
| `recenter/two_clusters/z` | Keeps two clusters, one below and one above a z plane, at two separate targets |

The fixes only shift coordinates. They do not change velocities or forces.

If you use them, please [cite our paper](#how-to-cite).

## Install

You need the LAMMPS source code and CMake. No extra packages are required.

1. Copy the four files into the LAMMPS `src` folder:

   ```bash
   cp fix_recenter_*.cpp fix_recenter_*.h /path/to/lammps/src/
   ```

2. Build LAMMPS as usual. If you already have a build folder, run the first command again so CMake sees the new files:

   ```bash
   cd /path/to/lammps
   cmake -S cmake -B build
   cmake --build build -j 4
   ```

3. Check that the fixes are included:

   ```bash
   ./build/lmp -h | grep recenter
   ```

   You should see `recenter/cluster` and `recenter/two_clusters/z` in the list.

## Quick start

Add these lines to your input script, after your integrator fix (for example `fix nvt`):

```lammps
compute clus all cluster/atom 1.5
fix keep all recenter/cluster NULL NULL 0.0 clus shift all minatoms 100
```

This finds clusters with a cutoff of 1.5. It then moves the whole system every step so that the largest cluster's centre of mass stays at z = 0. Clusters smaller than 100 atoms are ignored.

To hold two condensates at z = -30 and z = +30:

```lammps
compute clus all cluster/atom 1.5
fix hold all recenter/two_clusters/z NULL NULL -30.0 NULL NULL 30.0 clus zsplit 0.0 minatoms 100
```

## `fix recenter/cluster`

```
fix ID group-ID recenter/cluster x y z clusterComputeID [keyword value ...]
```

- `x y z` is the target position. Each value is a number, `NULL` to leave that direction alone, or `INIT` to use the cluster's starting position.
- `clusterComputeID` is the ID of a `compute cluster/atom`.

| Keyword | Values | Default | Meaning |
|---|---|---|---|
| `units` | `box`, `lattice`, `fraction` | `box` | Units of `x y z` |
| `shift` | `cluster`, `all`, or a group ID | `cluster` | Which atoms are moved |
| `minatoms` | N | 1 | Ignore clusters smaller than N atoms |
| `nevery` | N | 1 | Recenter every N steps |

Tip: use `shift all` to move the whole system together. With `shift cluster`, only the cluster moves, which can push it into nearby atoms.

Output for `thermo_style` (`f_ID`):

| Value | Meaning |
|---|---|
| `f_ID` | Size of the shift |
| `f_ID[1]` to `f_ID[3]` | Shift in x, y, z |
| `f_ID[4]` to `f_ID[6]` | Cluster centre of mass before the shift |
| `f_ID[7]` | Number of atoms in the cluster |

## `fix recenter/two_clusters/z`

```
fix ID group-ID recenter/two_clusters/z xneg yneg zneg xpos ypos zpos clusterComputeID [keyword value ...]
```

- `xneg yneg zneg` is the target for the cluster below `zsplit`.
- `xpos ypos zpos` is the target for the cluster above `zsplit`.
- x and y can be `NULL`. Both z values must be numbers.

| Keyword | Values | Default | Meaning |
|---|---|---|---|
| `units` | `box`, `lattice`, `fraction` | `box` | Units of the targets |
| `zsplit` | z | 0.0 | Plane that separates the two sides. Uses box units, even with `units fraction` |
| `minatoms` | N | 1 | Ignore clusters smaller than N atoms |
| `nevery` | N | 1 | Recenter every N steps |

Each cluster is moved to its own target. If one cluster is missing, or the two have merged into one, nothing is moved on that step.

Output for `thermo_style` (`f_ID`):

| Value | Meaning |
|---|---|
| `f_ID` | Total size of both shifts |
| `f_ID[1]` to `f_ID[3]` | Shift of the lower cluster |
| `f_ID[4]` to `f_ID[6]` | Shift of the upper cluster |
| `f_ID[7]` to `f_ID[9]` | Lower cluster centre of mass before the shift |
| `f_ID[10]` to `f_ID[12]` | Upper cluster centre of mass before the shift |
| `f_ID[13]`, `f_ID[14]` | Atom count of each cluster |
| `f_ID[15]`, `f_ID[16]` | Cluster ID of each cluster, 0 if none |

## Good to know

- Define these fixes after your integrator fix, the same as `fix recenter`.
- The pair cutoff must be at least the `cluster/atom` cutoff.
- Clusters that cross a periodic boundary are handled correctly.
- The shifts make unwrapped trajectories jump. Subtract the reported shifts before computing diffusion.

## How to cite

If you use these fixes in your work, please cite this paper:

> Self-organization and collective migration of active condensates. *Communications Physics* **9**, 327 (2026). https://doi.org/10.1038/s42005-026-02903-1

LAMMPS also prints this reference on screen and in `log.cite` when you run a simulation that uses either fix.

## License

GNU General Public License v2 or later, the same as LAMMPS. See [LICENSE](LICENSE).
