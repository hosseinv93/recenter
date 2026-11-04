// clang-format off
/* ----------------------------------------------------------------------
   Custom LAMMPS fix: recenter/two_clusters/z

   Syntax:

     fix ID group-ID recenter/two_clusters/z \
         xneg yneg zneg xpos ypos zpos clusterComputeID keyword values ...

   The selected negative cluster is the largest valid cluster mostly in z < zsplit.
   The selected positive cluster is the largest valid cluster mostly in z >= zsplit.

   Important safety behavior:
   - If one cluster is missing, no recentering is performed.
   - If the two sides are one merged cluster, no recentering is performed.
   - After shifting atoms, coordinates are remapped and a neighbor rebuild is forced.
------------------------------------------------------------------------- */

#include "fix_recenter_two_clusters_z.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "lattice.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "utils.h"

#include <cmath>
#include <cstring>
#include <map>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

static const char cite_fix_recenter_two_clusters_z[] =
    "fix recenter/two_clusters/z command: doi:10.1038/s42005-026-02903-1\n\n"
    "@Article{Vahid2026,\n"
    " author =  {H. Vahid and J.-U. Sommer and A. Sharma},\n"
    " title =   {Self-organization and collective migration of active condensates},\n"
    " journal = {Communications Physics},\n"
    " year =    2026,\n"
    " volume =  9,\n"
    " pages =   {327},\n"
    " doi =     {10.1038/s42005-026-02903-1}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixRecenterTwoClustersZ::FixRecenterTwoClustersZ(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_recenter_two_clusters_z);

  if (narg < 10)
    utils::missing_cmd_args(FLERR, "fix recenter/two_clusters/z", error);

  for (int d = 0; d < 3; d++) {
    negflag[d] = 1;
    posflag[d] = 1;

    target_neg_input[d] = 0.0;
    target_pos_input[d] = 0.0;

    target_neg[d] = 0.0;
    target_pos[d] = 0.0;

    shift_neg[d] = 0.0;
    shift_pos[d] = 0.0;

    selected_neg_xcm[d] = 0.0;
    selected_pos_xcm[d] = 0.0;
  }

  distance = 0.0;

  selected_neg_id = 0.0;
  selected_pos_id = 0.0;

  selected_neg_count = 0.0;
  selected_pos_count = 0.0;

  selected_neg_mass = 0.0;
  selected_pos_mass = 0.0;

  cluster_compute_id = nullptr;
  cluster_compute = nullptr;

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 16;

  extscalar = 0;
  extvector = 0;

  global_freq = 1;
  dynamic_group_allow = 1;

  // Important: this fix can move atoms between MPI subdomains.
  force_reneighbor = 1;
  next_reneighbor = 0;

  // ---------------- parse negative-side target ----------------

  if (strcmp(arg[3], "NULL") == 0) negflag[0] = 0;
  else target_neg_input[0] = utils::numeric(FLERR, arg[3], false, lmp);

  if (strcmp(arg[4], "NULL") == 0) negflag[1] = 0;
  else target_neg_input[1] = utils::numeric(FLERR, arg[4], false, lmp);

  if (strcmp(arg[5], "NULL") == 0) negflag[2] = 0;
  else target_neg_input[2] = utils::numeric(FLERR, arg[5], false, lmp);

  // ---------------- parse positive-side target ----------------

  if (strcmp(arg[6], "NULL") == 0) posflag[0] = 0;
  else target_pos_input[0] = utils::numeric(FLERR, arg[6], false, lmp);

  if (strcmp(arg[7], "NULL") == 0) posflag[1] = 0;
  else target_pos_input[1] = utils::numeric(FLERR, arg[7], false, lmp);

  if (strcmp(arg[8], "NULL") == 0) posflag[2] = 0;
  else target_pos_input[2] = utils::numeric(FLERR, arg[8], false, lmp);

  if (!negflag[2] || !posflag[2])
    error->all(FLERR, "fix recenter/two_clusters/z requires numeric z targets for both clusters");

  cluster_compute_id = utils::strdup(arg[9]);

  scaleflag = BOX;
  zsplit = 0.0;
  minatoms = 1;
  nevery = 1;

  int iarg = 10;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/two_clusters/z units", error);

      if (strcmp(arg[iarg + 1], "box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg + 1], "lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg + 1], "fraction") == 0) scaleflag = FRACTION;
      else error->all(FLERR, "Unknown fix recenter/two_clusters/z units argument {}", arg[iarg + 1]);

      iarg += 2;

    } else if (strcmp(arg[iarg], "zsplit") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/two_clusters/z zsplit", error);

      zsplit = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg], "minatoms") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/two_clusters/z minatoms", error);

      minatoms = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);

      if (minatoms < 1)
        error->all(FLERR, "Illegal fix recenter/two_clusters/z minatoms value");

      iarg += 2;

    } else if (strcmp(arg[iarg], "nevery") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/two_clusters/z nevery", error);

      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);

      if (nevery < 1)
        error->all(FLERR, "Illegal fix recenter/two_clusters/z nevery value");

      iarg += 2;

    } else {
      error->all(FLERR, "Unknown fix recenter/two_clusters/z keyword {}", arg[iarg]);
    }
  }

  if (scaleflag == LATTICE) {
    double scale[3];
    scale[0] = domain->lattice->xlattice;
    scale[1] = domain->lattice->ylattice;
    scale[2] = domain->lattice->zlattice;

    for (int d = 0; d < 3; d++) {
      if (negflag[d]) target_neg_input[d] *= scale[d];
      if (posflag[d]) target_pos_input[d] *= scale[d];
    }

    zsplit *= scale[2];
  }

  if (group->count(igroup) == 0)
    error->all(FLERR, "Fix recenter/two_clusters/z group {} has no atoms", arg[1]);
}

/* ---------------------------------------------------------------------- */

FixRecenterTwoClustersZ::~FixRecenterTwoClustersZ()
{
  delete[] cluster_compute_id;
}

/* ---------------------------------------------------------------------- */

int FixRecenterTwoClustersZ::setmask()
{
  int mask = 0;

  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRecenterTwoClustersZ::init()
{
  int after = 0;
  int flag = 0;

  for (const auto &ifix : modify->get_fix_list()) {
    if (strcmp(id, ifix->id) == 0) after = 1;
    else if ((modify->get_fix_mask(ifix) & INITIAL_INTEGRATE) && after) flag = 1;
  }

  if (flag && comm->me == 0)
    error->warning(FLERR, "Fix recenter/two_clusters/z should come after all other integration fixes");

  cluster_compute = modify->get_compute_by_id(cluster_compute_id);

  if (cluster_compute == nullptr)
    error->all(FLERR, "Fix recenter/two_clusters/z compute ID {} does not exist", cluster_compute_id);

  if (cluster_compute->peratom_flag == 0)
    error->all(FLERR, "Fix recenter/two_clusters/z compute {} does not compute per-atom values", cluster_compute_id);

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRecenterTwoClustersZ::convert_targets_to_box()
{
  for (int d = 0; d < 3; d++) {
    target_neg[d] = target_neg_input[d];
    target_pos[d] = target_pos_input[d];
  }

  if (scaleflag == FRACTION) {
    double *bboxlo;
    double *bboxhi;

    if (domain->triclinic == 0) {
      bboxlo = domain->boxlo;
      bboxhi = domain->boxhi;
    } else {
      bboxlo = domain->boxlo_bound;
      bboxhi = domain->boxhi_bound;
    }

    for (int d = 0; d < 3; d++) {
      if (negflag[d])
        target_neg[d] = bboxlo[d] + target_neg_input[d] * (bboxhi[d] - bboxlo[d]);

      if (posflag[d])
        target_pos[d] = bboxlo[d] + target_pos_input[d] * (bboxhi[d] - bboxlo[d]);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixRecenterTwoClustersZ::target_image_coord(int dim, double coord, double target)
{
  int periodic = 0;
  double prd = 0.0;

  if (dim == 0) {
    periodic = domain->xperiodic;
    prd = domain->xprd;
  } else if (dim == 1) {
    periodic = domain->yperiodic;
    prd = domain->yprd;
  } else if (dim == 2) {
    periodic = domain->zperiodic;
    prd = domain->zprd;
  }

  if (!periodic || prd <= 0.0) return coord;

  double d = coord - target;
  d -= prd * std::round(d / prd);

  return target + d;
}

/* ---------------------------------------------------------------------- */

double FixRecenterTwoClustersZ::pbc_correct_shift(int dim, double target, double current)
{
  int periodic = 0;
  double prd = 0.0;

  if (dim == 0) {
    periodic = domain->xperiodic;
    prd = domain->xprd;
  } else if (dim == 1) {
    periodic = domain->yperiodic;
    prd = domain->yprd;
  } else if (dim == 2) {
    periodic = domain->zperiodic;
    prd = domain->zprd;
  }

  double s = target - current;

  if (periodic && prd > 0.0) {
    s -= prd * std::round(s / prd);
  }

  return s;
}

/* ---------------------------------------------------------------------- */

void FixRecenterTwoClustersZ::find_side_clusters()
{
  struct ClusterStats {
    double mass;
    double count;

    double count_neg_side;
    double count_pos_side;

    double sum_neg[3];
    double sum_pos[3];

    ClusterStats()
    {
      mass = 0.0;
      count = 0.0;

      count_neg_side = 0.0;
      count_pos_side = 0.0;

      for (int d = 0; d < 3; d++) {
        sum_neg[d] = 0.0;
        sum_pos[d] = 0.0;
      }
    }
  };

  if (cluster_compute->invoked_peratom != update->ntimestep) {
    cluster_compute->compute_peratom();
    cluster_compute->invoked_peratom = update->ntimestep;
  }

  double *cid = cluster_compute->vector_atom;

  if (cid == nullptr)
    error->all(FLERR, "Fix recenter/two_clusters/z compute {} has null vector_atom", cluster_compute_id);

  double **x = atom->x;
  imageint *image = atom->image;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *mass = atom->mass;
  int *type = atom->type;
  double *rmass = atom->rmass;

  std::map<double, ClusterStats> local_clusters;

  double unwrap[3];

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double clid = cid[i];
    if (clid <= 0.0) continue;

    double mi = rmass ? rmass[i] : mass[type[i]];

    domain->unmap(x[i], image[i], unwrap);

    double xneg[3];
    double xpos[3];

    for (int d = 0; d < 3; d++) {
      xneg[d] = unwrap[d];
      xpos[d] = unwrap[d];
    }

    /*
       For COM calculation, map atoms to the periodic image closest
       to the requested target of each condensate.
    */
    for (int d = 0; d < 3; d++) {
      if (negflag[d])
        xneg[d] = target_image_coord(d, x[i][d], target_neg[d]);

      if (posflag[d])
        xpos[d] = target_image_coord(d, x[i][d], target_pos[d]);
    }

    /*
       For side classification, use wrapped z coordinate.
    */
    double zwrap = x[i][2];

    if (domain->zperiodic) {
      const double zlo = domain->boxlo[2];
      const double zhi = domain->boxhi[2];
      const double lz  = domain->zprd;

      while (zwrap < zlo) zwrap += lz;
      while (zwrap >= zhi) zwrap -= lz;
    }

    ClusterStats &s = local_clusters[clid];

    s.mass += mi;
    s.count += 1.0;

    if (zwrap < zsplit) s.count_neg_side += 1.0;
    else                s.count_pos_side += 1.0;

    for (int d = 0; d < 3; d++) {
      s.sum_neg[d] += mi * xneg[d];
      s.sum_pos[d] += mi * xpos[d];
    }
  }

  int nlocal_clusters = static_cast<int>(local_clusters.size());

  std::vector<double> sendbuf;
  sendbuf.reserve(11 * nlocal_clusters);

  for (auto &kv : local_clusters) {
    sendbuf.push_back(kv.first);                 // 0 cluster ID
    sendbuf.push_back(kv.second.mass);           // 1 mass
    sendbuf.push_back(kv.second.count);          // 2 total count

    sendbuf.push_back(kv.second.count_neg_side); // 3 count in z < zsplit
    sendbuf.push_back(kv.second.count_pos_side); // 4 count in z >= zsplit

    sendbuf.push_back(kv.second.sum_neg[0]);     // 5
    sendbuf.push_back(kv.second.sum_neg[1]);     // 6
    sendbuf.push_back(kv.second.sum_neg[2]);     // 7

    sendbuf.push_back(kv.second.sum_pos[0]);     // 8
    sendbuf.push_back(kv.second.sum_pos[1]);     // 9
    sendbuf.push_back(kv.second.sum_pos[2]);     // 10
  }

  int sendcount = static_cast<int>(sendbuf.size());

  std::vector<int> recvcounts(comm->nprocs);
  MPI_Allgather(&sendcount, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, world);

  std::vector<int> displs(comm->nprocs);

  int totalcount = 0;
  for (int iproc = 0; iproc < comm->nprocs; iproc++) {
    displs[iproc] = totalcount;
    totalcount += recvcounts[iproc];
  }

  if (totalcount % 11 != 0)
    error->all(FLERR, "Fix recenter/two_clusters/z internal MPI gather size error");

  double dummy = 0.0;

  std::vector<double> recvbuf(totalcount > 0 ? totalcount : 1);

  double *sendptr = sendcount > 0 ? sendbuf.data() : &dummy;
  double *recvptr = totalcount > 0 ? recvbuf.data() : &dummy;

  MPI_Allgatherv(
    sendptr, sendcount, MPI_DOUBLE,
    recvptr, recvcounts.data(), displs.data(), MPI_DOUBLE,
    world
  );

  std::map<double, ClusterStats> global_clusters;

  for (int i = 0; i < totalcount; i += 11) {
    double clid = recvbuf[i];

    ClusterStats &s = global_clusters[clid];

    s.mass += recvbuf[i + 1];
    s.count += recvbuf[i + 2];

    s.count_neg_side += recvbuf[i + 3];
    s.count_pos_side += recvbuf[i + 4];

    s.sum_neg[0] += recvbuf[i + 5];
    s.sum_neg[1] += recvbuf[i + 6];
    s.sum_neg[2] += recvbuf[i + 7];

    s.sum_pos[0] += recvbuf[i + 8];
    s.sum_pos[1] += recvbuf[i + 9];
    s.sum_pos[2] += recvbuf[i + 10];
  }

  selected_neg_id = 0.0;
  selected_pos_id = 0.0;

  selected_neg_count = 0.0;
  selected_pos_count = 0.0;

  selected_neg_mass = 0.0;
  selected_pos_mass = 0.0;

  for (int d = 0; d < 3; d++) {
    selected_neg_xcm[d] = 0.0;
    selected_pos_xcm[d] = 0.0;
  }

  /*
     Select the largest cluster mostly located in z < zsplit.
  */
  for (auto &kv : global_clusters) {
    double clid = kv.first;
    ClusterStats &s = kv.second;

    if (s.count < minatoms) continue;
    if (s.mass <= 0.0) continue;

    if (s.count_neg_side <= s.count_pos_side) continue;

    if (
      s.count > selected_neg_count ||
      (s.count == selected_neg_count && s.mass > selected_neg_mass)
    ) {
      selected_neg_id = clid;
      selected_neg_count = s.count;
      selected_neg_mass = s.mass;

      for (int d = 0; d < 3; d++)
        selected_neg_xcm[d] = s.sum_neg[d] / s.mass;
    }
  }

  /*
     Select the largest cluster mostly located in z >= zsplit.
     Do not select the same cluster twice.
  */
  for (auto &kv : global_clusters) {
    double clid = kv.first;
    ClusterStats &s = kv.second;

    if (s.count < minatoms) continue;
    if (s.mass <= 0.0) continue;

    if (selected_neg_id > 0.0 && clid == selected_neg_id) continue;

    if (s.count_pos_side <= s.count_neg_side) continue;

    if (
      s.count > selected_pos_count ||
      (s.count == selected_pos_count && s.mass > selected_pos_mass)
    ) {
      selected_pos_id = clid;
      selected_pos_count = s.count;
      selected_pos_mass = s.mass;

      for (int d = 0; d < 3; d++)
        selected_pos_xcm[d] = s.sum_pos[d] / s.mass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRecenterTwoClustersZ::initial_integrate(int /*vflag*/)
{
  if (nevery > 1 && update->ntimestep % nevery != 0) {
    for (int d = 0; d < 3; d++) {
      shift_neg[d] = 0.0;
      shift_pos[d] = 0.0;
    }

    distance = 0.0;
    return;
  }

  convert_targets_to_box();

  find_side_clusters();

  /*
     Safety condition:
     If the two condensates have merged, or if one side cannot be found,
     do not recenter anything. This avoids huge jumps and invalid atom ownership.
  */
  if (selected_neg_id <= 0.0 || selected_pos_id <= 0.0 ||
      selected_neg_id == selected_pos_id) {

    for (int d = 0; d < 3; d++) {
      shift_neg[d] = 0.0;
      shift_pos[d] = 0.0;
    }

    distance = 0.0;
    return;
  }

  for (int d = 0; d < 3; d++) {
    shift_neg[d] = 0.0;
    shift_pos[d] = 0.0;
  }

  if (selected_neg_count >= minatoms) {
    for (int d = 0; d < 3; d++) {
      if (negflag[d])
        shift_neg[d] = pbc_correct_shift(d, target_neg[d], selected_neg_xcm[d]);
    }
  }

  if (selected_pos_count >= minatoms) {
    for (int d = 0; d < 3; d++) {
      if (posflag[d])
        shift_pos[d] = pbc_correct_shift(d, target_pos[d], selected_pos_xcm[d]);
    }
  }

  distance = sqrt(
    shift_neg[0] * shift_neg[0] +
    shift_neg[1] * shift_neg[1] +
    shift_neg[2] * shift_neg[2] +
    shift_pos[0] * shift_pos[0] +
    shift_pos[1] * shift_pos[1] +
    shift_pos[2] * shift_pos[2]
  );

  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *cid = cluster_compute->vector_atom;

  if (cid == nullptr)
    error->all(FLERR, "Fix recenter/two_clusters/z compute {} has null vector_atom", cluster_compute_id);

  int moved = 0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    if (cid[i] == selected_neg_id) {
      x[i][0] += shift_neg[0];
      x[i][1] += shift_neg[1];
      x[i][2] += shift_neg[2];

      domain->remap(x[i], image[i]);
      moved = 1;

    } else if (cid[i] == selected_pos_id) {
      x[i][0] += shift_pos[0];
      x[i][1] += shift_pos[1];
      x[i][2] += shift_pos[2];

      domain->remap(x[i], image[i]);
      moved = 1;
    }
  }

  if (moved) {
    next_reneighbor = update->ntimestep;
  }
}

/* ---------------------------------------------------------------------- */

void FixRecenterTwoClustersZ::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

double FixRecenterTwoClustersZ::compute_scalar()
{
  return distance;
}

/* ---------------------------------------------------------------------- */

double FixRecenterTwoClustersZ::compute_vector(int n)
{
  if (n == 0) return shift_neg[0];
  if (n == 1) return shift_neg[1];
  if (n == 2) return shift_neg[2];

  if (n == 3) return shift_pos[0];
  if (n == 4) return shift_pos[1];
  if (n == 5) return shift_pos[2];

  if (n == 6) return selected_neg_xcm[0];
  if (n == 7) return selected_neg_xcm[1];
  if (n == 8) return selected_neg_xcm[2];

  if (n == 9) return selected_pos_xcm[0];
  if (n == 10) return selected_pos_xcm[1];
  if (n == 11) return selected_pos_xcm[2];

  if (n == 12) return selected_neg_count;
  if (n == 13) return selected_pos_count;

  if (n == 14) return selected_neg_id;
  if (n == 15) return selected_pos_id;

  return 0.0;
}
