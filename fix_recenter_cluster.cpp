// clang-format off
/* ----------------------------------------------------------------------
   Custom LAMMPS fix: recenter/cluster

   This fix recenters the largest cluster identified by compute cluster/atom.

   Syntax:

     fix ID group-ID recenter/cluster x y z clusterComputeID keyword values ...

   x,y,z can be numbers, NULL, or INIT.

   Required:
     clusterComputeID = ID of a compute cluster/atom

   Optional keywords:
     units box|lattice|fraction
     shift cluster|all|group-ID
     minatoms N
     nevery N

   Examples:

     compute myclus all cluster/atom 1.5
     fix keep all recenter/cluster NULL NULL -10.0 myclus units box shift all minatoms 300

   Output:
     f_ID      = magnitude of applied shift
     f_ID[1]   = shift_x
     f_ID[2]   = shift_y
     f_ID[3]   = shift_z
     f_ID[4]   = selected cluster COM_x before shifting
     f_ID[5]   = selected cluster COM_y before shifting
     f_ID[6]   = selected cluster COM_z before shifting
     f_ID[7]   = selected cluster atom count

   Safety updates:
     - Handles zero local clusters safely in MPI gather.
     - Checks for null compute vector.
     - Remaps shifted atoms back into the periodic box.
     - Forces neighbor rebuild after shifting atoms.
------------------------------------------------------------------------- */

#include "fix_recenter_cluster.h"

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

static const char cite_fix_recenter_cluster[] =
    "fix recenter/cluster command: doi:10.1038/s42005-026-02903-1\n\n"
    "@Article{ActiveCondensates2026,\n"
    " title =   {Self-organization and collective migration of active condensates},\n"
    " journal = {Communications Physics},\n"
    " year =    2026,\n"
    " volume =  9,\n"
    " pages =   {327},\n"
    " doi =     {10.1038/s42005-026-02903-1}\n"
    "}\n\n";

/* ---------------------------------------------------------------------- */

FixRecenterCluster::FixRecenterCluster(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (lmp->citeme) lmp->citeme->add(cite_fix_recenter_cluster);

  if (narg < 7) utils::missing_cmd_args(FLERR, "fix recenter/cluster", error);

  xcom = ycom = zcom = 0.0;
  xflag = yflag = zflag = 1;
  xinitflag = yinitflag = zinitflag = 0;

  xinit = yinit = zinit = 0.0;

  shift[0] = shift[1] = shift[2] = 0.0;
  distance = 0.0;

  selected_cluster_id = 0.0;
  selected_cluster_count = 0.0;
  selected_cluster_mass = 0.0;
  selected_xcm[0] = selected_xcm[1] = selected_xcm[2] = 0.0;

  cluster_compute_id = nullptr;
  cluster_compute = nullptr;

  scalar_flag = 1;
  vector_flag = 1;

  // 1-3: shift vector
  // 4-6: selected cluster COM before shift
  // 7: selected cluster atom count
  size_vector = 7;

  // Do not let thermo normalize these values by number of atoms.
  extscalar = 0;
  extvector = 0;

  global_freq = 1;
  dynamic_group_allow = 1;

  // Important: this fix can move atoms between MPI subdomains.
  force_reneighbor = 1;
  next_reneighbor = 0;

  if (strcmp(arg[3], "NULL") == 0) xflag = 0;
  else if (strcmp(arg[3], "INIT") == 0) xinitflag = 1;
  else xcom = utils::numeric(FLERR, arg[3], false, lmp);

  if (strcmp(arg[4], "NULL") == 0) yflag = 0;
  else if (strcmp(arg[4], "INIT") == 0) yinitflag = 1;
  else ycom = utils::numeric(FLERR, arg[4], false, lmp);

  if (strcmp(arg[5], "NULL") == 0) zflag = 0;
  else if (strcmp(arg[5], "INIT") == 0) zinitflag = 1;
  else zcom = utils::numeric(FLERR, arg[5], false, lmp);

  cluster_compute_id = utils::strdup(arg[6]);

  // Default to box units.
  scaleflag = BOX;

  // Default: shift only selected largest cluster.
  // Recommended for some analyses: use "shift all" in the input.
  shiftmode = SHIFT_CLUSTER;

  group2bit = groupbit;
  minatoms = 1;
  nevery = 1;

  int iarg = 7;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/cluster units", error);

      if (strcmp(arg[iarg + 1], "box") == 0) scaleflag = BOX;
      else if (strcmp(arg[iarg + 1], "lattice") == 0) scaleflag = LATTICE;
      else if (strcmp(arg[iarg + 1], "fraction") == 0) scaleflag = FRACTION;
      else error->all(FLERR, "Unknown fix recenter/cluster units argument {}", arg[iarg + 1]);

      iarg += 2;

    } else if (strcmp(arg[iarg], "shift") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/cluster shift", error);

      if (strcmp(arg[iarg + 1], "cluster") == 0) {
        shiftmode = SHIFT_CLUSTER;
      } else if (strcmp(arg[iarg + 1], "all") == 0) {
        shiftmode = SHIFT_ALL;
      } else {
        shiftmode = SHIFT_GROUP;
        group2bit = group->get_bitmask_by_id(FLERR, arg[iarg + 1], "fix recenter/cluster");
      }

      iarg += 2;

    } else if (strcmp(arg[iarg], "minatoms") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/cluster minatoms", error);

      minatoms = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);

      if (minatoms < 1)
        error->all(FLERR, "Illegal fix recenter/cluster minatoms value");

      iarg += 2;

    } else if (strcmp(arg[iarg], "nevery") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix recenter/cluster nevery", error);

      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);

      if (nevery < 1)
        error->all(FLERR, "Illegal fix recenter/cluster nevery value");

      iarg += 2;

    } else {
      error->all(FLERR, "Unknown fix recenter/cluster keyword {}", arg[iarg]);
    }
  }

  double xscale, yscale, zscale;

  if (scaleflag == LATTICE) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else {
    xscale = yscale = zscale = 1.0;
  }

  xcom *= xscale;
  ycom *= yscale;
  zcom *= zscale;

  if (group->count(igroup) == 0)
    error->all(FLERR, "Fix recenter/cluster group {} has no atoms", arg[1]);
}

/* ---------------------------------------------------------------------- */

FixRecenterCluster::~FixRecenterCluster()
{
  delete[] cluster_compute_id;
}

/* ---------------------------------------------------------------------- */

int FixRecenterCluster::setmask()
{
  int mask = 0;

  // Same stage as standard fix recenter.
  // This fix should be defined after the integration fix.
  mask |= INITIAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;

  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRecenterCluster::init()
{
  int after = 0;
  int flag = 0;

  for (const auto &ifix : modify->get_fix_list()) {
    if (strcmp(id, ifix->id) == 0) after = 1;
    else if ((modify->get_fix_mask(ifix) & INITIAL_INTEGRATE) && after) flag = 1;
  }

  if (flag && comm->me == 0)
    error->warning(FLERR, "Fix recenter/cluster should come after all other integration fixes");

  cluster_compute = modify->get_compute_by_id(cluster_compute_id);

  if (cluster_compute == nullptr)
    error->all(FLERR, "Fix recenter/cluster compute ID {} does not exist", cluster_compute_id);

  if (cluster_compute->peratom_flag == 0)
    error->all(FLERR, "Fix recenter/cluster compute {} does not compute per-atom values", cluster_compute_id);

  if (xinitflag || yinitflag || zinitflag) {
    find_largest_cluster();

    if (selected_cluster_count < minatoms)
      error->all(FLERR, "Fix recenter/cluster found no valid initial cluster");

    double xcm[3];
    compute_cluster_com(xcm, nullptr);

    xinit = xcm[0];
    yinit = xcm[1];
    zinit = xcm[2];
  }

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixRecenterCluster::find_largest_cluster()
{
  if (cluster_compute->invoked_peratom != update->ntimestep) {
    cluster_compute->compute_peratom();
    cluster_compute->invoked_peratom = update->ntimestep;
  }

  double *cid = cluster_compute->vector_atom;

  if (cid == nullptr)
    error->all(FLERR, "Fix recenter/cluster compute {} has null vector_atom", cluster_compute_id);

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *mass = atom->mass;
  int *type = atom->type;
  double *rmass = atom->rmass;

  std::map<double, std::pair<double, double>> local_clusters;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    double clid = cid[i];

    if (clid <= 0.0) continue;

    double mi = rmass ? rmass[i] : mass[type[i]];

    local_clusters[clid].first += mi;
    local_clusters[clid].second += 1.0;
  }

  int nlocal_clusters = static_cast<int>(local_clusters.size());

  std::vector<double> sendbuf;
  sendbuf.reserve(3 * nlocal_clusters);

  for (auto &kv : local_clusters) {
    sendbuf.push_back(kv.first);
    sendbuf.push_back(kv.second.first);
    sendbuf.push_back(kv.second.second);
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

  if (totalcount % 3 != 0)
    error->all(FLERR, "Fix recenter/cluster internal MPI gather size error");

  double dummy = 0.0;

  std::vector<double> recvbuf(totalcount > 0 ? totalcount : 1);

  double *sendptr = sendcount > 0 ? sendbuf.data() : &dummy;
  double *recvptr = totalcount > 0 ? recvbuf.data() : &dummy;

  MPI_Allgatherv(
    sendptr, sendcount, MPI_DOUBLE,
    recvptr, recvcounts.data(), displs.data(), MPI_DOUBLE,
    world
  );

  std::map<double, std::pair<double, double>> global_clusters;

  for (int i = 0; i < totalcount; i += 3) {
    double clid = recvbuf[i];
    double msum = recvbuf[i + 1];
    double csum = recvbuf[i + 2];

    global_clusters[clid].first += msum;
    global_clusters[clid].second += csum;
  }

  selected_cluster_id = 0.0;
  selected_cluster_mass = 0.0;
  selected_cluster_count = 0.0;

  for (auto &kv : global_clusters) {
    double clid = kv.first;
    double msum = kv.second.first;
    double csum = kv.second.second;

    if (csum < minatoms) continue;

    if (
      csum > selected_cluster_count ||
      (csum == selected_cluster_count && msum > selected_cluster_mass)
    ) {
      selected_cluster_id = clid;
      selected_cluster_mass = msum;
      selected_cluster_count = csum;
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixRecenterCluster::target_image_coord(int dim, double coord, double target)
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

  // Map coordinate to the periodic image closest to target.
  d -= prd * std::round(d / prd);

  return target + d;
}

/* ---------------------------------------------------------------------- */

double FixRecenterCluster::pbc_correct_shift(int dim, double target, double current)
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

void FixRecenterCluster::compute_cluster_com(double *xcm, const double *target)
{
  if (selected_cluster_id <= 0.0 || selected_cluster_mass <= 0.0) {
    xcm[0] = xcm[1] = xcm[2] = 0.0;
    return;
  }

  if (cluster_compute->invoked_peratom != update->ntimestep) {
    cluster_compute->compute_peratom();
    cluster_compute->invoked_peratom = update->ntimestep;
  }

  double *cid = cluster_compute->vector_atom;

  if (cid == nullptr)
    error->all(FLERR, "Fix recenter/cluster compute {} has null vector_atom", cluster_compute_id);

  double **x = atom->x;
  imageint *image = atom->image;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *mass = atom->mass;
  int *type = atom->type;
  double *rmass = atom->rmass;

  double local[4];
  local[0] = local[1] = local[2] = local[3] = 0.0;

  double unwrap[3];

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (cid[i] != selected_cluster_id) continue;

    double mi = rmass ? rmass[i] : mass[type[i]];

    domain->unmap(x[i], image[i], unwrap);

    double xi = unwrap[0];
    double yi = unwrap[1];
    double zi = unwrap[2];

    // For recentered directions, use wrapped coordinate mapped to the
    // periodic image closest to the target. This prevents a droplet/slab
    // across PBC from being interpreted as far away from the target.
    if (target != nullptr) {
      if (xflag) xi = target_image_coord(0, x[i][0], target[0]);
      if (yflag) yi = target_image_coord(1, x[i][1], target[1]);
      if (zflag) zi = target_image_coord(2, x[i][2], target[2]);
    }

    local[0] += mi * xi;
    local[1] += mi * yi;
    local[2] += mi * zi;
    local[3] += mi;
  }

  double global[4];
  MPI_Allreduce(local, global, 4, MPI_DOUBLE, MPI_SUM, world);

  if (global[3] <= 0.0)
    error->all(FLERR, "Fix recenter/cluster selected cluster has zero mass");

  xcm[0] = global[0] / global[3];
  xcm[1] = global[1] / global[3];
  xcm[2] = global[2] / global[3];
}

/* ---------------------------------------------------------------------- */

void FixRecenterCluster::initial_integrate(int /*vflag*/)
{
  if (nevery > 1 && update->ntimestep % nevery != 0) {
    shift[0] = shift[1] = shift[2] = 0.0;
    distance = 0.0;
    return;
  }

  double xtarget, ytarget, ztarget;
  double *bboxlo = nullptr;
  double *bboxhi = nullptr;

  if (scaleflag == FRACTION) {
    if (domain->triclinic == 0) {
      bboxlo = domain->boxlo;
      bboxhi = domain->boxhi;
    } else {
      bboxlo = domain->boxlo_bound;
      bboxhi = domain->boxhi_bound;
    }
  }

  if (xinitflag) xtarget = xinit;
  else if (scaleflag == FRACTION) xtarget = bboxlo[0] + xcom * (bboxhi[0] - bboxlo[0]);
  else xtarget = xcom;

  if (yinitflag) ytarget = yinit;
  else if (scaleflag == FRACTION) ytarget = bboxlo[1] + ycom * (bboxhi[1] - bboxlo[1]);
  else ytarget = ycom;

  if (zinitflag) ztarget = zinit;
  else if (scaleflag == FRACTION) ztarget = bboxlo[2] + zcom * (bboxhi[2] - bboxlo[2]);
  else ztarget = zcom;

  find_largest_cluster();

  if (selected_cluster_count < minatoms || selected_cluster_id <= 0.0) {
    shift[0] = shift[1] = shift[2] = 0.0;
    distance = 0.0;

    selected_xcm[0] = selected_xcm[1] = selected_xcm[2] = 0.0;

    return;
  }

  double target[3];
  target[0] = xtarget;
  target[1] = ytarget;
  target[2] = ztarget;

  double xcm[3];
  compute_cluster_com(xcm, target);

  selected_xcm[0] = xcm[0];
  selected_xcm[1] = xcm[1];
  selected_xcm[2] = xcm[2];

  shift[0] = xflag ? pbc_correct_shift(0, xtarget, xcm[0]) : 0.0;
  shift[1] = yflag ? pbc_correct_shift(1, ytarget, xcm[1]) : 0.0;
  shift[2] = zflag ? pbc_correct_shift(2, ztarget, xcm[2]) : 0.0;

  distance = sqrt(
    shift[0] * shift[0] +
    shift[1] * shift[1] +
    shift[2] * shift[2]
  );

  double **x = atom->x;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *cid = cluster_compute->vector_atom;

  if (cid == nullptr)
    error->all(FLERR, "Fix recenter/cluster compute {} has null vector_atom", cluster_compute_id);

  int moved = 0;

  for (int i = 0; i < nlocal; i++) {
    int do_shift = 0;

    if (shiftmode == SHIFT_CLUSTER) {
      if ((mask[i] & groupbit) && cid[i] == selected_cluster_id) do_shift = 1;
    } else if (shiftmode == SHIFT_ALL) {
      do_shift = 1;
    } else if (shiftmode == SHIFT_GROUP) {
      if (mask[i] & group2bit) do_shift = 1;
    }

    if (do_shift) {
      x[i][0] += shift[0];
      x[i][1] += shift[1];
      x[i][2] += shift[2];

      domain->remap(x[i], image[i]);
      moved = 1;
    }
  }

  if (moved) {
    next_reneighbor = update->ntimestep;
  }
}

/* ---------------------------------------------------------------------- */

void FixRecenterCluster::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) initial_integrate(vflag);
}

/* ---------------------------------------------------------------------- */

double FixRecenterCluster::compute_scalar()
{
  return distance;
}

/* ---------------------------------------------------------------------- */

double FixRecenterCluster::compute_vector(int n)
{
  if (n == 0) return shift[0];
  if (n == 1) return shift[1];
  if (n == 2) return shift[2];

  if (n == 3) return selected_xcm[0];
  if (n == 4) return selected_xcm[1];
  if (n == 5) return selected_xcm[2];

  if (n == 6) return selected_cluster_count;

  return 0.0;
}
