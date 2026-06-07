/* -*- c++ -*- ----------------------------------------------------------
   Custom LAMMPS fix: recenter/cluster

   Recenter the largest cluster, identified from a compute cluster/atom.

   Safety updates:
   - Handles zero local clusters safely in MPI gather.
   - Checks for null compute vector.
   - Remaps shifted atoms back into the periodic box.
   - Forces neighbor rebuild after shifting atoms.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(recenter/cluster,FixRecenterCluster);
// clang-format on
#else

#ifndef LMP_FIX_RECENTER_CLUSTER_H
#define LMP_FIX_RECENTER_CLUSTER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRecenterCluster : public Fix {
 public:
  FixRecenterCluster(class LAMMPS *, int, char **);
  ~FixRecenterCluster() override;

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void initial_integrate_respa(int, int, int) override;

  double compute_scalar() override;
  double compute_vector(int) override;

 protected:
  enum { BOX, LATTICE, FRACTION };
  enum { SHIFT_CLUSTER, SHIFT_ALL, SHIFT_GROUP };

  int scaleflag;
  int shiftmode;
  int group2bit;

  int xflag, yflag, zflag;
  int xinitflag, yinitflag, zinitflag;

  int nlevels_respa;

  double xcom, ycom, zcom;
  double xinit, yinit, zinit;

  double distance;
  double shift[3];

  char *cluster_compute_id;
  class Compute *cluster_compute;

  int minatoms;
  int nevery;

  double selected_cluster_id;
  double selected_cluster_count;
  double selected_cluster_mass;

  // COM of selected cluster before applying the recentering shift.
  double selected_xcm[3];

  void find_largest_cluster();

  double target_image_coord(int dim, double coord, double target);
  double pbc_correct_shift(int dim, double target, double current);

  void compute_cluster_com(double *xcm, const double *target);
};

}    // namespace LAMMPS_NS

#endif
#endif
