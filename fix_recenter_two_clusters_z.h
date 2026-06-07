/* -*- c++ -*- ----------------------------------------------------------
   Custom LAMMPS fix: recenter/two_clusters/z

   Recenter one cluster in z < zsplit and one cluster in z >= zsplit.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(recenter/two_clusters/z,FixRecenterTwoClustersZ);
// clang-format on
#else

#ifndef LMP_FIX_RECENTER_TWO_CLUSTERS_Z_H
#define LMP_FIX_RECENTER_TWO_CLUSTERS_Z_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRecenterTwoClustersZ : public Fix {
 public:
  FixRecenterTwoClustersZ(class LAMMPS *, int, char **);
  ~FixRecenterTwoClustersZ() override;

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void initial_integrate_respa(int, int, int) override;

  double compute_scalar() override;
  double compute_vector(int) override;

 protected:
  enum { BOX, LATTICE, FRACTION };

  int scaleflag;

  int negflag[3];
  int posflag[3];

  double target_neg_input[3];
  double target_pos_input[3];

  double target_neg[3];
  double target_pos[3];

  double zsplit;

  int nlevels_respa;

  double distance;

  double shift_neg[3];
  double shift_pos[3];

  char *cluster_compute_id;
  class Compute *cluster_compute;

  int minatoms;
  int nevery;

  double selected_neg_id;
  double selected_pos_id;

  double selected_neg_count;
  double selected_pos_count;

  double selected_neg_mass;
  double selected_pos_mass;

  double selected_neg_xcm[3];
  double selected_pos_xcm[3];

  void find_side_clusters();

  double target_image_coord(int dim, double coord, double target);
  double pbc_correct_shift(int dim, double target, double current);
  void convert_targets_to_box();
};

}    // namespace LAMMPS_NS

#endif
#endif
