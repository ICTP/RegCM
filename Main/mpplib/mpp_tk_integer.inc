!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Instantiate the mpp_tk_collectives.inc, mpp_tk_distribute.inc and
! mpp_tk_linear.inc templates for integer(ik4).

#include "mpp_tk_clear.inc"

#define TK_TYPE                                integer(ik4)
#define TK_MPI                                 mpi_integer4
#define TK_NAME                                'integer'
#define TK_VEC1                                i4vector1
#define TK_VEC2                                i4vector2
#define TK_SUBGRID                             i4subgrid
#define TK_GLOBAL_SUBGRID                      global_i4subgrid

! mpp_tk_collectives.inc
#define BCAST_TK                               bcast_int4
#define BCAST_ARR_TK                           bcast_arr_int4
#define SUMALL_TK                              sumall_int4
#define MAXALL_TK                              maxall_integer4
#define MINALL_TK                              minall_integer4

! mpp_tk_distribute.inc
#define TK_2D_DO_DISTRIBUTE                    integer4_2d_do_distribute
#define TK_2D_DISTRIBUTE                       integer_2d_distribute
#define TK_3D_DISTRIBUTE                       integer_3d_distribute
#define TK_4D_DISTRIBUTE                       integer_4d_distribute
#define TK_2D_DO_SUB_DISTRIBUTE                integer4_2d_do_sub_distribute
#define TK_2D_SUB_DISTRIBUTE                   integer_2d_sub_distribute
#define TK_3D_SUB_DISTRIBUTE                   integer_3d_sub_distribute
#define TK_2D_DO_COLLECT                       integer4_2d_do_collect
#define TK_COLLECT_PLANES                      integer4_collect_planes
#define TK_2D_COLLECT                          integer_2d_collect
#define TK_3D_COLLECT                          integer_3d_collect
#define TK_4D_COLLECT                          integer_4d_collect
#define TK_2D_DO_SUB_COLLECT                   integer4_2d_do_sub_collect
#define TK_2D_SUB_COLLECT                      integer_2d_sub_collect
#define TK_3D_SUB_COLLECT                      integer_3d_sub_collect

! mpp_tk_linear.inc
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID cartesian_to_linear_integer_subgrid_subgrid
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID linear_to_cartesian_integer_subgrid_subgrid
#define CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID    cartesian_to_linear_integer_grid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID    global_to_linear_integer_subgrid_subgrid
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID    linear_to_global_integer_subgrid_subgrid
#define GLOBAL_TO_LINEAR_TK_GRID_SUBGRID       global_to_linear_integer_grid_subgrid
#define MYPACK_TK_GRID                         mypack_integer_grid
#define MYPACK_TK_SUBGRID                      mypack_integer_subgrid
#define MYUNPACK_TK_GRID                       myunpack_integer_grid
#define MYUNPACK_TK_SUBGRID                    myunpack_integer_subgrid
#define MYPACK_GLOBAL_TK_GRID                  mypack_global_integer_grid
#define MYPACK_GLOBAL_TK_SUBGRID               mypack_global_integer_subgrid
#define MYUNPACK_GLOBAL_TK_GRID                myunpack_global_integer_grid
#define MYUNPACK_GLOBAL_TK_SUBGRID             myunpack_global_integer_subgrid

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
