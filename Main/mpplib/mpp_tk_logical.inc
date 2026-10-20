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
! mpp_tk_linear.inc templates for logical.

#include "mpp_tk_clear.inc"

#define TK_LOGICAL

#define TK_TYPE                                logical
#define TK_MPI                                 mpi_logical
#define TK_NAME                                'logical'
#define TK_VEC1                                lvector1
#define TK_VEC2                                lvector2
#define TK_SUBGRID                             lsubgrid
#define TK_GLOBAL_SUBGRID                      global_lsubgrid

! mpp_tk_collectives.inc
#define BCAST_TK                               bcast_logical
#define BCAST_ARR_TK                           bcast_arr_logical

! mpp_tk_distribute.inc
#define TK_2D_DO_DISTRIBUTE                    logical_2d_do_distribute
#define TK_2D_DISTRIBUTE                       logical_2d_distribute
#define TK_3D_DISTRIBUTE                       logical_3d_distribute
#define TK_4D_DISTRIBUTE                       logical_4d_distribute
#define TK_2D_DO_SUB_DISTRIBUTE                logical_2d_do_sub_distribute
#define TK_2D_SUB_DISTRIBUTE                   logical_2d_sub_distribute
#define TK_2D_DO_COLLECT                       logical_2d_do_collect
#define TK_2D_COLLECT                          logical_2d_collect
#define TK_2D_DO_SUB_COLLECT                   logical_2d_do_sub_collect
#define TK_2D_SUB_COLLECT                      logical_2d_sub_collect
#define TK_3D_SUB_COLLECT                      logical_3d_sub_collect

! mpp_tk_linear.inc
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID cartesian_to_linear_logical_subgrid_subgrid
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID linear_to_cartesian_logical_subgrid_subgrid
#define CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID    cartesian_to_linear_logical_grid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID    global_to_linear_logical_subgrid_subgrid
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID    linear_to_global_logical_subgrid_subgrid
#define GLOBAL_TO_LINEAR_TK_GRID_SUBGRID       global_to_linear_logical_grid_subgrid
#define MYPACK_TK_GRID                         mypack_logical_grid
#define MYPACK_TK_SUBGRID                      mypack_logical_subgrid
#define MYUNPACK_TK_GRID                       myunpack_logical_grid
#define MYUNPACK_TK_SUBGRID                    myunpack_logical_subgrid
#define MYPACK_GLOBAL_TK_GRID                  mypack_global_logical_grid
#define MYPACK_GLOBAL_TK_SUBGRID               mypack_global_logical_subgrid
#define MYUNPACK_GLOBAL_TK_GRID                myunpack_global_logical_grid
#define MYUNPACK_GLOBAL_TK_SUBGRID             myunpack_global_logical_subgrid

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
