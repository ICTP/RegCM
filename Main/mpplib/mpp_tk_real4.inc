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

! Instantiate the mpp_tk_*.inc templates for real(rk4).
! It must be included BEFORE the type/kind methods!

#include "mpp_tk_clear.inc"

#define TK_REAL

#define TK_TYPE                                   real(rk4)
#define TK_KIND                                   rk4
#define TK_MPI                                    mpi_real4
#define TK_NAME                                   'real4'
#define TK_VEC1                                   r4vector1
#define TK_VEC2                                   r4vector2
#define TK_SUBGRID                                r4subgrid
#define TK_GLOBAL_SUBGRID                         global_r4subgrid
#define TK_MISSVAL                                smissval_r4

! mpp_tk_collectives.inc
#define BCAST_TK                                  bcast_real4
#define BCAST_ARR_TK                              bcast_arr_real4
#define BCAST_MATR_TK                             bcast_matr_real4
#define SUMALL_TK                                 sumall_real4
#define MAXALL_TK                                 maxall_real4
#define MEANALL_1D_TK                             meanall_1D_real4
#define MEANALL_TK                                meanall_real4
#define MINALL_TK                                 minall_real4
#define EXCHANGE_ARRAY_TK                         exchange_array_r4
#define CYCLIC_EXCHANGE_ARRAY_TK                  cyclic_exchange_array_r4
#define TK_COLUMN_REDUCE                          real4_column_reduce
#define TK_ROW_REDUCE                             real4_row_reduce

! mpp_tk_distribute.inc
#define TK_2D_DO_DISTRIBUTE                       real4_2d_do_distribute
#define TK_2D_DISTRIBUTE                          real4_2d_distribute
#define TK_3D_DISTRIBUTE                          real4_3d_distribute
#define TK_4D_DISTRIBUTE                          real4_4d_distribute
#define TK_2D_DO_SUB_DISTRIBUTE                   real4_2d_do_sub_distribute
#define TK_2D_SUB_DISTRIBUTE                      real4_2d_sub_distribute
#define TK_3D_SUB_DISTRIBUTE                      real4_3d_sub_distribute
#define TK_2D_DO_COLLECT                          real4_2d_do_collect
#define TK_COLLECT_PLANES                         real4_collect_planes
#define TK_2D_COLLECT                             real4_2d_collect
#define TK_2D_3D_COLLECT                          real4_2d_3d_collect
#define TK_3D_COLLECT                             real4_3d_collect
#define TK_3D_2D_COLLECT                          real4_3d_2d_collect
#define TK_4D_COLLECT                             real4_4d_collect
#define TK_4D_2D_COLLECT                          real4_4d_2d_collect
#define TK_2D_DO_SUB_COLLECT                      real4_2d_do_sub_collect
#define TK_2D_SUB_COLLECT                         real4_2d_sub_collect
#define TK_3D_SUB_COLLECT                         real4_3d_sub_collect
#define TK_2D_GRID_FILL_EXTEND1                   real4_2d_grid_fill_extend1
#define TK_2D_GRID_FILL_EXTEND2                   real4_2d_grid_fill_extend2

! mpp_tk_exchange.inc
#define TK_2D_EXCHANGE                            real4_2d_exchange
#define TK_3D_EXCHANGE                            real4_3d_exchange
#define TK_4D_EXCHANGE                            real4_4d_exchange
#define TK_2D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real4_2d_exchange_left_right_bottom_top
#define TK_3D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real4_3d_exchange_left_right_bottom_top
#define TK_4D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real4_4d_exchange_left_right_bottom_top
#define TK_2D_EXCHANGE_LEFT_RIGHT                 real4_2d_exchange_left_right
#define TK_3D_EXCHANGE_LEFT_RIGHT                 real4_3d_exchange_left_right
#define TK_4D_EXCHANGE_LEFT_RIGHT                 real4_4d_exchange_left_right
#define TK_2D_EXCHANGE_BOTTOM_TOP                 real4_2d_exchange_bottom_top
#define TK_3D_EXCHANGE_BOTTOM_TOP                 real4_3d_exchange_bottom_top
#define TK_4D_EXCHANGE_BOTTOM_TOP                 real4_4d_exchange_bottom_top
#define TK_2D_EXCHANGE_LEFT_BOTTOM                real4_2d_exchange_left_bottom
#define TK_3D_EXCHANGE_LEFT_BOTTOM                real4_3d_exchange_left_bottom
#define TK_4D_EXCHANGE_LEFT_BOTTOM                real4_4d_exchange_left_bottom
#define TK_2D_EXCHANGE_RIGHT_TOP                  real4_2d_exchange_right_top
#define TK_3D_EXCHANGE_RIGHT_TOP                  real4_3d_exchange_right_top
#define TK_4D_EXCHANGE_RIGHT_TOP                  real4_4d_exchange_right_top
#define TK_BDY_EXCHANGE_LEFT_RIGHT                real4_bdy_exchange_left_right
#define TK_BDY_EXCHANGE_BOTTOM_TOP                real4_bdy_exchange_bottom_top

! mpp_tk_subgrid.inc
#define REORDER_ADD_SUBGRID_2D_TK                 reorder_add_subgrid_2d_real4
#define REORDER_SUBGRID_2D_TK                     reorder_subgrid_2d_real4
#define REORDER_ADD_SUBGRID_2D3D_TK               reorder_add_subgrid_2d3d_real4
#define REORDER_SUBGRID_2D3D_TK                   reorder_subgrid_2d3d_real4
#define REORDER_ADD_SUBGRID_3D_TK                 reorder_add_subgrid_3d_real4
#define REORDER_SUBGRID_3D_TK                     reorder_subgrid_3d_real4
#define REORDER_SUBGRID_4D_TK                     reorder_subgrid_4d_real4
#define INPUT_REORDER_TK                          input_reorder_real4

! mpp_tk_linear.inc
#define CL_SETUP_TK                               cl_setup_real4
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID_4D cartesian_to_linear_real4_subgrid_subgrid_4d
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID_4D linear_to_cartesian_real4_subgrid_subgrid_4d
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID    cartesian_to_linear_real4_subgrid_subgrid
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID    linear_to_cartesian_real4_subgrid_subgrid
#define CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID       cartesian_to_linear_real4_grid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID       global_to_linear_real4_subgrid_subgrid
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID       linear_to_global_real4_subgrid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID_4D    global_to_linear_real4_subgrid_subgrid_4d
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID_4D    linear_to_global_real4_subgrid_subgrid_4d
#define GLOBAL_TO_LINEAR_TK_GRID_SUBGRID          global_to_linear_real4_grid_subgrid
#define MYPACK_TK_GRID                            mypack_real4_grid
#define MYPACK_TK_SUBGRID                         mypack_real4_subgrid
#define MYPACK_TK_SUBGRID_4D                      mypack_real4_subgrid_4d
#define MYPACK_TK_SUBGRID_SLICE                   mypack_real4_subgrid_slice
#define MYUNPACK_TK_GRID                          myunpack_real4_grid
#define MYUNPACK_TK_SUBGRID                       myunpack_real4_subgrid
#define MYUNPACK_TK_SUBGRID_4D                    myunpack_real4_subgrid_4d
#define MYUNPACK_TK_SUBGRID_SLICE                 myunpack_real4_subgrid_slice
#define MYPACK_GLOBAL_TK_GRID                     mypack_global_real4_grid
#define MYPACK_GLOBAL_TK_SUBGRID                  mypack_global_real4_subgrid
#define MYPACK_GLOBAL_TK_SUBGRID_4D               mypack_global_real4_subgrid_4d
#define MYPACK_GLOBAL_TK_SUBGRID_SLICE            mypack_global_real4_subgrid_slice
#define MYUNPACK_GLOBAL_TK_GRID                   myunpack_global_real4_grid
#define MYUNPACK_GLOBAL_TK_SUBGRID                myunpack_global_real4_subgrid
#define MYUNPACK_GLOBAL_TK_SUBGRID_4D             myunpack_global_real4_subgrid_4d
#define MYUNPACK_GLOBAL_TK_SUBGRID_SLICE          myunpack_global_real4_subgrid_slice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
