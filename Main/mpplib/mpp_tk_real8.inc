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

! Instantiate the mpp_tk_*.inc templates for real(rk8).
! It must be included BEFORE the type/kind methods!

#include "mpp_tk_clear.inc"

#define TK_REAL

#define TK_TYPE                                   real(rk8)
#define TK_KIND                                   rk8
#define TK_MPI                                    mpi_real8
#define TK_NAME                                   'real8'
#define TK_VEC1                                   r8vector1
#define TK_VEC2                                   r8vector2
#define TK_SUBGRID                                r8subgrid
#define TK_GLOBAL_SUBGRID                         global_r8subgrid
#define TK_MISSVAL                                dmissval

! mpp_tk_collectives.inc
#define BCAST_TK                                  bcast_real8
#define BCAST_ARR_TK                              bcast_arr_real8
#define BCAST_MATR_TK                             bcast_matr_real8
#define SUMALL_TK                                 sumall_real8
#define MAXALL_TK                                 maxall_real8
#define MEANALL_1D_TK                             meanall_1D_real8
#define MEANALL_TK                                meanall_real8
#define MINALL_TK                                 minall_real8
#define EXCHANGE_ARRAY_TK                         exchange_array_r8
#define CYCLIC_EXCHANGE_ARRAY_TK                  cyclic_exchange_array_r8
#define TK_COLUMN_REDUCE                          real8_column_reduce
#define TK_ROW_REDUCE                             real8_row_reduce

! mpp_tk_distribute.inc
#define TK_2D_DO_DISTRIBUTE                       real8_2d_do_distribute
#define TK_2D_DISTRIBUTE                          real8_2d_distribute
#define TK_3D_DISTRIBUTE                          real8_3d_distribute
#define TK_4D_DISTRIBUTE                          real8_4d_distribute
#define TK_2D_DO_SUB_DISTRIBUTE                   real8_2d_do_sub_distribute
#define TK_2D_SUB_DISTRIBUTE                      real8_2d_sub_distribute
#define TK_3D_SUB_DISTRIBUTE                      real8_3d_sub_distribute
#define TK_2D_DO_COLLECT                          real8_2d_do_collect
#define TK_COLLECT_PLANES                         real8_collect_planes
#define TK_2D_COLLECT                             real8_2d_collect
#define TK_2D_3D_COLLECT                          real8_2d_3d_collect
#define TK_3D_COLLECT                             real8_3d_collect
#define TK_3D_2D_COLLECT                          real8_3d_2d_collect
#define TK_4D_COLLECT                             real8_4d_collect
#define TK_4D_2D_COLLECT                          real8_4d_2d_collect
#define TK_2D_DO_SUB_COLLECT                      real8_2d_do_sub_collect
#define TK_2D_SUB_COLLECT                         real8_2d_sub_collect
#define TK_3D_SUB_COLLECT                         real8_3d_sub_collect
#define TK_2D_GRID_FILL_EXTEND1                   real8_2d_grid_fill_extend1
#define TK_2D_GRID_FILL_EXTEND2                   real8_2d_grid_fill_extend2

! mpp_tk_exchange.inc
#define TK_2D_EXCHANGE                            real8_2d_exchange
#define TK_3D_EXCHANGE                            real8_3d_exchange
#define TK_4D_EXCHANGE                            real8_4d_exchange
#define TK_2D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real8_2d_exchange_left_right_bottom_top
#define TK_3D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real8_3d_exchange_left_right_bottom_top
#define TK_4D_EXCHANGE_LEFT_RIGHT_BOTTOM_TOP      real8_4d_exchange_left_right_bottom_top
#define TK_2D_EXCHANGE_LEFT_RIGHT                 real8_2d_exchange_left_right
#define TK_3D_EXCHANGE_LEFT_RIGHT                 real8_3d_exchange_left_right
#define TK_4D_EXCHANGE_LEFT_RIGHT                 real8_4d_exchange_left_right
#define TK_2D_EXCHANGE_BOTTOM_TOP                 real8_2d_exchange_bottom_top
#define TK_3D_EXCHANGE_BOTTOM_TOP                 real8_3d_exchange_bottom_top
#define TK_4D_EXCHANGE_BOTTOM_TOP                 real8_4d_exchange_bottom_top
#define TK_2D_EXCHANGE_LEFT_BOTTOM                real8_2d_exchange_left_bottom
#define TK_3D_EXCHANGE_LEFT_BOTTOM                real8_3d_exchange_left_bottom
#define TK_4D_EXCHANGE_LEFT_BOTTOM                real8_4d_exchange_left_bottom
#define TK_2D_EXCHANGE_RIGHT_TOP                  real8_2d_exchange_right_top
#define TK_3D_EXCHANGE_RIGHT_TOP                  real8_3d_exchange_right_top
#define TK_4D_EXCHANGE_RIGHT_TOP                  real8_4d_exchange_right_top
#define TK_BDY_EXCHANGE_LEFT_RIGHT                real8_bdy_exchange_left_right
#define TK_BDY_EXCHANGE_BOTTOM_TOP                real8_bdy_exchange_bottom_top

! mpp_tk_subgrid.inc
#define REORDER_ADD_SUBGRID_2D_TK                 reorder_add_subgrid_2d_real8
#define REORDER_SUBGRID_2D_TK                     reorder_subgrid_2d_real8
#define REORDER_ADD_SUBGRID_2D3D_TK               reorder_add_subgrid_2d3d_real8
#define REORDER_SUBGRID_2D3D_TK                   reorder_subgrid_2d3d_real8
#define REORDER_ADD_SUBGRID_3D_TK                 reorder_add_subgrid_3d_real8
#define REORDER_SUBGRID_3D_TK                     reorder_subgrid_3d_real8
#define REORDER_SUBGRID_4D_TK                     reorder_subgrid_4d_real8
#define INPUT_REORDER_TK                          input_reorder_real8

! mpp_tk_linear.inc
#define CL_SETUP_TK                               cl_setup_real8
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID_4D cartesian_to_linear_real8_subgrid_subgrid_4d
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID_4D linear_to_cartesian_real8_subgrid_subgrid_4d
#define CARTESIAN_TO_LINEAR_TK_SUBGRID_SUBGRID    cartesian_to_linear_real8_subgrid_subgrid
#define LINEAR_TO_CARTESIAN_TK_SUBGRID_SUBGRID    linear_to_cartesian_real8_subgrid_subgrid
#define CARTESIAN_TO_LINEAR_TK_GRID_SUBGRID       cartesian_to_linear_real8_grid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID       global_to_linear_real8_subgrid_subgrid
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID       linear_to_global_real8_subgrid_subgrid
#define GLOBAL_TO_LINEAR_TK_SUBGRID_SUBGRID_4D    global_to_linear_real8_subgrid_subgrid_4d
#define LINEAR_TO_GLOBAL_TK_SUBGRID_SUBGRID_4D    linear_to_global_real8_subgrid_subgrid_4d
#define GLOBAL_TO_LINEAR_TK_GRID_SUBGRID          global_to_linear_real8_grid_subgrid
#define MYPACK_TK_GRID                            mypack_real8_grid
#define MYPACK_TK_SUBGRID                         mypack_real8_subgrid
#define MYPACK_TK_SUBGRID_4D                      mypack_real8_subgrid_4d
#define MYPACK_TK_SUBGRID_SLICE                   mypack_real8_subgrid_slice
#define MYUNPACK_TK_GRID                          myunpack_real8_grid
#define MYUNPACK_TK_SUBGRID                       myunpack_real8_subgrid
#define MYUNPACK_TK_SUBGRID_4D                    myunpack_real8_subgrid_4d
#define MYUNPACK_TK_SUBGRID_SLICE                 myunpack_real8_subgrid_slice
#define MYPACK_GLOBAL_TK_GRID                     mypack_global_real8_grid
#define MYPACK_GLOBAL_TK_SUBGRID                  mypack_global_real8_subgrid
#define MYPACK_GLOBAL_TK_SUBGRID_4D               mypack_global_real8_subgrid_4d
#define MYPACK_GLOBAL_TK_SUBGRID_SLICE            mypack_global_real8_subgrid_slice
#define MYUNPACK_GLOBAL_TK_GRID                   myunpack_global_real8_grid
#define MYUNPACK_GLOBAL_TK_SUBGRID                myunpack_global_real8_subgrid
#define MYUNPACK_GLOBAL_TK_SUBGRID_4D             myunpack_global_real8_subgrid_4d
#define MYUNPACK_GLOBAL_TK_SUBGRID_SLICE          myunpack_global_real8_subgrid_slice

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
