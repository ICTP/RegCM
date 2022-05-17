!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_cmip6

  use mod_intkinds
  use mod_realkinds
  use mod_cmip6_helper
  use mod_message
  use mod_date
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_kdinterp
  use mod_write
  use mod_vectutil
  use mod_mksst
  use mod_humid
  use mod_hgt
  use mod_vertint
  use netcdf

  implicit none

  private

  character(len=*) , parameter , public :: mpihr_version = 'v20190710'
  character(len=*) , parameter , public :: mpihr_version1 = 'v20190815'
  character(len=*) , parameter , public :: hadmm_version = 'v20191207'
  character(len=*) , parameter , public :: hadmm_version1 = 'v20210331'
  character(len=*) , parameter , public :: hadmm_version2 = 'v20201019'
  character(len=*) , parameter , public :: hadmm_version3 = 'v20201113'
  character(len=*) , parameter , public :: normm_version = 'v20191108'
  character(len=*) , parameter , public :: normm_version1 = 'v20210319'
  character(len=*) , parameter , public :: normm_version2 = 'v20200702'
  character(len=*) , parameter , public :: normm_version3 = 'v20200218'
  character(len=*) , parameter , public :: normm_version4 = 'v20210319'
  character(len=*) , parameter , public :: cnrm_version = 'v20181205'
  character(len=*) , parameter , public :: cnrm_version1 = 'v20181206'
  character(len=*) , parameter , public :: cesm_version = 'v20190514'
  character(len=*) , parameter , public :: cesm_version1 = 'v20200528'
  character(len=*) , parameter , public :: cesm_version2 = 'v20200210'
  character(len=*) , parameter , public :: ecea_version = 'v20210324'
  character(len=*) , parameter , public :: ecea_version1 = 'v20200310'

  public :: init_cmip6 , get_cmip6 , conclude_cmip6

  ! Pressure levels to interpolate to if dataset is on model sigma levels.
  integer(ik4) , parameter :: nipl = 38
  real(rkx) , target , dimension(nipl) :: fplev = &
    [ 1000.0_rkx, 975.0_rkx, 950.0_rkx, 925.0_rkx, 900.0_rkx, 875.0_rkx, &
       850.0_rkx, 825.0_rkx, 800.0_rkx, 775.0_rkx, 750.0_rkx, 700.0_rkx, &
       650.0_rkx, 600.0_rkx, 550.0_rkx, 500.0_rkx, 450.0_rkx, 425.0_rkx, &
       400.0_rkx, 350.0_rkx, 300.0_rkx, 250.0_rkx, 225.0_rkx, 200.0_rkx, &
       175.0_rkx, 150.0_rkx, 125.0_rkx, 100.0_rkx,  70.0_rkx,  50.0_rkx, &
        30.0_rkx,  20.0_rkx,  10.0_rkx,   7.0_rkx,   5.0_rkx,   3.0_rkx, &
        2.0_rkx,    1.0_rkx ]

!    [  1.0_rkx,   2.0_rkx,   3.0_rkx,   5.0_rkx,   7.0_rkx,  10.0_rkx, &
!      20.0_rkx,  30.0_rkx,  50.0_rkx,  70.0_rkx, 100.0_rkx, 125.0_rkx, &
!     150.0_rkx, 175.0_rkx, 200.0_rkx, 225.0_rkx, 250.0_rkx, 300.0_rkx, &
!     350.0_rkx, 400.0_rkx, 425.0_rkx, 450.0_rkx, 500.0_rkx, 550.0_rkx, &
!     600.0_rkx, 650.0_rkx, 700.0_rkx, 750.0_rkx, 775.0_rkx, 800.0_rkx, &
!     825.0_rkx, 850.0_rkx, 875.0_rkx, 900.0_rkx, 925.0_rkx, 950.0_rkx, &
!    975.0_rkx, 1000.0_rkx ]

  type(cmip6_2d_var) , pointer :: orog => null( )
  type(cmip6_2d_var) , pointer :: ps => null( )
  type(cmip6_3d_var) , pointer :: ta => null( )
  type(cmip6_3d_var) , pointer :: qa => null( )
  type(cmip6_3d_var) , pointer :: ua => null( )
  type(cmip6_3d_var) , pointer :: va => null( )
  type(cmip6_3d_var) , pointer :: zg => null( )

  real(rkx) , dimension(:) , pointer :: sigmar
  real(rkx) :: pss , pst

  real(rkx) , dimension(:,:,:) , pointer :: pa_in , zp_in
  real(rkx) , dimension(:,:,:) , pointer :: tvar
  real(rkx) , dimension(:,:,:) , pointer :: uvar
  real(rkx) , dimension(:,:,:) , pointer :: vvar
  real(rkx) , dimension(:,:,:) , pointer :: qvar
  real(rkx) , dimension(:,:,:) , pointer :: zvar

  real(rkx) , dimension(:,:,:) , pointer :: tah
  real(rkx) , dimension(:,:,:) , pointer :: uah
  real(rkx) , dimension(:,:,:) , pointer :: vah
  real(rkx) , dimension(:,:,:) , pointer :: qah
  real(rkx) , dimension(:,:,:) , pointer :: zgh

  real(rkx) , pointer , dimension(:,:,:) :: dv , du
  real(rkx) , pointer , dimension(:,:,:) :: hv , hu
  real(rkx) , pointer , dimension(:,:) :: topou , topov

  integer(ik4) :: nkin

  logical , parameter :: only_coord = .true.

  contains

    subroutine init_cmip6(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      integer(ik4) :: k
      select case (cmip6_model)
        case ( 'MPI-ESM1-2-HR' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_mpihr(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          call read_2d_mpihr(idate,ps,only_coord)
          call read_3d_mpihr(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_mpihr(idate,ua,only_coord)
          call read_3d_mpihr(idate,va,only_coord)
          call read_3d_mpihr(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:mpihr:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:mpihr:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:mpihr:zp_in')
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:mpihr:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:mpihr:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:mpihr:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:mpihr:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:mpihr:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:mpihr:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:mpihr:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:mpihr:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:mpihr:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:mpihr:zgh')
        case ( 'HadGEM3-GC31-MM' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_hadmm(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          call read_2d_hadmm(idate,ps,only_coord)
          call read_3d_hadmm(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_hadmm(idate,ua,only_coord)
          call read_3d_hadmm(idate,va,only_coord)
          call read_3d_hadmm(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:hadmm:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:hadmm:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:hadmm:zp_in')
          do k = 1 , ta%nk
            zp_in(:,:,k) = ta%vcoord%ak(k) + ta%vcoord%bk(k) * orog%var(:,:)
          end do
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:hadmm:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:hadmm:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:hadmm:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:hadmm:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:hadmm:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:hadmm:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:hadmm:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:hadmm:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:hadmm:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:hadmm:zgh')
        case ( 'CNRM-ESM2-1' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_cnrm(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          call read_2d_cnrm(idate,ps,only_coord)
          call read_3d_cnrm(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_cnrm(idate,ua,only_coord)
          call read_3d_cnrm(idate,va,only_coord)
          call read_3d_cnrm(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:cnrm:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cnrm:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cnrm:zp_in')
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cnrm:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:cnrm:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:cnrm:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:cnrm:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cnrm:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:cnrm:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:cnrm:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:cnrm:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:cnrm:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:cnrm:zgh')
        case ( 'EC-Earth3' )
          allocate(ps,ua,va,ta,qa,zg,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          zg%vname = 'zg'
          orog%vname = 'orog'
          call read_fx_ecea(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          zg%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          zg%hcoord => orog%hcoord
          call read_2d_ecea(idate,ps,only_coord)
          stop
          call read_3d_cnrm(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_cnrm(idate,ua,only_coord)
          call read_3d_cnrm(idate,va,only_coord)
          call read_3d_cnrm(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:cnrm:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cnrm:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cnrm:zp_in')
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cnrm:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:cnrm:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:cnrm:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:cnrm:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cnrm:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:cnrm:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:cnrm:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:cnrm:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:cnrm:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:cnrm:zgh')
        case ( 'CESM2' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_cesm(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          call read_2d_cesm(idate,ps,only_coord)
          call read_3d_cesm(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_cesm(idate,ua,only_coord)
          call read_3d_cesm(idate,va,only_coord)
          call read_3d_cesm(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:cesm:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cesm:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:cesm:zp_in')
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cesm:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:cesm:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:cesm:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:cesm:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:cesm:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:cesm:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:cesm:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:cesm:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:cesm:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:cesm:zgh')
        case ( 'NorESM2-MM' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_normm(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          call read_2d_normm(idate,ps,only_coord)
          call read_3d_normm(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_normm(idate,ua,only_coord)
          call read_3d_normm(idate,va,only_coord)
          call read_3d_normm(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:normm:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:normm:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:normm:zp_in')
          do k = 1 , ta%nk
            zp_in(:,:,k) = ta%vcoord%ak(k) + ta%vcoord%bk(k) * orog%var(:,:)
          end do
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:normm:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:normm:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:normm:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:normm:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:normm:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:normm:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:normm:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:normm:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:normm:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:normm:zgh')
        case default
          call die(__FILE__,'__LINE__ : Unsupported cmip6 model.',-1)
      end select

      if ( idynamic == 3 ) then
        call getmem2d(topou,1,jx,1,iy,'cmip6:topou')
        call getmem2d(topov,1,jx,1,iy,'cmip6:topov')
        call getmem3d(du,1,jx,1,iy,1,nkin,'cmip6:du')
        call getmem3d(dv,1,jx,1,iy,1,nkin,'cmip6:dv')
        call getmem3d(hu,1,jx,1,iy,1,nkin,'cmip6:hu')
        call getmem3d(hv,1,jx,1,iy,1,nkin,'cmip6:hv')
        call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
        call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
        call ucrs2dot(topou,topogm,jx,iy,i_band)
        call vcrs2dot(topov,topogm,jx,iy,i_crm)
      end if

      write (stdout,*) 'Read in Static fields OK'
    end subroutine init_cmip6

    subroutine get_cmip6(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      integer(ik4) :: i , j , k

      select case (cmip6_model)
        case ( 'MPI-ESM1-2-HR' )
!$OMP SECTIONS
!$OMP SECTION
          call read_2d_mpihr(idate,ps)
!$OMP SECTION
          call read_3d_mpihr(idate,ua)
!$OMP SECTION
          call read_3d_mpihr(idate,va)
!$OMP SECTION
          call read_3d_mpihr(idate,ta)
!$OMP SECTION
          call read_3d_mpihr(idate,qa)
!$OMP END SECTIONS
          call sph2mxr(qa%var,qa%ni,qa%nj,qa%nk)
          do k = 1 , ta%nk
            do j = 1 , ta%nj
              do i = 1 , ta%ni
                pa_in(i,j,k) = ta%vcoord%ak(k) + ta%vcoord%bk(k) * ps%var(i,j)
              end do
            end do
          end do
          pa_in = pa_in * 0.01_rkx
!$OMP SECTIONS
!$OMP SECTION
          call intlin(uvar,ua%var,pa_in,ua%ni,ua%nj,ua%nk,fplev,nkin)
!$OMP SECTION
          call intlin(vvar,va%var,pa_in,va%ni,va%nj,va%nk,fplev,nkin)
!$OMP SECTION
          call intlog(tvar,ta%var,pa_in,ta%ni,ta%nj,ta%nk,fplev,nkin)
!$OMP SECTION
          call intlin(qvar,qa%var,pa_in,qa%ni,qa%nj,qa%nk,fplev,nkin)
!$OMP END SECTIONS
          ps%var = ps%var * 0.01_rkx
          call htsig(zp_in,ta%var,pa_in,qa%var,ps%var,orog%var)
          call height(zvar,zp_in,ta%var,ps%var,pa_in,orog%var, &
                      ta%ni,ta%nj,ta%nk,fplev,nkin)
        case ( 'HadGEM3-GC31-MM' )
!$OMP SECTIONS
!$OMP SECTION
          call read_2d_hadmm(idate,ps)
!$OMP SECTION
          call read_3d_hadmm(idate,ua)
!$OMP SECTION
          call read_3d_hadmm(idate,va)
!$OMP SECTION
          call read_3d_hadmm(idate,ta)
!$OMP SECTION
          call read_3d_hadmm(idate,qa)
          call sph2mxr(qa%var,qa%ni,qa%nj,qa%nk)
!$OMP END SECTIONS
          call psig(ta%var,zp_in,pa_in,ps%var,orog%var,ta%ni,ta%nj,ta%nk)
          pa_in = pa_in * 0.01_rkx
!$OMP SECTIONS
!$OMP SECTION
          call intlin(uvar,ua%var,pa_in,ua%ni,ua%nj,ua%nk,fplev,nkin)
!$OMP SECTION
          call intlin(vvar,va%var,pa_in,va%ni,va%nj,va%nk,fplev,nkin)
!$OMP SECTION
          call intlog(tvar,ta%var,pa_in,ta%ni,ta%nj,ta%nk,fplev,nkin)
!$OMP SECTION
          call intlin(qvar,qa%var,pa_in,qa%ni,qa%nj,qa%nk,fplev,nkin)
!$OMP END SECTIONS
          ps%var = ps%var * 0.01_rkx
          call htsig(zp_in,ta%var,pa_in,qa%var,ps%var,orog%var)
          call height(zvar,zp_in,ta%var,ps%var,pa_in,orog%var, &
                      ta%ni,ta%nj,ta%nk,fplev,nkin)
        case ( 'CNRM-ESM2-1' )
!$OMP SECTIONS
!$OMP SECTION
          call read_2d_cnrm(idate,ps)
!$OMP SECTION
          call read_3d_cnrm(idate,ua)
!$OMP SECTION
          call read_3d_cnrm(idate,va)
!$OMP SECTION
          call read_3d_cnrm(idate,ta)
!$OMP SECTION
          call read_3d_cnrm(idate,qa)
!$OMP END SECTIONS
          call sph2mxr(qa%var,qa%ni,qa%nj,qa%nk)
          do k = 1 , ta%nk
            do j = 1 , ta%nj
              do i = 1 , ta%ni
                pa_in(i,j,k) = ta%vcoord%ak(k) + ta%vcoord%bk(k) * ps%var(i,j)
              end do
            end do
          end do
          pa_in = pa_in * 0.01_rkx
!$OMP SECTIONS
!$OMP SECTION
          call intlin(uvar,ua%var,pa_in,ua%ni,ua%nj,ua%nk,fplev,nkin)
!$OMP SECTION
          call intlin(vvar,va%var,pa_in,va%ni,va%nj,va%nk,fplev,nkin)
!$OMP SECTION
          call intlog(tvar,ta%var,pa_in,ta%ni,ta%nj,ta%nk,fplev,nkin)
!$OMP SECTION
          call intlin(qvar,qa%var,pa_in,qa%ni,qa%nj,qa%nk,fplev,nkin)
!$OMP END SECTIONS
          ps%var = ps%var * 0.01_rkx
          call htsig(zp_in,ta%var,pa_in,qa%var,ps%var,orog%var)
          call height(zvar,zp_in,ta%var,ps%var,pa_in,orog%var, &
                      ta%ni,ta%nj,ta%nk,fplev,nkin)
        case ( 'CESM2' , 'NorESM2-MM' )
!$OMP SECTIONS
!$OMP SECTION
          call read_2d_cesm(idate,ps)
!$OMP SECTION
          call read_3d_cesm(idate,ua)
!$OMP SECTION
          call read_3d_cesm(idate,va)
!$OMP SECTION
          call read_3d_cesm(idate,ta)
!$OMP SECTION
          call read_3d_cesm(idate,qa)
!$OMP END SECTIONS
          call sph2mxr(qa%var,qa%ni,qa%nj,qa%nk)
          do k = 1 , ta%nk
            do j = 1 , ta%nj
              do i = 1 , ta%ni
                pa_in(i,j,k) = ta%vcoord%ak(k) * ta%vcoord%p0 + &
                       ta%vcoord%bk(k) * ps%var(i,j)
              end do
            end do
          end do
          pa_in = pa_in * 0.01_rkx
!$OMP SECTIONS
!$OMP SECTION
          call intlin(uvar,ua%var,pa_in,ua%ni,ua%nj,ua%nk,fplev,nkin)
!$OMP SECTION
          call intlin(vvar,va%var,pa_in,va%ni,va%nj,va%nk,fplev,nkin)
!$OMP SECTION
          call intlog(tvar,ta%var,pa_in,ta%ni,ta%nj,ta%nk,fplev,nkin)
!$OMP SECTION
          call intlin(qvar,qa%var,pa_in,qa%ni,qa%nj,qa%nk,fplev,nkin)
!$OMP END SECTIONS
          ps%var = ps%var * 0.01_rkx
          call htsig(zp_in,ta%var,pa_in,qa%var,ps%var,orog%var)
          call height(zvar,zp_in,ta%var,ps%var,pa_in,orog%var, &
                      ta%ni,ta%nj,ta%nk,fplev,nkin)
        case default
          call die(__FILE__,'__LINE__ : Unsupported cmip6 model.',-1)
      end select

      write (stdout,*) 'Read in fields at Date: ', tochar(idate)

!$OMP SECTIONS
!$OMP SECTION
      call h_interpolate_cont(ta%hint(1),tvar,tah)
!$OMP SECTION
      call h_interpolate_cont(qa%hint(1),qvar,qah)
!$OMP SECTION
      call h_interpolate_cont(ta%hint(1),zvar,zgh)
!$OMP END SECTIONS
      if ( idynamic == 3 ) then
!$OMP SECTIONS
        call h_interpolate_cont(ua%hint(2),uvar,uah)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),vvar,dv)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(3),uvar,du)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(3),vvar,vah)
        call pju%wind_rotate(uah,dv)
        call pjv%wind_rotate(du,vah)
!$OMP END SECTIONS
      else
!$OMP SECTIONS
        call h_interpolate_cont(ua%hint(2),uvar,uah)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),vvar,vah)
        call pjd%wind_rotate(uah,vah)
!$OMP END SECTIONS
      end if

      if ( idynamic == 3 ) then
        call ucrs2dot(hu,zgh,jx,iy,nkin,i_band)
        call vcrs2dot(hv,zgh,jx,iy,nkin,i_crm)
        call intzps(ps4,topogm,tah,zgh,pss,sigmar,pst, &
                    xlat,yeardayfrac(idate),jx,iy,nkin)
        call intz3(ts4,tah,zgh,topogm,jx,iy,nkin,0.6_rkx,0.5_rkx,0.85_rkx)
      else
        call intgtb(pa,za,tlayer,topogm,tah,zgh,pss,sigmar,pst,jx,iy,nkin)
        call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
        call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
        call intv3(ts4,tah,ps4,pss,sigmar,ptop,pst,jx,iy,nkin)
      end if

      call readsst(ts4,idate)

      if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
        call intz1(u4,uah,zud4,hu,topou,jx,iy,kz,nkin,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
        call intz1(v4,vah,zvd4,hv,topov,jx,iy,kz,nkin,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
        call intz1(t4,tah,z0,zgh,topogm,jx,iy,kz,nkin,0.6_rkx,0.5_rkx,0.85_rkx)
!$OMP SECTION
        call intz1(q4,qah,z0,zgh,topogm,jx,iy,kz,nkin,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
      else
!$OMP SECTIONS
!$OMP SECTION
        call intv1(u4,uah,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP SECTION
        call intv1(v4,vah,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP SECTION
        call intv2(t4,tah,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin)
!$OMP SECTION
        call intv1(q4,qah,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP END SECTIONS
      end if
    end subroutine get_cmip6

    subroutine conclude_cmip6( )
      implicit none
      if ( idynamic == 3 ) then
        call h_interpolator_destroy(orog%hint(1))
        call h_interpolator_destroy(orog%hint(2))
        call h_interpolator_destroy(orog%hint(3))
      else
        call h_interpolator_destroy(orog%hint(1))
        call h_interpolator_destroy(orog%hint(2))
      end if
      deallocate(orog%hint)
    end subroutine conclude_cmip6

    subroutine read_hcoord_mpihr(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_mpi:lon')
      call getmem1d(lat,1,nlat,'cmip6_mpi:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_mpihr

    subroutine read_vcoord_mpihr(ncid,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(ap,1,nlev,'cmip6_mpi:ap')
      call getmem1d(b,1,nlev,'cmip6_mpi:b')
      istatus = nf90_inq_varid(ncid,'ap',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_vcoord_mpihr

    recursive subroutine read_3d_mpihr(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = year
        if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 1
        end if
        if ( y < 2015 .and. ( v%vname == 'ua' .or. v%vname == 'va' ) ) then
          ver = mpihr_version1
        else
          ver = mpihr_version
        end if
        write(v%filename,'(a,i4,a,i4,a)') &
          trim(cmip6_path(y,'6hrLev',ver,v%vname)), &
          y, '01010600-', y+1, '01010000.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_mpihr(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_mpihr(v%ncid,v%vcoord%ak,v%vcoord%bk)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:mpihr:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_mpihr(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_mpihr

    recursive subroutine read_2d_mpihr(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = (year / 5) * 5
        if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 5
        end if
        write(v%filename,'(a,i4,a,i4,a)') &
          trim(cmip6_path(y,'6hrLev',mpihr_version,v%vname)), &
          y, '01010600-', y+5, '01010000.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_mpihr(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:mpihr:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_mpihr(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_mpihr

    recursive subroutine read_fx_mpihr(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(mpihr_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_mpihr(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:mpihr:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_mpihr

    subroutine read_hcoord_hadmm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_mpi:lon')
      call getmem1d(lat,1,nlat,'cmip6_mpi:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_hadmm

    subroutine read_vcoord_hadmm(ncid,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(ap,1,nlev,'cmip6_mpi:lev')
      call getmem1d(b,1,nlev,'cmip6_mpi:b')
      istatus = nf90_inq_varid(ncid,'lev',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lev var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_vcoord_hadmm

    recursive subroutine read_3d_hadmm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y , yp , m , mp
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        y = year
        m = month / 3 * 3 + 1
        if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
          y = y - 1
          m = 10
        end if
        yp = y
        mp = m + 3
        if ( mp > 12 ) then
          yp = yp + 1
          mp = 1
        end if
        if ( y > 2015 ) then
          ver = hadmm_version3
        else
          ver = hadmm_version2
        end if
        write(v%filename,'(a,i4,i0.2,a,i4,i0.2,a)') &
          trim(cmip6_path(y,'6hrLev',ver,v%vname)), &
          y, m, '010600-', yp, mp, '010000.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_hadmm(v%ncid,v%vcoord%ak,v%vcoord%bk)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:hadmm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_hadmm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_hadmm

    recursive subroutine read_2d_hadmm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year >= 2015 .and. year < 2020 ) then
          if ( year == 2015 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2010,'6hrLev',hadmm_version2,v%vname)), &
              '201001010600-', '201501010000.nc'
          else
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2015,'6hrLev',hadmm_version3,v%vname)), &
              '201501010600-', '202001010000.nc'
          end if
        else if ( year == 2020 .and. month == 1 .and. &
                 day == 1 .and. hour == 0 ) then
          write(v%filename,'(a,a,a)') &
            trim(cmip6_path(2015,'6hrLev',hadmm_version3,v%vname)), &
            '201501010600-', '202001010000.nc'
        else if ( year >= 2010 .and. year < 2015 ) then
          if ( year == 2010 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2000,'6hrLev',hadmm_version2,v%vname)), &
              '200001010600-', '201001010000.nc'
          else
            write(v%filename,'(a,a,a)') &
              trim(cmip6_path(2010,'6hrLev',hadmm_version2,v%vname)), &
              '201001010600-', '201501010000.nc'
          end if
        else
          y = (year / 10) * 10
          if ( y == year .and. month == 1 .and. day == 1 .and. hour == 0 ) then
            y = y - 10
          end if
          if ( y > 2015 ) then
            write(v%filename,'(a,i4,a,i4,a)') &
              trim(cmip6_path(y,'6hrLev',hadmm_version3,v%vname)), &
              y, '01010600-', y+10, '01010000.nc'
          else
            write(v%filename,'(a,i4,a,i4,a)') &
              trim(cmip6_path(y,'6hrLev',hadmm_version2,v%vname)), &
              y, '01010600-', y+10, '01010000.nc'
          end if
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:hadmm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_hadmm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_hadmm

    recursive subroutine read_fx_hadmm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(hadmm_version1,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:hadmm:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_hadmm

    subroutine read_hcoord_normm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_nor:lon')
      call getmem1d(lat,1,nlat,'cmip6_nor:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_normm

    subroutine read_vcoord_normm(ncid,a,b,p0)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: a , b
      real(rkx) , intent(out) :: p0
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(a,1,nlev,'cmip6_nor:a')
      call getmem1d(b,1,nlev,'cmip6_nor:b')
      istatus = nf90_inq_varid(ncid,'a',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,a)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(ncid,'p0',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(ncid,ivarid,p0)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read p0 var')
    end subroutine read_vcoord_normm

    recursive subroutine read_3d_normm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( idate < 2010010106 ) then
          y = (year / 10) * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version3,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        else if ( idate < 2015010106 ) then
          y = 2010
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version3,v%vname)) // &
            '201001010600-201501010000.nc'
        else if ( idate < 2022010106 ) then
          y = 2015
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)) // &
            '201501010600-202101010000.nc'
        else
          y = (year-2021)/10*10 + 2021
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_normm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_normm(v%ncid,v%vcoord%ak,v%vcoord%bk,v%vcoord%p0)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_normm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_normm

    recursive subroutine read_2d_normm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( idate < 2010010106 ) then
          y = (year / 10) * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version1,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        else if ( idate < 2015010106 ) then
          y = 2010
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version1,v%vname)) // &
            '201001010600-201501010000.nc'
        else if ( idate < 2021010106 ) then
          y = 2015
          v%filename = trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)) // &
            '201501010600-202101010000.nc'
        else
          y = (year-2021)/10*10 + 2021
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(y,'6hrLev',normm_version2,v%vname)), &
            y, '01010600-', y+10, '01010000.nc'
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_normm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_normm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_normm

    recursive subroutine read_fx_normm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(normm_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:normm:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_normm

    subroutine read_hcoord_cnrm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_cnrm:lon')
      call getmem1d(lat,1,nlat,'cmip6_cnrm:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_cnrm

    subroutine read_vcoord_cnrm(ncid,ap,b)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: ap , b
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(ap,1,nlev,'cmip6_cnrm:ap')
      call getmem1d(b,1,nlev,'cmip6_cnrm:b')
      istatus = nf90_inq_varid(ncid,'ap',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,ap)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
    end subroutine read_vcoord_cnrm

    recursive subroutine read_3d_cnrm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y1 , y2 , m1 , m2
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( v%vname == 'ua' .or. v%vname == 'va' ) then
          y1 = year
          m1 = ((month-1)/4) * 4 + 1
          if ( m1 == month .and. day == 1 .and. hour == 0 ) then
            m1 = m1 - 4
            if ( m1 < 1 ) then
              m1 = 12 + m1
              y1 = y1 - 1
            end if
          end if
          m2 = m1 + 4
          if ( m2 > 12 ) then
            m2 = 1
            y2 = y1 + 1
          else
            y2 = y1
          end if
        else
          y1 = year
          m1 = ((month-1)/6) * 6 + 1
          if ( m1 == month .and. day == 1 .and. hour == 0 ) then
            m1 = m1 - 6
            if ( m1 < 1 ) then
              m1 = 12 + m1
              y1 = y1 - 1
            end if
          end if
          m2 = m1 + 6
          if ( m2 > 12 ) then
            m2 = 1
            y2 = y1 + 1
          else
            y2 = y1
          end if
        end if
        ver = cnrm_version1
        write(v%filename,'(a,i4,i0.2,a,i4,i0.2,a)') &
          trim(cmip6_path(year,'6hrLev',ver,v%vname)), &
          y1, m1, '010600-', y2, m2, '010000.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_cnrm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_cnrm(v%ncid,v%vcoord%ak,v%vcoord%bk)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:cnrm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_cnrm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_cnrm

    recursive subroutine read_2d_cnrm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year < 2000 ) then
          write(v%filename,'(a,a)') &
            trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
            '195001010600-200001010000.nc'
        else if ( year < 2015 ) then
          if ( year == 2000 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '195001010600-200001010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '200001010600-201501010000.nc'
          end if
        else if ( year < 2065 ) then
          if ( year == 2015 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '200001010600-201501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '201501010600-206501010000.nc'
          end if
        else
          if ( year == 2065 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '201501010600-206501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cnrm_version1,v%vname)), &
              '206501010600-210101010000.nc'
          end if
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_cnrm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cnrm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_cnrm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_cnrm

    recursive subroutine read_fx_cnrm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(cnrm_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cnrm:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_cnrm

    subroutine read_hcoord_cesm(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_cesm:lon')
      call getmem1d(lat,1,nlat,'cmip6_cesm:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_cesm

    subroutine read_vcoord_cesm(ncid,a,b,p0)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: a , b
      real(rkx) , intent(out) :: p0
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(a,1,nlev,'cmip6_cesm:a')
      call getmem1d(b,1,nlev,'cmip6_cesm:b')
      istatus = nf90_inq_varid(ncid,'a',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,a)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(ncid,'p0',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(ncid,ivarid,p0)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read p0 var')
    end subroutine read_vcoord_cesm

    recursive subroutine read_3d_cesm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour , y1 , y2
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif
      character(len=16) :: ver

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year < 2010 ) then
          if ( v%vname /= 'ta' ) then
            ver = cesm_version
          else
            ver = cesm_version2
          end if
          y1 = year/10*10
          y2 = y1+9
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',ver,v%vname)), &
            y1, '01010000-', y2, '12311800.nc'
        else if ( year >= 2015 .and. year < 2095 ) then
          y1 = (year-2015)/10*10 + 2015
          y2 = y1+9
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
            y1, '01010000-', y2, '12311800.nc'
        else
          if ( idate < 2015010100 ) then
            if ( v%vname /= 'ta' ) then
              ver = cesm_version
            else
              ver = cesm_version2
            end if
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',ver,v%vname)), &
              '201001010000-201501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
              '209501010000-210101010000.nc'
          end if
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_cesm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_cesm(v%ncid,v%vcoord%ak,v%vcoord%bk,v%vcoord%p0)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:cesm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_cesm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_cesm

    recursive subroutine read_2d_cesm(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: y1
      integer(ik4) :: year , month , day , hour
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        if ( year < 2010 ) then
          y1 = year/10 * 10
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version,v%vname)), &
            y1,'01010000-',y1+9,'12311800.nc'
        else if ( year >= 2015 .and. year < 2095 ) then
          y1 = (year-2015)/10*10 + 2015
          write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
            y1,'01010000-',y1+9,'12311800.nc'
        else
          if ( idate < 2015010100 ) then
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version,v%vname)), &
              '201001010000-201501010000.nc'
          else
            write(v%filename,'(a,a)') &
              trim(cmip6_path(year,'6hrLev',cesm_version1,v%vname)), &
              '209501010000-210101010000.nc'
          end if
        end if
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_cesm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cesm:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_cesm(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_cesm

    recursive subroutine read_fx_cesm(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(cesm_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:cesm:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_cesm

    subroutine read_hcoord_ecea(ncid,lon,lat)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: lon , lat
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlon , nlat
      istatus = nf90_inq_dimid(ncid,'lon',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lon dim')
      istatus = nf90_inq_dimid(ncid,'lat',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lat dim')
      call getmem1d(lon,1,nlon,'cmip6_ecea:lon')
      call getmem1d(lat,1,nlat,'cmip6_ecea:lat')
      istatus = nf90_inq_varid(ncid,'lon',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lon var')
      istatus = nf90_get_var(ncid,ivarid,lon)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lon var')
      istatus = nf90_inq_varid(ncid,'lat',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lat var')
      istatus = nf90_get_var(ncid,ivarid,lat)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read lat var')
    end subroutine read_hcoord_ecea

    subroutine read_vcoord_ecea(ncid,a,b,p0)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , pointer , dimension(:) , intent(inout) :: a , b
      real(rkx) , intent(out) :: p0
      integer(ik4) :: istatus , idimid , ivarid
      integer(ik4) :: nlev
      istatus = nf90_inq_dimid(ncid,'lev',idimid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find lev dim')
      istatus = nf90_inquire_dimension(ncid,idimid,len=nlev)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire lev dim')
      call getmem1d(a,1,nlev,'cmip6_ecea:a')
      call getmem1d(b,1,nlev,'cmip6_ecea:b')
      istatus = nf90_inq_varid(ncid,'a',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(ncid,ivarid,a)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(ncid,'b',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(ncid,ivarid,b)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(ncid,'p0',ivarid)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(ncid,ivarid,p0)
      call cmip6_error(istatus,__FILE__,__LINE__,'Error read p0 var')
    end subroutine read_vcoord_ecea

    recursive subroutine read_3d_ecea(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_3d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(4) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',ecea_version1,v%vname)), &
            year, '01010000-', year, '12311800.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_ecea(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%vcoord) ) then
        allocate(v%vcoord)
        call read_vcoord_ecea(v%ncid,v%vcoord%ak,v%vcoord%bk,v%vcoord%p0)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        v%nk = size(v%vcoord%ak)
        call getmem3d(v%var,1,v%ni,1,v%nj,1,v%nk,'cmip6:ecea:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ', &
          v%ni,'x',v%nj,'x',v%nk
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_3d_ecea(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = size(v%var,3)
      icount(4) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_3d_ecea

    recursive subroutine read_2d_ecea(idate,v,lonlyc)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      logical , optional , intent(in) :: lonlyc
      integer(ik4) :: istatus , idimid , it , irec
      integer(ik4) :: year , month , day , hour
      character(len=32) :: timecal , timeunit
      integer(ik4) , dimension(3) :: istart , icount
      real(rk8) , dimension(2) :: times
      type(rcm_time_interval) :: tdif

      if ( v%ncid == -1 ) then
        call split_idate(idate, year, month, day, hour)
        write(v%filename,'(a,i4,a,i4,a)') &
            trim(cmip6_path(year,'6hrLev',ecea_version1,v%vname)), &
            year,'01010000-',year,'12311800.nc'
#ifdef DEBUG
        write(stderr,*) 'Opening ',trim(v%filename)
#endif
        istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      end if
      if ( .not. associated(v%hcoord) ) then
        allocate(v%hcoord)
        call read_hcoord_ecea(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      end if
      if ( .not. associated(v%var) ) then
        v%ni = size(v%hcoord%lon1d)
        v%nj = size(v%hcoord%lat1d)
        call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:ecea:'//trim(v%vname))
#ifdef DEBUG
        write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
        if ( present(lonlyc) ) then
          if ( lonlyc ) return
        end if
      end if
      if ( v%ivar == -1 ) then
        istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searchong '//trim(v%vname)//' var in file '// &
          trim(v%filename)//'.')
      end if
      if ( v%nrec == -1 ) then
        istatus = nf90_inq_dimid(v%ncid,'time',idimid)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error find time dim')
        istatus = nf90_inquire_dimension(v%ncid,idimid,len=v%nrec)
        call cmip6_error(istatus,__FILE__,__LINE__,'Error inquire time dim')
        istatus = nf90_inq_varid(v%ncid,'time',it)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error searching time var in file '//trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"calendar",timecal)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute calendar in file '//&
          trim(v%filename)//'.')
        istatus = nf90_get_att(v%ncid,it,"units",timeunit)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time attribute units in file '//trim(v%filename)//'.')
        istart(1) = 1
        icount(1) = 2
        istatus = nf90_get_var(v%ncid,it,times,istart(1:1),icount(1:1))
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error reading time from file '//trim(v%filename)//'.')
        v%first_date = timeval2date(times(1),timeunit,timecal)
      end if

      tdif = idate - v%first_date
      irec = nint(tohours(tdif)/6.0) + 1

      if ( irec > v%nrec ) then
        istatus = nf90_close(v%ncid)
        call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
        v%ncid = -1
        v%ivar = -1
        v%nrec = -1
        irec = 1
#ifdef DEBUG
        write(stderr, *) 'Closed file, switching to the next in series...'
#endif
        call read_2d_ecea(idate,v)
      end if
      istart(1) = 1
      istart(2) = 1
      istart(3) = irec
      icount(1) = size(v%var,1)
      icount(2) = size(v%var,2)
      icount(3) = 1
      istatus = nf90_get_var(v%ncid,v%ivar,v%var,istart,icount)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
    end subroutine read_2d_ecea

    recursive subroutine read_fx_ecea(v)
      implicit none
      type(cmip6_2d_var) , pointer , intent(inout) :: v
      integer(ik4) :: istatus

      v%filename = trim(cmip6_fxpath(ecea_version,v%vname))
#ifdef DEBUG
      write(stderr,*) 'Opening ',trim(v%filename)
#endif
      istatus = nf90_open(v%filename,nf90_nowrite,v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error opening file '//trim(v%filename)//'.')
      allocate(v%hcoord)
      call read_hcoord_hadmm(v%ncid,v%hcoord%lon1d,v%hcoord%lat1d)
      v%ni = size(v%hcoord%lon1d)
      v%nj = size(v%hcoord%lat1d)
      call getmem2d(v%var,1,v%ni,1,v%nj,'cmip6:ecea:'//trim(v%vname))
#ifdef DEBUG
      write(stderr,*) 'Input shape for ',trim(v%vname),' = ',v%ni,'x',v%nj
#endif
      istatus = nf90_inq_varid(v%ncid,v%vname,v%ivar)
      call cmip6_error(istatus,__FILE__,__LINE__, &
        'Error searchong '//trim(v%vname)//' var in file '// &
        trim(v%filename)//'.')
      istatus = nf90_get_var(v%ncid,v%ivar,v%var)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error read variable '//v%vname//' from '//trim(v%filename)//'.')
      istatus = nf90_close(v%ncid)
      call cmip6_error(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(v%filename)//'.')
    end subroutine read_fx_ecea

end module mod_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
