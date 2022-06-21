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
  use mod_cmip6_ecea
  use mod_cmip6_cesm
  use mod_cmip6_cnrm
  use mod_cmip6_hadmm
  use mod_cmip6_miroc6
  use mod_cmip6_normm
  use mod_cmip6_mpihr

  implicit none

  private

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
        case ( 'MIROC6' )
          allocate(ps,ua,va,ta,qa,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          orog%vname = 'orog'
          call read_fx_miroc6(orog)
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
          call read_2d_miroc6(idate,ps,only_coord)
          call read_3d_miroc6(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_miroc6(idate,ua,only_coord)
          call read_3d_miroc6(idate,va,only_coord)
          call read_3d_miroc6(idate,qa,only_coord)
          nkin = nipl
          call getmem1d(sigmar,1,nkin,'cmip6:miroc6:sigmar')
          do k = 1 , nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem3d(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:miroc6:pa_in')
          call getmem3d(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip6:miroc6:zp_in')
          call getmem3d(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:miroc6:tvar')
          call getmem3d(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip6:miroc6:uvar')
          call getmem3d(vvar,1,va%ni,1,va%nj,1,nkin,'cmip6:miroc6:vvar')
          call getmem3d(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip6:miroc6:qvar')
          call getmem3d(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip6:miroc6:zvar')
          call getmem3d(tah,1,jx,1,iy,1,nkin,'cmip6:miroc6:tah')
          call getmem3d(qah,1,jx,1,iy,1,nkin,'cmip6:miroc6:qah')
          call getmem3d(uah,1,jx,1,iy,1,nkin,'cmip6:miroc6:uah')
          call getmem3d(vah,1,jx,1,iy,1,nkin,'cmip6:miroc6:vah')
          call getmem3d(zgh,1,jx,1,iy,1,nkin,'cmip6:miroc6:zgh')
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
          call read_2d_cnrm(idate,ps,only_coord)
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
          call die(__FILE__,'Unsupported cmip6 model. Stop at line ',__LINE__)
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
        case ( 'MIROC6' )
!$OMP SECTIONS
!$OMP SECTION
          call read_2d_miroc6(idate,ps)
!$OMP SECTION
          call read_3d_miroc6(idate,ua)
!$OMP SECTION
          call read_3d_miroc6(idate,va)
!$OMP SECTION
          call read_3d_miroc6(idate,ta)
!$OMP SECTION
          call read_3d_miroc6(idate,qa)
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
          call die(__FILE__,'Unsupported cmip6 model. Stop at line ',__LINE__)
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

end module mod_cmip6

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
