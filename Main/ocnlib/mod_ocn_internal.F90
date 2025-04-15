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

module mod_ocn_internal
!
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_runparams , only : idcsst , lakemod , iseaice , dtsrf , iocnflx
  use mod_regcm_types

  implicit none

  public

  real(rkx) :: dtocn , dtlake , dtsst

  real(rkx) , parameter :: aarea = 0.02_rkx
  real(rkx) , parameter :: age3 = 0.3_rkx
  real(rkx) , parameter :: threedays = 86400.0_rkx*3.0_rkx  ! 3 days

  integer(ik4) :: nocnp
  integer(ik4) :: iocnbeg , iocnend

  logical :: llake = .false.
  logical :: ldcsst = .false.
  logical :: lseaice = .false.

  real(rkx) , pointer , dimension(:) :: cprate => null( )
  real(rkx) , pointer , dimension(:) :: czenith => null( )
  real(rkx) , pointer , dimension(:) :: deltas => null( )
  real(rkx) , pointer , dimension(:) :: dhlake => null( )
  real(rkx) , pointer , dimension(:) :: drag => null( )
  real(rkx) , pointer , dimension(:) :: tskin => null( )
  real(rkx) , pointer , dimension(:) :: dwrlwf => null( )
  real(rkx) , pointer , dimension(:) :: emiss => null( )
  real(rkx) , pointer , dimension(:) :: evpr => null( )
  real(rkx) , pointer , dimension(:) :: ht => null( )      ! hgt
  real(rkx) , pointer , dimension(:) :: hpbl => null( )    ! hpbl
  real(rkx) , pointer , dimension(:) :: lat => null( )     ! xlat
  real(rkx) , pointer , dimension(:) :: ncprate => null( )
  real(rkx) , pointer , dimension(:) :: prcp => null( )
  real(rkx) , pointer , dimension(:) :: q2m => null( )
  real(rkx) , pointer , dimension(:) :: qv => null( )      ! qvatm
  real(rkx) , pointer , dimension(:) :: rhox => null( )    ! rhox
  real(rkx) , pointer , dimension(:) :: rlwf => null( )    ! rlwf
  real(rkx) , pointer , dimension(:) :: rswf => null( )    ! rswf
  real(rkx) , pointer , dimension(:) :: sent => null( )
  real(rkx) , pointer , dimension(:) :: sfice => null( )
  real(rkx) , pointer , dimension(:) :: tatm => null( )    ! tatm
  real(rkx) , pointer , dimension(:) :: sfps => null( )    ! sfps
  real(rkx) , pointer , dimension(:) :: snag => null( )
  real(rkx) , pointer , dimension(:) :: sncv => null( )
  real(rkx) , pointer , dimension(:) :: sm => null( )
  real(rkx) , pointer , dimension(:) :: sst => null( )
  real(rkx) , pointer , dimension(:) :: t2m => null( )
  real(rkx) , pointer , dimension(:) :: sfta => null( )
  real(rkx) , pointer , dimension(:) :: patm => null( )
  real(rkx) , pointer , dimension(:) :: taux => null( )
  real(rkx) , pointer , dimension(:) :: tauy => null( )
  real(rkx) , pointer , dimension(:) :: tdeltas => null( )
  real(rkx) , pointer , dimension(:) :: tgb => null( )     ! tground2
  real(rkx) , pointer , dimension(:) :: tgbrd => null( )
  real(rkx) , pointer , dimension(:) :: tgrd => null( )
  real(rkx) , pointer , dimension(:) :: u10m => null( )
  real(rkx) , pointer , dimension(:) :: v10m => null( )
  real(rkx) , pointer , dimension(:) :: um10 => null( )
  real(rkx) , pointer , dimension(:) :: usw => null( )     ! uatm
  real(rkx) , pointer , dimension(:) :: vsw => null( )     ! vatm
  real(rkx) , pointer , dimension(:) :: ustr => null( )    ! ustar
  real(rkx) , pointer , dimension(:) :: zoo => null( )     ! zo
  real(rkx) , pointer , dimension(:) :: rhoa => null( )    ! xdens
  real(rkx) , pointer , dimension(:) :: ram1 => null( )
  real(rkx) , pointer , dimension(:) :: rah1 => null( )
  real(rkx) , pointer , dimension(:) :: br => null( )

  real(rkx) , pointer , dimension(:) :: laketa => null( )
  real(rkx) , pointer , dimension(:) :: lakhi => null( )
  real(rkx) , pointer , dimension(:) :: lakaveice => null( )
  real(rkx) , pointer , dimension(:,:) :: laktlake => null( )
  integer(ik4) , pointer , dimension(:) :: ilake => null( )

  real(rkx) , pointer , dimension(:) :: swdiral => null( )
  real(rkx) , pointer , dimension(:) :: lwdiral => null( )
  real(rkx) , pointer , dimension(:) :: swdifal => null( )
  real(rkx) , pointer , dimension(:) :: lwdifal => null( )

  integer(ik4) , pointer , dimension(:) :: mask => null( )
  integer(ik4) , pointer , dimension(:) :: icpl => null( )
  integer(ik4) , pointer , dimension(:) :: omask => null( )

  contains

  subroutine allocate_mod_ocn_internal(co)
    implicit none
    type (masked_comm) , intent(in) :: co
    nocnp = co%linear_npoint_sg(myid+1)
    iocnbeg = 1
    iocnend = nocnp
    call getmem1d(dhlake,1,nocnp,'ocn_internal:dhlake')
    call getmem1d(drag,1,nocnp,'ocn_internal:drag')
    call getmem1d(dwrlwf,1,nocnp,'ocn_internal:dwrlwf')
    call getmem1d(emiss,1,nocnp,'ocn_internal:emiss')
    call getmem1d(evpr,1,nocnp,'ocn_internal:evpr')
    call getmem1d(ht,1,nocnp,'ocn_internal:ht')
    call getmem1d(hpbl,1,nocnp,'ocn_internal:hpbl')
    call getmem1d(lat,1,nocnp,'ocn_internal:lat')
    call getmem1d(omask,1,nocnp,'ocn_internal:omask')
    call getmem1d(q2m,1,nocnp,'ocn_internal:q2m')
    call getmem1d(qv,1,nocnp,'ocn_internal:qv')
    call getmem1d(rhox,1,nocnp,'ocn_internal:rhox')
    call getmem1d(rlwf,1,nocnp,'ocn_internal:rlwf')
    call getmem1d(rswf,1,nocnp,'ocn_internal:rswf')
    call getmem1d(sent,1,nocnp,'ocn_internal:sent')
    call getmem1d(sfps,1,nocnp,'ocn_internal:sfps')
    call getmem1d(tatm,1,nocnp,'ocn_internal:tatm')
    call getmem1d(t2m,1,nocnp,'ocn_internal:t2m')
    call getmem1d(sfta,1,nocnp,'ocn_internal:sfta')
    call getmem1d(patm,1,nocnp,'ocn_internal:patm')
    call getmem1d(taux,1,nocnp,'ocn_internal:taux')
    call getmem1d(tauy,1,nocnp,'ocn_internal:tauy')
    call getmem1d(tgb,1,nocnp,'ocn_internal:tgb')
    call getmem1d(tgrd,1,nocnp,'ocn_internal:tgrd')
    call getmem1d(tgbrd,1,nocnp,'ocn_internal:tgbrd')
    call getmem1d(u10m,1,nocnp,'ocn_internal:u10m')
    call getmem1d(v10m,1,nocnp,'ocn_internal:v10m')
    call getmem1d(usw,1,nocnp,'ocn_internal:usw')
    call getmem1d(vsw,1,nocnp,'ocn_internal:vsw')
    call getmem1d(ustr,1,nocnp,'ocn_internal:ustr')
    call getmem1d(zoo,1,nocnp,'ocn_internal:zoo')
    call getmem1d(rhoa,1,nocnp,'ocn_internal:rhoa')
    call getmem1d(mask,1,nocnp,'ocn_internal:mask')
    call getmem1d(icpl,1,nocnp,'ocn_internal:icpl')
    call getmem1d(czenith,1,nocnp,'ocn_internal:czenith')
    call getmem1d(cprate,1,nocnp,'ocn_internal:cprate')
    call getmem1d(ncprate,1,nocnp,'ocn_internal:ncprate')
    call getmem1d(prcp,1,nocnp,'ocn_internal:prcp')
    call getmem1d(um10,1,nocnp,'ocn_internal:um10')
    call getmem1d(ram1,1,nocnp,'ocn_internal:ram1')
    call getmem1d(rah1,1,nocnp,'ocn_internal:rah1')
    call getmem1d(br,1,nocnp,'ocn_internal:br')
    if ( lakemod == 1 ) llake = .true.
    if ( idcsst == 1 ) ldcsst = .true.
    if ( iseaice == 1 ) lseaice = .true.
    icpl(:) = 0
    dtocn = dtsrf
    if ( ldcsst ) then
      call getmem1d(deltas,1,nocnp,'ocn_internal:deltas')
      call getmem1d(tskin,1,nocnp,'ocn_internal:tskin')
      call getmem1d(tdeltas,1,nocnp,'ocn_internal:tdeltas')
      call getmem1d(sst,1,nocnp,'ocn_internal:sst')
      dtsst = dtsrf
    end if
    if ( lseaice .or. llake ) then
      call getmem1d(ilake,1,nocnp,'ocn_internal:ilake')
      call getmem1d(sfice,1,nocnp,'ocn_internal:sfice')
      call getmem1d(snag,1,nocnp,'ocn_internal:snag')
      call getmem1d(sncv,1,nocnp,'ocn_internal:sncv')
      call getmem1d(sm,1,nocnp,'ocn_internal:sm')
      if ( llake ) then
        call getmem1d(laketa,1,nocnp,'ocn_internal:laketa')
        call getmem1d(lakhi,1,nocnp,'ocn_internal:lakhi')
        call getmem1d(lakaveice,1,nocnp,'ocn_internal:lakaveice')
        call getmem2d(laktlake,1,nocnp,1,ndpmax,'ocn_internal:laktlake')
        dtlake = dtsrf
      end if
    end if
    call getmem1d(swdiral,1,nocnp,'ocn_internal:swdiral')
    call getmem1d(lwdiral,1,nocnp,'ocn_internal:lwdiral')
    call getmem1d(swdifal,1,nocnp,'ocn_internal:swdifal')
    call getmem1d(lwdifal,1,nocnp,'ocn_internal:lwdifal')
  end subroutine allocate_mod_ocn_internal
!
end module mod_ocn_internal
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
