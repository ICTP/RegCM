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
  use mod_runparams, only : idcsst, lakemod, iseaice, dtsrf, iocnflx
  use mod_regcm_types

  implicit none (type, external)

  public

  real(rkx) :: dtocn, dtlake, dtsst

  real(rkx), parameter :: aarea = 0.02_rkx
  real(rkx), parameter :: age3 = 0.3_rkx
  real(rkx), parameter :: threedays = 86400.0_rkx*3.0_rkx  ! 3 days

  integer(ik4) :: nocnp
  integer(ik4) :: iocnbeg, iocnend

  logical :: llake = .false.
  logical :: ldcsst = .false.
  logical :: lseaice = .false.

  real(rkx), pointer, contiguous, dimension(:) :: cprate => null( )
  real(rkx), pointer, contiguous, dimension(:) :: czenith => null( )
  real(rkx), pointer, contiguous, dimension(:) :: deltas => null( )
  real(rkx), pointer, contiguous, dimension(:) :: dhlake => null( )
  real(rkx), pointer, contiguous, dimension(:) :: drag => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tskin => null( )
  real(rkx), pointer, contiguous, dimension(:) :: dwrlwf => null( )
  real(rkx), pointer, contiguous, dimension(:) :: emiss => null( )
  real(rkx), pointer, contiguous, dimension(:) :: evpr => null( )
  real(rkx), pointer, contiguous, dimension(:) :: ht => null( )      ! hgt
  real(rkx), pointer, contiguous, dimension(:) :: hpbl => null( )    ! hpbl
  real(rkx), pointer, contiguous, dimension(:) :: lat => null( )     ! xlat
  real(rkx), pointer, contiguous, dimension(:) :: ncprate => null( )
  real(rkx), pointer, contiguous, dimension(:) :: prcp => null( )
  real(rkx), pointer, contiguous, dimension(:) :: q2m => null( )
  real(rkx), pointer, contiguous, dimension(:) :: qv => null( )      ! qvatm
  real(rkx), pointer, contiguous, dimension(:) :: rhox => null( )    ! rhox
  real(rkx), pointer, contiguous, dimension(:) :: rlwf => null( )    ! rlwf
  real(rkx), pointer, contiguous, dimension(:) :: rswf => null( )    ! rswf
  real(rkx), pointer, contiguous, dimension(:) :: sent => null( )
  real(rkx), pointer, contiguous, dimension(:) :: sfice => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tatm => null( )    ! tatm
  real(rkx), pointer, contiguous, dimension(:) :: sfps => null( )    ! sfps
  real(rkx), pointer, contiguous, dimension(:) :: snag => null( )
  real(rkx), pointer, contiguous, dimension(:) :: sncv => null( )
  real(rkx), pointer, contiguous, dimension(:) :: sm => null( )
  real(rkx), pointer, contiguous, dimension(:) :: sst => null( )
  real(rkx), pointer, contiguous, dimension(:) :: t2m => null( )
  real(rkx), pointer, contiguous, dimension(:) :: sfta => null( )
  real(rkx), pointer, contiguous, dimension(:) :: patm => null( )
  real(rkx), pointer, contiguous, dimension(:) :: taux => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tauy => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tdeltas => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tgb => null( )     ! tground2
  real(rkx), pointer, contiguous, dimension(:) :: tgbrd => null( )
  real(rkx), pointer, contiguous, dimension(:) :: tgrd => null( )
  real(rkx), pointer, contiguous, dimension(:) :: u10m => null( )
  real(rkx), pointer, contiguous, dimension(:) :: v10m => null( )
  real(rkx), pointer, contiguous, dimension(:) :: um10 => null( )
  real(rkx), pointer, contiguous, dimension(:) :: usw => null( )     ! uatm
  real(rkx), pointer, contiguous, dimension(:) :: vsw => null( )     ! vatm
  real(rkx), pointer, contiguous, dimension(:) :: ustr => null( )    ! ustar
  real(rkx), pointer, contiguous, dimension(:) :: zoo => null( )     ! zo
  real(rkx), pointer, contiguous, dimension(:) :: rhoa => null( )    ! xdens
  real(rkx), pointer, contiguous, dimension(:) :: ram1 => null( )
  real(rkx), pointer, contiguous, dimension(:) :: rah1 => null( )
  real(rkx), pointer, contiguous, dimension(:) :: br => null( )

  real(rkx), pointer, contiguous, dimension(:) :: laketa => null( )
  real(rkx), pointer, contiguous, dimension(:) :: lakhi => null( )
  real(rkx), pointer, contiguous, dimension(:) :: lakaveice => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: laktlake => null( )
  integer(ik4), pointer, contiguous, dimension(:) :: ilake => null( )

  real(rkx), pointer, contiguous, dimension(:) :: swdiral => null( )
  real(rkx), pointer, contiguous, dimension(:) :: lwdiral => null( )
  real(rkx), pointer, contiguous, dimension(:) :: swdifal => null( )
  real(rkx), pointer, contiguous, dimension(:) :: lwdifal => null( )

  integer(ik4), pointer, contiguous, dimension(:) :: mask => null( )
  integer(ik4), pointer, contiguous, dimension(:) :: icpl => null( )
  integer(ik4), pointer, contiguous, dimension(:) :: omask => null( )

  contains

  subroutine allocate_mod_ocn_internal(co)
    implicit none (type, external)
    type (masked_comm), intent(in) :: co
    nocnp = co%linear_npoint_sg(myid+1)
    iocnbeg = 1
    iocnend = nocnp
    call getmem(dhlake,1,nocnp,'ocn_internal:dhlake')
    call getmem(drag,1,nocnp,'ocn_internal:drag')
    call getmem(dwrlwf,1,nocnp,'ocn_internal:dwrlwf')
    call getmem(emiss,1,nocnp,'ocn_internal:emiss')
    call getmem(evpr,1,nocnp,'ocn_internal:evpr')
    call getmem(ht,1,nocnp,'ocn_internal:ht')
    call getmem(hpbl,1,nocnp,'ocn_internal:hpbl')
    call getmem(lat,1,nocnp,'ocn_internal:lat')
    call getmem(omask,1,nocnp,'ocn_internal:omask')
    call getmem(q2m,1,nocnp,'ocn_internal:q2m')
    call getmem(qv,1,nocnp,'ocn_internal:qv')
    call getmem(rhox,1,nocnp,'ocn_internal:rhox')
    call getmem(rlwf,1,nocnp,'ocn_internal:rlwf')
    call getmem(rswf,1,nocnp,'ocn_internal:rswf')
    call getmem(sent,1,nocnp,'ocn_internal:sent')
    call getmem(sfps,1,nocnp,'ocn_internal:sfps')
    call getmem(tatm,1,nocnp,'ocn_internal:tatm')
    call getmem(t2m,1,nocnp,'ocn_internal:t2m')
    call getmem(sfta,1,nocnp,'ocn_internal:sfta')
    call getmem(patm,1,nocnp,'ocn_internal:patm')
    call getmem(taux,1,nocnp,'ocn_internal:taux')
    call getmem(tauy,1,nocnp,'ocn_internal:tauy')
    call getmem(tgb,1,nocnp,'ocn_internal:tgb')
    call getmem(tgrd,1,nocnp,'ocn_internal:tgrd')
    call getmem(tgbrd,1,nocnp,'ocn_internal:tgbrd')
    call getmem(u10m,1,nocnp,'ocn_internal:u10m')
    call getmem(v10m,1,nocnp,'ocn_internal:v10m')
    call getmem(usw,1,nocnp,'ocn_internal:usw')
    call getmem(vsw,1,nocnp,'ocn_internal:vsw')
    call getmem(ustr,1,nocnp,'ocn_internal:ustr')
    call getmem(zoo,1,nocnp,'ocn_internal:zoo')
    call getmem(rhoa,1,nocnp,'ocn_internal:rhoa')
    call getmem(mask,1,nocnp,'ocn_internal:mask')
    call getmem(icpl,1,nocnp,'ocn_internal:icpl')
    call getmem(czenith,1,nocnp,'ocn_internal:czenith')
    call getmem(cprate,1,nocnp,'ocn_internal:cprate')
    call getmem(ncprate,1,nocnp,'ocn_internal:ncprate')
    call getmem(prcp,1,nocnp,'ocn_internal:prcp')
    call getmem(um10,1,nocnp,'ocn_internal:um10')
    call getmem(ram1,1,nocnp,'ocn_internal:ram1')
    call getmem(rah1,1,nocnp,'ocn_internal:rah1')
    call getmem(br,1,nocnp,'ocn_internal:br')
    if ( lakemod == 1 ) llake = .true.
    if ( idcsst == 1 ) ldcsst = .true.
    if ( iseaice == 1 ) lseaice = .true.
    icpl(:) = 0
    dtocn = dtsrf
    if ( ldcsst ) then
      call getmem(deltas,1,nocnp,'ocn_internal:deltas')
      call getmem(tskin,1,nocnp,'ocn_internal:tskin')
      call getmem(tdeltas,1,nocnp,'ocn_internal:tdeltas')
      call getmem(sst,1,nocnp,'ocn_internal:sst')
      dtsst = dtsrf
    end if
    if ( lseaice .or. llake ) then
      call getmem(ilake,1,nocnp,'ocn_internal:ilake')
      call getmem(sfice,1,nocnp,'ocn_internal:sfice')
      call getmem(snag,1,nocnp,'ocn_internal:snag')
      call getmem(sncv,1,nocnp,'ocn_internal:sncv')
      call getmem(sm,1,nocnp,'ocn_internal:sm')
      if ( llake ) then
        call getmem(laketa,1,nocnp,'ocn_internal:laketa')
        call getmem(lakhi,1,nocnp,'ocn_internal:lakhi')
        call getmem(lakaveice,1,nocnp,'ocn_internal:lakaveice')
        call getmem(laktlake,1,nocnp,1,ndpmax,'ocn_internal:laktlake')
        dtlake = dtsrf
      end if
    end if
    call getmem(swdiral,1,nocnp,'ocn_internal:swdiral')
    call getmem(lwdiral,1,nocnp,'ocn_internal:lwdiral')
    call getmem(swdifal,1,nocnp,'ocn_internal:swdifal')
    call getmem(lwdifal,1,nocnp,'ocn_internal:lwdifal')
  end subroutine allocate_mod_ocn_internal
!
end module mod_ocn_internal
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
