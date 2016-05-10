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

  integer(ik4) :: nocnp
  integer(ik4) :: iocnbeg , iocnend

  logical :: llake = .false.
  logical :: ldcsst = .false.
  logical :: lseaice = .false.

  real(rkx) , pointer , dimension(:) :: cprate
  real(rkx) , pointer , dimension(:) :: czenith
  real(rkx) , pointer , dimension(:) :: deltas
  real(rkx) , pointer , dimension(:) :: dhlake
  real(rkx) , pointer , dimension(:) :: drag
  real(rkx) , pointer , dimension(:) :: tskin
  real(rkx) , pointer , dimension(:) :: dwrlwf
  real(rkx) , pointer , dimension(:) :: emiss
  real(rkx) , pointer , dimension(:) :: evpr
  real(rkx) , pointer , dimension(:) :: ht      ! hgt
  real(rkx) , pointer , dimension(:) :: hpbl    ! hpbl
  real(rkx) , pointer , dimension(:) :: lat     ! xlat
  real(rkx) , pointer , dimension(:) :: ncprate
  real(rkx) , pointer , dimension(:) :: prcp
  real(rkx) , pointer , dimension(:) :: q2m
  real(rkx) , pointer , dimension(:) :: qv      ! qvatm
  real(rkx) , pointer , dimension(:) :: rhox    ! rhox
  real(rkx) , pointer , dimension(:) :: rlwf    ! rlwf
  real(rkx) , pointer , dimension(:) :: rswf    ! rswf
  real(rkx) , pointer , dimension(:) :: scvk
  real(rkx) , pointer , dimension(:) :: sent
  real(rkx) , pointer , dimension(:) :: sfice
  real(rkx) , pointer , dimension(:) :: sts     ! tatm
  real(rkx) , pointer , dimension(:) :: sfps    ! sfps
  real(rkx) , pointer , dimension(:) :: snag
  real(rkx) , pointer , dimension(:) :: sncv
  real(rkx) , pointer , dimension(:) :: sm
  real(rkx) , pointer , dimension(:) :: sst
  real(rkx) , pointer , dimension(:) :: t2m
  real(rkx) , pointer , dimension(:) :: tatm
  real(rkx) , pointer , dimension(:) :: taux
  real(rkx) , pointer , dimension(:) :: tauy
  real(rkx) , pointer , dimension(:) :: tdeltas
  real(rkx) , pointer , dimension(:) :: tgb     ! tground2
  real(rkx) , pointer , dimension(:) :: tgbrd
  real(rkx) , pointer , dimension(:) :: tgrd
  real(rkx) , pointer , dimension(:) :: u10m
  real(rkx) , pointer , dimension(:) :: v10m
  real(rkx) , pointer , dimension(:) :: um10
  real(rkx) , pointer , dimension(:) :: usw     ! uatm
  real(rkx) , pointer , dimension(:) :: vsw     ! vatm
  real(rkx) , pointer , dimension(:) :: ustr    ! ustar
  real(rkx) , pointer , dimension(:) :: zoo     ! zo
  real(rkx) , pointer , dimension(:) :: rhoa    ! xdens

  real(rkx) , pointer , dimension(:) :: laketa
  real(rkx) , pointer , dimension(:) :: lakhi
  real(rkx) , pointer , dimension(:) :: lakaveice
  real(rkx) , pointer , dimension(:,:) :: laktlake
  integer(ik4) , pointer , dimension(:) :: ilake

  real(rkx) , pointer , dimension(:) :: swdiral
  real(rkx) , pointer , dimension(:) :: lwdiral
  real(rkx) , pointer , dimension(:) :: swdifal
  real(rkx) , pointer , dimension(:) :: lwdifal

  integer(ik4) , pointer , dimension(:) :: mask
  integer(ik4) , pointer , dimension(:) :: xmask
  integer(ik4) , pointer , dimension(:) :: icpl

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
    call getmem1d(q2m,1,nocnp,'ocn_internal:q2m')
    call getmem1d(qv,1,nocnp,'ocn_internal:qv')
    call getmem1d(rhox,1,nocnp,'ocn_internal:rhox')
    call getmem1d(rlwf,1,nocnp,'ocn_internal:rlwf')
    call getmem1d(rswf,1,nocnp,'ocn_internal:rswf')
    call getmem1d(sent,1,nocnp,'ocn_internal:sent')
    call getmem1d(sfps,1,nocnp,'ocn_internal:sfps')
    call getmem1d(sts,1,nocnp,'ocn_internal:sts')
    call getmem1d(t2m,1,nocnp,'ocn_internal:t2m')
    call getmem1d(tatm,1,nocnp,'ocn_internal:tatm')
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
      call getmem1d(xmask,1,nocnp,'ocn_internal:xmask')
      call getmem1d(ilake,1,nocnp,'ocn_internal:ilake')
      call getmem1d(scvk,1,nocnp,'ocn_internal:scvk')
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
