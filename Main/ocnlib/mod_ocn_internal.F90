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

  real(rk8) :: dtocn , dtlake , dtsst

  real(rk8) , parameter :: aarea = 0.02D0
  real(rk8) , parameter :: age3 = 0.3D0

  integer(ik4) :: nocnp
  integer(ik4) :: iocnbeg , iocnend

  logical :: llake = .false.
  logical :: ldcsst = .false.
  logical :: lseaice = .false.

  real(rk8) , pointer , dimension(:) :: cprate
  real(rk8) , pointer , dimension(:) :: czenith
  real(rk8) , pointer , dimension(:) :: deltas
  real(rk8) , pointer , dimension(:) :: dhlake
  real(rk8) , pointer , dimension(:) :: drag
  real(rk8) , pointer , dimension(:) :: tskin
  real(rk8) , pointer , dimension(:) :: dwrlwf
  real(rk8) , pointer , dimension(:) :: emiss
  real(rk8) , pointer , dimension(:) :: evpr
  real(rk8) , pointer , dimension(:) :: ht      ! hgt
  real(rk8) , pointer , dimension(:) :: hpbl    ! hpbl
  real(rk8) , pointer , dimension(:) :: lat     ! xlat
  real(rk8) , pointer , dimension(:) :: ncprate
  real(rk8) , pointer , dimension(:) :: prcp
  real(rk8) , pointer , dimension(:) :: q2m
  real(rk8) , pointer , dimension(:) :: qv      ! qvatm
  real(rk8) , pointer , dimension(:) :: rhox    ! rhox
  real(rk8) , pointer , dimension(:) :: rlwf    ! rlwf
  real(rk8) , pointer , dimension(:) :: rswf    ! rswf
  real(rk8) , pointer , dimension(:) :: scvk
  real(rk8) , pointer , dimension(:) :: sent
  real(rk8) , pointer , dimension(:) :: sfice
  real(rk8) , pointer , dimension(:) :: sts     ! tatm
  real(rk8) , pointer , dimension(:) :: sfps    ! sfps
  real(rk8) , pointer , dimension(:) :: snag
  real(rk8) , pointer , dimension(:) :: sncv
  real(rk8) , pointer , dimension(:) :: sm
  real(rk8) , pointer , dimension(:) :: sst
  real(rk8) , pointer , dimension(:) :: t2m
  real(rk8) , pointer , dimension(:) :: taux
  real(rk8) , pointer , dimension(:) :: tauy
  real(rk8) , pointer , dimension(:) :: tdeltas
  real(rk8) , pointer , dimension(:) :: tgb     ! tground2
  real(rk8) , pointer , dimension(:) :: tgbrd
  real(rk8) , pointer , dimension(:) :: tgrd
  real(rk8) , pointer , dimension(:) :: u10m
  real(rk8) , pointer , dimension(:) :: v10m
  real(rk8) , pointer , dimension(:) :: usw     ! uatm
  real(rk8) , pointer , dimension(:) :: vsw     ! vatm

  real(rk8) , pointer , dimension(:) :: laketa
  real(rk8) , pointer , dimension(:) :: lakhi
  real(rk8) , pointer , dimension(:) :: lakaveice
  real(rk8) , pointer , dimension(:,:) :: laktlake
  integer(ik4) , pointer , dimension(:) :: ilake

  real(rk8) , pointer , dimension(:) :: swdiral
  real(rk8) , pointer , dimension(:) :: lwdiral
  real(rk8) , pointer , dimension(:) :: swdifal
  real(rk8) , pointer , dimension(:) :: lwdifal

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
    call getmem1d(taux,1,nocnp,'ocn_internal:taux')
    call getmem1d(tauy,1,nocnp,'ocn_internal:tauy')
    call getmem1d(tgb,1,nocnp,'ocn_internal:tgb')
    call getmem1d(tgrd,1,nocnp,'ocn_internal:tgrd')
    call getmem1d(tgbrd,1,nocnp,'ocn_internal:tgbrd')
    call getmem1d(u10m,1,nocnp,'ocn_internal:u10m')
    call getmem1d(v10m,1,nocnp,'ocn_internal:v10m')
    call getmem1d(usw,1,nocnp,'ocn_internal:usw')
    call getmem1d(vsw,1,nocnp,'ocn_internal:vsw')
    call getmem1d(mask,1,nocnp,'ocn_internal:mask')
    call getmem1d(icpl,1,nocnp,'ocn_internal:icpl')
    call getmem1d(czenith,1,nocnp,'ocn_internal:czenith')
    call getmem1d(cprate,1,nocnp,'ocn_internal:cprate')
    call getmem1d(ncprate,1,nocnp,'ocn_internal:ncprate')
    call getmem1d(prcp,1,nocnp,'ocn_internal:prcp')
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
