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

      module mod_bdycod

      use mod_dynparam

      implicit none
!
#ifndef BAND
      real(8) , allocatable , dimension(:,:) :: uj1 , uj2 , ujl , ujlx ,&
                      &  vj1 , vj2 , vjl , vjlx
#endif
      real(8) , allocatable , dimension(:,:) :: ps0 , ps1
      real(8) , allocatable , dimension(:,:,:) :: qb0 , qb1 , so0 ,     &
                      & so1 , tb0 , tb1 , ub0 , ub1 , vb0 , vb1
      real(8) , allocatable , dimension(:,:) :: ts0 , ts1
!
#ifndef BAND
      real(8) , allocatable , dimension(:,:) :: peb , pebt , pwb , pwbt
#endif
      real(8) , allocatable , dimension(:,:) :: pnb , pnbt , psbt , pss
#ifndef BAND
      real(8) , allocatable , dimension(:,:,:) :: qeb , qebt , qwb ,    &
           & qwbt , teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , &
           & veb , vebt , vwb , vwbt
#endif
      real(8) , allocatable , dimension(:,:,:) :: qnb , qnbt , qsb ,    &
           & qsbt , tnb , tnbt , tsb , tsbt
      real(8) , allocatable , dimension(:,:) :: ui1 , ui2 , uil , uilx ,&
           & vi1 , vi2 , vil , vilx
      real(8) , allocatable , dimension(:,:,:) :: unb , unbt , usb ,    &
           & usbt , vnb , vnbt , vsb , vsbt

      contains

        subroutine allocate_mod_bdycon
        implicit none
#ifdef MPP1
        allocate(ps0(iy,0:jxp+1))
        allocate(ps1(iy,0:jxp+1))
        allocate(ts0(iy,jxp))
        allocate(ts1(iy,jxp))

        allocate(qb0(iy,kz,jxp))
        allocate(qb1(iy,kz,jxp))
        allocate(so0(iy,kz,jxp))
        allocate(so1(iy,kz,jxp))
        allocate(tb0(iy,kz,jxp))
        allocate(tb1(iy,kz,jxp))
        allocate(ub0(iy,kz,jxp))
        allocate(ub1(iy,kz,jxp))
        allocate(vb0(iy,kz,jxp))
        allocate(vb1(iy,kz,jxp))

#ifndef BAND
        allocate(peb(iy,0:jxp+1))
        allocate(pebt(iy,0:jxp+1))
        allocate(pwb(iy,0:jxp+1))
        allocate(pwbt(iy,0:jxp+1))
#endif

        allocate(pnb(nspgx,0:jxp+1))
        allocate(pnbt(nspgx,0:jxp+1))
        allocate(psbt(nspgx,0:jxp+1))
        allocate(pss(nspgx,0:jxp+1))

#ifndef BAND
        allocate(qeb(iy,kz,0:jxp+1))
        allocate(qebt(iy,kz,0:jxp+1))
        allocate(qwb(iy,kz,0:jxp+1))
        allocate(qwbt(iy,kz,0:jxp+1))
        allocate(teb(iy,kz,0:jxp+1))
        allocate(tebt(iy,kz,0:jxp+1))
        allocate(twb(iy,kz,0:jxp+1))
        allocate(twbt(iy,kz,0:jxp+1))
        allocate(ueb(iy,kz,0:jxp+1))
        allocate(uebt(iy,kz,0:jxp+1))
        allocate(uwb(iy,kz,0:jxp+1))
        allocate(uwbt(iy,kz,0:jxp+1))
        allocate(veb(iy,kz,0:jxp+1))
        allocate(vebt(iy,kz,0:jxp+1))
        allocate(vwb(iy,kz,0:jxp+1))
        allocate(vwbt(iy,kz,0:jxp+1))
#endif

        allocate(qnb(nspgx,kz,0:jxp+1))
        allocate(qnbt(nspgx,kz,0:jxp+1))
        allocate(qsb(nspgx,kz,0:jxp+1))
        allocate(qsbt(nspgx,kz,0:jxp+1))
        allocate(tnb(nspgx,kz,0:jxp+1))
        allocate(tnbt(nspgx,kz,0:jxp+1))
        allocate(tsb(nspgx,kz,0:jxp+1))
        allocate(tsbt(nspgx,kz,0:jxp+1))

        allocate(ui1(kz,0:jxp+1))
        allocate(ui2(kz,0:jxp+1))
        allocate(uil(kz,0:jxp+1))
        allocate(uilx(kz,0:jxp+1))
        allocate(vi1(kz,0:jxp+1))
        allocate(vi2(kz,0:jxp+1))
        allocate(vil(kz,0:jxp+1))
        allocate(vilx(kz,0:jxp+1))

        allocate(unb(nspgd,kz,0:jxp+1))
        allocate(unbt(nspgd,kz,0:jxp+1))
        allocate(usb(nspgd,kz,0:jxp+1))
        allocate(usbt(nspgd,kz,0:jxp+1))
        allocate(vnb(nspgd,kz,0:jxp+1))
        allocate(vnbt(nspgd,kz,0:jxp+1))
        allocate(vsb(nspgd,kz,0:jxp+1))
        allocate(vsbt(nspgd,kz,0:jxp+1))
#else
        allocate(ps0(iy,jx))
        allocate(ps1(iy,jx))
        allocate(ts0(iy,jx))
        allocate(ts1(iy,jx))
        allocate(qb0(iy,kz,jx))
        allocate(qb1(iy,kz,jx))
        allocate(so0(iy,kz,jx))
        allocate(so1(iy,kz,jx))
        allocate(tb0(iy,kz,jx))
        allocate(tb1(iy,kz,jx))
        allocate(ub0(iy,kz,jx))
        allocate(ub1(iy,kz,jx))
        allocate(vb0(iy,kz,jx))
        allocate(vb1(iy,kz,jx))

#ifndef BAND
        allocate(peb(iy,nspgx))
        allocate(pebt(iy,nspgx))
        allocate(pwb(iy,nspgx))
        allocate(pwbt(iy,nspgx))
#endif

        allocate(pnb(nspgx,jx))
        allocate(pnbt(nspgx,jx))
        allocate(psbt(nspgx,jx))
        allocate(pss(nspgx,jx))

#ifndef BAND
        allocate(qeb(iy,kz,nspgx))
        allocate(qebt(iy,kz,nspgx))
        allocate(qwb(iy,kz,nspgx))
        allocate(qwbt(iy,kz,nspgx))
        allocate(teb(iy,kz,nspgx))
        allocate(tebt(iy,kz,nspgx))
        allocate(twb(iy,kz,nspgx))
        allocate(twbt(iy,kz,nspgx))
#endif

        allocate(qnb(nspgx,kz,jx))
        allocate(qnbt(nspgx,kz,jx))
        allocate(qsb(nspgx,kz,jx))
        allocate(qsbt(nspgx,kz,jx))
        allocate(tnb(nspgx,kz,jx))
        allocate(tnbt(nspgx,kz,jx))
        allocate(tsb(nspgx,kz,jx))
        allocate(tsbt(nspgx,kz,jx))

#ifndef BAND
        allocate(ueb(iy,kz,nspgd))
        allocate(uebt(iy,kz,nspgd))
        allocate(uwb(iy,kz,nspgd))
        allocate(uwbt(iy,kz,nspgd))
        allocate(veb(iy,kz,nspgd))
        allocate(vebt(iy,kz,nspgd))
        allocate(vwb(iy,kz,nspgd))
        allocate(vwbt(iy,kz,nspgd))
#endif

        allocate(ui1(kz,jx))
        allocate(ui2(kz,jx))
        allocate(uil(kz,jx))
        allocate(uilx(kz,jx))
        allocate(vi1(kz,jx))
        allocate(vi2(kz,jx))
        allocate(vil(kz,jx))
        allocate(vilx(kz,jx))
        allocate(unb(nspgd,kz,jx))
        allocate(unbt(nspgd,kz,jx))
        allocate(usb(nspgd,kz,jx))
        allocate(usbt(nspgd,kz,jx))
        allocate(vnb(nspgd,kz,jx))
        allocate(vnbt(nspgd,kz,jx))
        allocate(vsb(nspgd,kz,jx))
        allocate(vsbt(nspgd,kz,jx))
#endif 
#ifndef BAND
        allocate(uj1(iy,kz))
        allocate(uj2(iy,kz))
        allocate(ujl(iy,kz))
        allocate(ujlx(iy,kz))
        allocate(vj1(iy,kz))
        allocate(vj2(iy,kz))
        allocate(vjl(iy,kz))
        allocate(vjlx(iy,kz))
#endif

        end subroutine allocate_mod_bdycon

      end module mod_bdycod
