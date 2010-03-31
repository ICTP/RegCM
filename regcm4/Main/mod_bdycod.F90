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

      use mod_regcm_param

      implicit none
!
#ifdef MPP1
      real(8) , dimension(iy,0:jxp+1) :: ps0 , ps1
      real(8) , dimension(iy,kz,jxp) :: qb0 , qb1 , so0 , so1 , tb0 ,   &
                                      & tb1 , ub0 , ub1 , vb0 , vb1
      real(8) , dimension(iy,jxp) :: ts0 , ts1
#else
      real(8) , dimension(iy,jx) :: ps0 , ps1 , ts0 , ts1
      real(8) , dimension(iy,kz,jx) :: qb0 , qb1 , so0 , so1 , tb0 ,    &
                                     & tb1 , ub0 , ub1 , vb0 , vb1
#endif
!
#ifdef MPP1
      real(8) , dimension(iy,0:jxp+1) :: peb , pebt , pwb , pwbt
      real(8) , dimension(nspgx,0:jxp+1) :: pnb , pnbt , psbt , pss
      real(8) , dimension(iy,kz,0:jxp+1) :: qeb , qebt , qwb , qwbt ,   &
           & teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , veb ,  &
           & vebt , vwb , vwbt
      real(8) , dimension(nspgx,kz,0:jxp+1) :: qnb , qnbt , qsb , qsbt ,&
           & tnb , tnbt , tsb , tsbt
      real(8) , dimension(kz,0:jxp+1) :: ui1 , ui2 , uil , uilx , vi1 , &
           & vi2 , vil , vilx
      real(8) , dimension(iy,kz) :: uj1 , uj2 , ujl , ujlx , vj1 , vj2 ,&
                                  & vjl , vjlx
      real(8) , dimension(nspgd,kz,0:jxp+1) :: unb , unbt , usb , usbt ,&
           & vnb , vnbt , vsb , vsbt
#else
      real(8) , dimension(iy,nspgx) :: peb , pebt , pwb , pwbt
      real(8) , dimension(nspgx,jx) :: pnb , pnbt , psbt , pss
      real(8) , dimension(iy,kz,nspgx) :: qeb , qebt , qwb , qwbt ,     &
           & teb , tebt , twb , twbt
      real(8) , dimension(nspgx,kz,jx) :: qnb , qnbt , qsb , qsbt ,     &
           & tnb , tnbt , tsb , tsbt
      real(8) , dimension(iy,kz,nspgd) :: ueb , uebt , uwb , uwbt ,     &
           & veb , vebt , vwb , vwbt
      real(8) , dimension(kz,jx) :: ui1 , ui2 , uil , uilx , vi1 , vi2 ,&
                                  & vil , vilx
      real(8) , dimension(iy,kz) :: uj1 , uj2 , ujl , ujlx , vj1 , vj2 ,&
                                  & vjl , vjlx
      real(8) , dimension(nspgd,kz,jx) :: unb , unbt , usb , usbt ,     &
           & vnb , vnbt , vsb , vsbt
#endif

      end module mod_bdycod
