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
      real(8) ,allocatable, dimension(:,:) :: ps0 , ps1
      real(8) ,allocatable, dimension(:,:,:) :: qb0 , qb1 , so0 , so1 , tb0 ,   &
                                      & tb1 , ub0 , ub1 , vb0 , vb1
      real(8) ,allocatable, dimension(:,:) :: ts0 , ts1
#else
      real(8) , dimension(iy,jx) :: ps0 , ps1 , ts0 , ts1
      real(8) , dimension(iy,kz,jx) :: qb0 , qb1 , so0 , so1 , tb0 ,    &
                                     & tb1 , ub0 , ub1 , vb0 , vb1
#endif
!
#ifdef MPP1
      real(8) ,allocatable, dimension(:,:) :: peb , pebt , pwb , pwbt
      real(8) ,allocatable, dimension(:,:) :: pnb , pnbt , psbt , pss
      real(8) ,allocatable, dimension(:,:,:) :: qeb , qebt , qwb , qwbt ,   &
           & teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , veb ,  &
           & vebt , vwb , vwbt
      real(8) ,allocatable, dimension(:,:,:) :: qnb , qnbt , qsb , qsbt ,&
           & tnb , tnbt , tsb , tsbt
      real(8) ,allocatable, dimension(:,:) :: ui1 , ui2 , uil , uilx , vi1 , &
           & vi2 , vil , vilx
      real(8) , dimension(iy,kz) :: uj1 , uj2 , ujl , ujlx , vj1 , vj2 ,&
                                  & vjl , vjlx
      real(8) ,allocatable, dimension(:,:,:) :: unb , unbt , usb , usbt ,&
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

contains
	subroutine allocate_mod_bdycon

	allocate(ps0(iy,0:jxp+1))
	allocate(ps1(iy,0:jxp+1))
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

	allocate(ts0(iy,jxp))
	allocate(ts1(iy,jxp))

	allocate(peb(iy,0:jxp+1))
	allocate(pebt(iy,0:jxp+1))
	allocate(pwb(iy,0:jxp+1))
	allocate(pwbt(iy,0:jxp+1))

	allocate(pnb(nspgx,0:jxp+1))
	allocate(pnbt(nspgx,0:jxp+1))
	allocate(psbt(nspgx,0:jxp+1))
	allocate(pss(nspgx,0:jxp+1))
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

	end subroutine allocate_mod_bdycon

      end module mod_bdycod
