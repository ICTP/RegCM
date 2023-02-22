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

module mod_pbl_interface

  use mod_realkinds
  use mod_service
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_mppparam
  use mod_regcm_types
  use mod_pbl_common , only : ricr , uwstate , kmxpbl
  use mod_pbl_holtbl , only : holtbl , allocate_mod_pbl_holtbl
  use mod_pbl_uwtcm , only : allocate_tcm_state
  use mod_pbl_uwtcm , only : uwtcm , init_mod_pbl_uwtcm , uwtkemin
  use mod_pbl_gfs , only : init_pbl_gfs , pbl_gfs
  use mod_pbl_myj , only : init_myjpbl , myjpbl , myjtkemin
  use mod_runparams , only : ibltyp , pc_physic
  use mod_runparams , only : iqc , iqv , dt , rdt , ichem , hsigma , dsigma

  implicit none

  private

  type(mod_2_pbl) :: m2p
  type(pbl_2_mod) :: p2m

  public :: allocate_pblscheme
  public :: init_pblscheme
  public :: pblscheme

  public :: uwstate
  public :: ricr
  public :: kmxpbl

  real(rkx) , public :: tkemin = 0.0_rkx
  real(rkx) , pointer , dimension(:,:,:) :: utenx , vtenx
  real(rkx) , pointer , dimension(:,:,:) :: utend , vtend

  contains

  subroutine allocate_pblscheme
    use mod_atm_interface
    implicit none
    if ( ibltyp == 1 ) then
      call getmem2d(ricr,jci1,jci2,ici1,ici2,'pbl_common:ricr')
      call getmem2d(kmxpbl,jci1,jci2,ici1,ici2,'pbl_common:kmxpbl')
!$acc enter data create(kmxpbl)
      call allocate_mod_pbl_holtbl
    else if ( ibltyp == 2 ) then
      call allocate_tcm_state(uwstate)
      tkemin = uwtkemin
      call init_mod_pbl_uwtcm
    else if ( ibltyp == 3 ) then
      call init_pbl_gfs
    else if ( ibltyp == 4 ) then
      tkemin = myjtkemin
      call init_myjpbl
    end if
    if ( ibltyp > 1 ) then
      if ( idynamic == 3 ) then
        call getmem3d(utenx,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'pbl_common:utenx')
        call getmem3d(vtenx,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'pbl_common:vtenx')
        call getmem3d(utend,jdi1,jdi2,ici1,ici2,1,kz,'pbl_common:utend')
        call getmem3d(vtend,jci1,jci2,idi1,idi2,1,kz,'pbl_common:vtend')
      else
        call getmem3d(utenx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'pbl_common:utenx')
        call getmem3d(vtenx,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'pbl_common:vtenx')
        call getmem3d(utend,jdi1,jdi2,idi1,idi2,1,kz,'pbl_common:utend')
        call getmem3d(vtend,jdi1,jdi2,idi1,idi2,1,kz,'pbl_common:vtend')
      end if
    end if
  end subroutine allocate_pblscheme

  subroutine init_pblscheme
    use mod_atm_interface
    use mod_che_interface
    use mod_che_common
    implicit none

    ! INPUT to PBL
    call assignpnt(mddom%coriol,m2p%coriol)
    call assignpnt(mddom%ldmsk,m2p%ldmsk)
    call assignpnt(mddom%ht,m2p%ht)
    if ( idynamic < 3 ) then
      call assignpnt(sfs%psb,m2p%psb)
      call assignpnt(sfs%psdotb,m2p%psdotb)
    end if
    call assignpnt(sfs%tgbb,m2p%tg)
    call assignpnt(sfs%qfx,m2p%qfx)
    call assignpnt(sfs%hfx,m2p%hfx)
    call assignpnt(sfs%zo,m2p%zo)
    call assignpnt(sfs%uvdrag,m2p%uvdrag)
    call assignpnt(sfs%ram1,m2p%ram1)
    call assignpnt(sfs%rah1,m2p%rah1)
    call assignpnt(sfs%br,m2p%br)
    call assignpnt(sfs%q2m,m2p%q2m)
    call assignpnt(sfs%u10m,m2p%u10m)
    call assignpnt(sfs%v10m,m2p%v10m)
    call assignpnt(sfs%ustar,m2p%ustar)
    call assignpnt(atms%ubx3d,m2p%uxatm)
    call assignpnt(atms%vbx3d,m2p%vxatm)
    call assignpnt(atms%ubd3d,m2p%udatm)
    call assignpnt(atms%vbd3d,m2p%vdatm)
    call assignpnt(atms%tb3d,m2p%tatm)
    call assignpnt(atms%tv3d,m2p%tvatm)
    call assignpnt(atms%pb3d,m2p%patm)
    call assignpnt(atms%pf3d,m2p%patmf)
    call assignpnt(atms%rhob3d,m2p%rhoatm)
    call assignpnt(atms%qxb3d,m2p%qxatm)
    call assignpnt(atms%chib3d,m2p%chib)
    call assignpnt(atms%th3d,m2p%thatm)
    call assignpnt(atms%za,m2p%za)
    call assignpnt(atms%zq,m2p%zq)
    call assignpnt(atms%dzq,m2p%dzq)
    call assignpnt(atms%rhox2d,m2p%rhox2d)
    if ( ibltyp == 2 ) then
      if ( idynamic == 3 ) then
        call assignpnt(mo_atm%tke,m2p%tkests)
      else
        call assignpnt(atm2%tke,m2p%tkests)
      end if
    else if ( ibltyp == 4 ) then
      call assignpnt(atms%tkepbl,m2p%tkests)
      call assignpnt(sfs%uz0,m2p%uz0)
      call assignpnt(sfs%vz0,m2p%vz0)
      call assignpnt(sfs%thz0,m2p%thz0)
      call assignpnt(sfs%qz0,m2p%qz0)
    end if
    call assignpnt(drydepv,m2p%drydepv)
    call assignpnt(chifxuw,m2p%chifxuw)
    call assignpnt(heatrt,m2p%heatrt)

    ! OUTPUT FROM PBL
    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%tten,p2m%tten)
      call assignpnt(mo_atm%uten,p2m%uten)
      call assignpnt(mo_atm%vten,p2m%vten)
      call assignpnt(mo_atm%qxten,p2m%qxten)
      if ( ibltyp == 2 ) then
        call assignpnt(mo_atm%tketen,p2m%tketen)
      else if ( ibltyp == 4 ) then
        call assignpnt(atms%tkepbl,p2m%tkepbl)
      end if
      call assignpnt(mo_atm%chiten,p2m%chiten)
    else
      call assignpnt(aten%t,p2m%tten,pc_physic)
      call assignpnt(aten%u,p2m%uten,pc_physic)
      call assignpnt(aten%v,p2m%vten,pc_physic)
      call assignpnt(aten%qx,p2m%qxten,pc_physic)
      if ( ibltyp == 2 ) then
        call assignpnt(aten%tke,p2m%tketen,pc_physic)
      else if ( ibltyp == 4 ) then
        call assignpnt(atms%tkepbl,p2m%tkepbl)
      end if
      call assignpnt(aten%chi,p2m%chiten,pc_physic)
    end if
    call assignpnt(utenx,p2m%uxten)
    call assignpnt(vtenx,p2m%vxten)
    call assignpnt(remdrd,p2m%remdrd)
    call assignpnt(zpbl,p2m%zpbl)
    call assignpnt(kpbl,p2m%kpbl)

  end subroutine init_pblscheme

  subroutine pblscheme
    use mod_atm_interface
    implicit none
    integer(ik4) :: i , j , k
    select case ( ibltyp )
      case (1)
        call holtbl(m2p,p2m)
      case (2)
        utenx = d_zero
        vtenx = d_zero
        utend = d_zero
        vtend = d_zero
        call uwtcm(m2p,p2m)
        if ( idynamic == 3 ) then
          call tenxtouvten(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)
              end do
            end do
            do i = ici1 , ici2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)
              end do
            end do
          end do
        else
          call uvcross2dot(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)*m2p%psdotb(j,i)
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)*m2p%psdotb(j,i)
              end do
            end do
          end do
        end if
      case (3)
        utenx = d_zero
        vtenx = d_zero
        utend = d_zero
        vtend = d_zero
        call pbl_gfs(m2p,p2m)
        if ( idynamic == 3 ) then
          call tenxtouvten(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)
              end do
            end do
            do i = ici1 , ici2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)
              end do
            end do
          end do
        else
          call uvcross2dot(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)*m2p%psdotb(j,i)
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)*m2p%psdotb(j,i)
              end do
            end do
          end do
        end if
      case (4)
        utenx = d_zero
        vtenx = d_zero
        utend = d_zero
        vtend = d_zero
        call myjpbl(m2p,p2m)
        if ( idynamic == 3 ) then
          call tenxtouvten(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jci1 , jci2
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)
              end do
            end do
            do i = ici1 , ici2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)
              end do
            end do
          end do
        else
          call uvcross2dot(utenx,vtenx,utend,vtend)
          do k = 1 , kz
            do i = idi1 , idi2
              do j = jdi1 , jdi2
                p2m%uten(j,i,k) = p2m%uten(j,i,k)+utend(j,i,k)*m2p%psdotb(j,i)
                p2m%vten(j,i,k) = p2m%vten(j,i,k)+vtend(j,i,k)*m2p%psdotb(j,i)
              end do
            end do
          end do
        end if
      case default
        return
    end select
  end subroutine pblscheme

end module mod_pbl_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
