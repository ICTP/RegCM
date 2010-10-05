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

      module mod_main
!
      use mod_dynparam
      use mod_runparams
!
! Storage for all the 3d prognostic variables in two
!     timesteps and all the 2d variables and constants
!
      type domain
        real(8) , allocatable , dimension(:,:) :: ht
        real(8) , allocatable , dimension(:,:) :: satbrt
        real(8) , allocatable , dimension(:,:) :: xlat
        real(8) , allocatable , dimension(:,:) :: xlong
        real(8) , allocatable , dimension(:,:) :: msfx
        real(8) , allocatable , dimension(:,:) :: msfd
        real(8) , allocatable , dimension(:,:) :: f
      end type domain

      type atmstate
        real(8) , allocatable , dimension(:,:,:) :: u
        real(8) , allocatable , dimension(:,:,:) :: v
        real(8) , allocatable , dimension(:,:,:) :: t
        real(8) , allocatable , dimension(:,:,:) :: qv
        real(8) , allocatable , dimension(:,:,:) :: qc
      end type atmstate

      type surfpstate
        real(8) , allocatable , dimension(:,:) :: ps
        real(8) , allocatable , dimension(:,:) :: pdot
      end type surfpstate

      type surftstate
        real(8) , allocatable , dimension(:,:) :: tg
      end type surftstate

      type domfact
        real(8) , allocatable , dimension(:,:) :: hgfact
      end type domfact

      type surfstate
        real(8) , allocatable , dimension(:,:) :: rainc
        real(8) , allocatable , dimension(:,:) :: rainnc
        real(8) , allocatable , dimension(:,:) :: hfx
        real(8) , allocatable , dimension(:,:) :: qfx
        real(8) , allocatable , dimension(:,:) :: tgbb
        real(8) , allocatable , dimension(:,:) :: zpbl
        real(8) , allocatable , dimension(:,:) :: uvdrag
      end type surfstate

      type sulfate
        real(8) , allocatable , dimension(:,:,:) :: so4
      end type sulfate

      type splitsave
        real(8) , allocatable, dimension(:,:,:) :: dstor
        real(8) , allocatable, dimension(:,:,:) :: hstor
      end type splitsave

      type grellwinds
        real(8) , allocatable , dimension(:,:) :: usk
        real(8) , allocatable , dimension(:,:) :: vsk
      end type grellwinds

      type(domain) , public :: mddom
      type(atmstate) , public :: atm1 , atm2
      type(surfpstate) , public :: sps1 , sps2
      type(surftstate) , public :: sts1 , sts2
      type(surfstate) , public :: sfsta
      type(sulfate) , public :: sulf
      type(splitsave) , public :: spsav
      type(domfact) , public :: domfc
      type(grellwinds) , public :: gwnd

      public :: atmstate , domain
      public :: allocate_mod_main , allocate_atmstate , allocate_domain

      private
!
      contains
!
        subroutine allocate_atmstate(atm,lmpi,ib,jb)
          implicit none
          logical , intent(in) :: lmpi
          integer , intent(in) :: ib , jb
          type(atmstate) , intent(out) :: atm
          integer :: is , ie , js , je
          if (lmpi) then
            is = 1-ib
            ie = iy+ib
            js = 1-jb
            je = jxp+jb
            allocate(atm%u(is:ie,kz,js:je))
            allocate(atm%v(is:ie,kz,js:je))
            allocate(atm%t(is:ie,kz,js:je))
            allocate(atm%qv(is:ie,kz,js:je))
            allocate(atm%qc(is:ie,kz,js:je))
          else
            allocate(atm%u(iy,kz,jx))
            allocate(atm%v(iy,kz,jx))
            allocate(atm%t(iy,kz,jx))
            allocate(atm%qv(iy,kz,jx))
            allocate(atm%qc(iy,kz,jx))
          end if
          atm%u  = 0.0D0
          atm%v  = 0.0D0
          atm%t  = 0.0D0
          atm%qv = 0.0D0
          atm%qc = 0.0D0
        end subroutine allocate_atmstate
!
        subroutine allocate_surfpstate(sps,lmpi)
          implicit none
          logical , intent(in) :: lmpi
          type(surfpstate) , intent(out) :: sps
          if (lmpi) then
            allocate(sps%ps(iy,-1:jxp+2))
            allocate(sps%pdot(iy,-1:jxp+2))
          else
            allocate(sps%ps(iy,jx))
            allocate(sps%pdot(iy,jx))
          endif
          sps%ps = 0.0D0
          sps%pdot = 0.0D0
        end subroutine allocate_surfpstate
!
        subroutine allocate_surftstate(sts,lmpi)
          implicit none
          logical , intent(in) :: lmpi
          type(surftstate) , intent(out) :: sts
          if (lmpi) then
            allocate(sts%tg(iy,jxp+1))
          else
            allocate(sts%tg(iy,jx))
          endif
          sts%tg = 0.0D0
        end subroutine allocate_surftstate
!
        subroutine allocate_domain(dom,lmpi)
          implicit none
          logical , intent(in) :: lmpi
          type(domain) , intent(out) :: dom

          if (lmpi) then
            allocate(dom%ht(iy,0:jxp+1))
            allocate(dom%satbrt(iy,jxp+1)) 
            allocate(dom%xlat(iy,jxp))
            allocate(dom%xlong(iy,jxp))
            allocate(dom%msfx(iy,-1:jxp+2))
            allocate(dom%msfd(iy,-1:jxp+2))
            allocate(dom%f(iy,jxp))
          else
            allocate(dom%ht(iy,jx))
            allocate(dom%satbrt(iy,jx)) 
            allocate(dom%xlat(iy,jx))
            allocate(dom%xlong(iy,jx))
            allocate(dom%msfx(iy,jx))
            allocate(dom%msfd(iy,jx))
            allocate(dom%f(iy,jx))
          end if
          dom%ht = 0.0D0
          dom%satbrt = 0.0D0
          dom%xlat = 0.0D0
          dom%xlong = 0.0D0
          dom%msfx = 0.0D0
          dom%msfd = 0.0D0
          dom%f = 0.0D0
        end subroutine allocate_domain
!
        subroutine allocate_domfact(dfa,lmpi)
          implicit none
          logical , intent(in) :: lmpi
          type(domfact) , intent(out) :: dfa
          if (lmpi) then
            allocate(dfa%hgfact(iy,jxp))
          else
            allocate(dfa%hgfact(iy,jx))
          end if
          dfa%hgfact = 0.0D0
        end subroutine allocate_domfact
!
        subroutine allocate_grellwinds
          implicit none
          allocate(gwnd%usk(iy,kz))
          allocate(gwnd%vsk(iy,kz))
          gwnd%usk = 0.0D0
          gwnd%vsk = 0.0D0
        end subroutine allocate_grellwinds

        subroutine allocate_mod_main(lmpi)
        implicit none
        logical , intent(in) :: lmpi

        call allocate_domain(mddom,lmpi)
        call allocate_domfact(domfc,lmpi)

        call allocate_atmstate(atm1,lmpi,0,2)
        call allocate_surfpstate(sps1,lmpi)
        call allocate_surftstate(sts1,lmpi)
        call allocate_atmstate(atm2,lmpi,0,2)
        call allocate_surfpstate(sps2,lmpi)
        call allocate_surftstate(sts2,lmpi)
        if (icup == 2) then
          call allocate_grellwinds
        end if

        if (lmpi) then
          allocate(sfsta%hfx(iy,jxp))
          allocate(sfsta%qfx(iy,jxp))
          allocate(sfsta%rainc(iy,jxp))
          allocate(sfsta%rainnc(iy,jxp))
          allocate(sfsta%tgbb(iy,jxp))
          allocate(sfsta%zpbl(iy,jxp))
          allocate(sfsta%uvdrag(iy,0:jxp))
          allocate(sulf%so4(iy,kz,jxp))
          allocate(spsav%dstor(iy,0:jxp+1,nsplit))
          allocate(spsav%hstor(iy,0:jxp+1,nsplit))
        else
          allocate(sfsta%hfx(iy,jx))
          allocate(sfsta%qfx(iy,jx))
          allocate(sfsta%rainc(iy,jx))
          allocate(sfsta%rainnc(iy,jx))
          allocate(sfsta%tgbb(iy,jx))
          allocate(sfsta%zpbl(iy,jx))
          allocate(sfsta%uvdrag(iy,jx))
          allocate(sulf%so4(iy,kz,jx))
          allocate(spsav%dstor(iy,jx,nsplit))
          allocate(spsav%hstor(iy,jx,nsplit))
        end if

        sfsta%hfx = 0.0D0
        sfsta%qfx = 0.0D0
        sfsta%rainc = 0.0D0
        sfsta%rainnc = 0.0D0
        sfsta%tgbb = 0.0D0
        sfsta%zpbl = 0.0D0
        sfsta%uvdrag = 0.0D0
        sulf%so4 = 0.0D0
        spsav%dstor = 0.0D0
        spsav%hstor = 0.0D0
!
        end subroutine allocate_mod_main 
!
      end module mod_main
