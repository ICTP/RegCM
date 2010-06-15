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

      module mod_maps

      implicit none

      real(4) , allocatable , dimension(:,:) :: coriol , dlat , dlon ,  &
                     & dmap , htgrid , htsdgrid , lndout , mask ,       &
                     & snowam , texout , xlat , xlon , xmap
      real(4) , allocatable , dimension(:,:,:) :: frac_tex
      integer , allocatable , dimension(:,:) :: intext , lnduse , nsc
      real(4) , allocatable , dimension(:,:) :: corc , hscr1 , htsavc , &
                     & sumc , wtmaxc
      character(1) , allocatable , dimension(:,:) :: ch
      real(4) , allocatable , dimension(:,:,:) :: itex , land

      real(4) , allocatable , dimension(:,:) :: coriol_s , dlat_s ,     &
           & dlon_s , dmap_s , htgrid_s , htsdgrid_s , lndout_s ,       &
           & mask_s , snowam_s , texout_s , xlat_s , xlon_s , xmap_s
      real(4) , allocatable , dimension(:,:,:) :: frac_tex_s
      integer , allocatable , dimension(:,:) :: intext_s , lnduse_s ,   &
                     & nsc_s
      real(4) , allocatable , dimension(:,:) :: corc_s , hscr1_s ,      &
                     & htsavc_s , sumc_s , wtmaxc_s
      character(1) , allocatable , dimension(:,:) :: ch_s
      real(4) , allocatable , dimension(:,:,:) :: itex_s , land_s

      real(4) , allocatable , dimension(:) :: sigma
      real(4) :: xn

      contains

      subroutine allocate_grid(iy,jx,kz,ntex)
        implicit none
        integer , intent(in) :: iy,jx,kz,ntex
        allocate(sigma(kz+1))
        allocate(ch(iy,jx))
        allocate(corc(iy,jx))
        allocate(nsc(iy,jx))
        allocate(hscr1(iy,jx))
        allocate(htsavc(iy,jx))
        allocate(sumc(iy,jx))
        allocate(wtmaxc(iy,jx))
        allocate(itex(iy,jx,2))
        allocate(land(iy,jx,2))
        allocate(coriol(iy,jx))
        allocate(xlat(iy,jx))
        allocate(xlon(iy,jx))
        allocate(dlat(iy,jx))
        allocate(dlon(iy,jx))
        allocate(xmap(iy,jx))
        allocate(dmap(iy,jx))
        allocate(htgrid(iy,jx))
        allocate(htsdgrid(iy,jx))
        allocate(lndout(iy,jx))
        allocate(mask(iy,jx))
        allocate(snowam(iy,jx))
        allocate(texout(iy,jx))
        allocate(intext(iy,jx))
        allocate(lnduse(iy,jx))
        allocate(frac_tex(iy,jx,ntex))
      end subroutine allocate_grid

      subroutine allocate_subgrid(iysg,jxsg,ntex)
        implicit none
        integer , intent(in) :: iysg,jxsg,ntex
        allocate(ch_s(iysg,jxsg))
        allocate(corc_s(iysg,jxsg))
        allocate(nsc_s(iysg,jxsg))
        allocate(hscr1_s(iysg,jxsg))
        allocate(htsavc_s(iysg,jxsg))
        allocate(sumc_s(iysg,jxsg))
        allocate(wtmaxc_s(iysg,jxsg))
        allocate(itex_s(iysg,jxsg,2))
        allocate(land_s(iysg,jxsg,2))
        allocate(coriol_s(iysg,jxsg))
        allocate(xlat_s(iysg,jxsg))
        allocate(xlon_s(iysg,jxsg))
        allocate(dlat_s(iysg,jxsg))
        allocate(dlon_s(iysg,jxsg))
        allocate(xmap_s(iysg,jxsg))
        allocate(dmap_s(iysg,jxsg))
        allocate(htgrid_s(iysg,jxsg))
        allocate(htsdgrid_s(iysg,jxsg))
        allocate(lndout_s(iysg,jxsg))
        allocate(mask_s(iysg,jxsg))
        allocate(snowam_s(iysg,jxsg))
        allocate(texout_s(iysg,jxsg))
        allocate(intext_s(iysg,jxsg))
        allocate(lnduse_s(iysg,jxsg))
        allocate(frac_tex_s(iysg,jxsg,ntex))
      end subroutine allocate_subgrid

      subroutine free_grid
        implicit none
        deallocate(sigma)
        deallocate(ch)
        deallocate(corc)
        deallocate(nsc)
        deallocate(hscr1)
        deallocate(htsavc)
        deallocate(sumc)
        deallocate(wtmaxc)
        deallocate(itex)
        deallocate(land)
        deallocate(coriol)
        deallocate(xlat)
        deallocate(xlon)
        deallocate(dlat)
        deallocate(dlon)
        deallocate(xmap)
        deallocate(dmap)
        deallocate(htgrid)
        deallocate(htsdgrid)
        deallocate(lndout)
        deallocate(mask)
        deallocate(snowam)
        deallocate(texout)
        deallocate(intext)
        deallocate(lnduse)
        deallocate(frac_tex)
      end subroutine free_grid

      subroutine free_subgrid
        implicit none
        deallocate(ch_s)
        deallocate(corc_s)
        deallocate(nsc_s)
        deallocate(hscr1_s)
        deallocate(htsavc_s)
        deallocate(sumc_s)
        deallocate(wtmaxc_s)
        deallocate(itex_s)
        deallocate(land_s)
        deallocate(coriol_s)
        deallocate(xlat_s)
        deallocate(xlon_s)
        deallocate(dlat_s)
        deallocate(dlon_s)
        deallocate(xmap_s)
        deallocate(dmap_s)
        deallocate(htgrid_s)
        deallocate(htsdgrid_s)
        deallocate(lndout_s)
        deallocate(mask_s)
        deallocate(snowam_s)
        deallocate(texout_s)
        deallocate(intext_s)
        deallocate(lnduse_s)
        deallocate(frac_tex_s)
      end subroutine free_subgrid

      end module mod_maps
