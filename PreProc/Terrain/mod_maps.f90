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

      use m_realkinds
      use m_mall
      use m_die

      real(SP) , allocatable , dimension(:,:) :: coriol , dlat , dlon ,  &
                     & dmap , htgrid , lndout , mask , dpth , snowam ,  &
                     & texout , xlat , xlon , xmap
      real(SP) , allocatable , dimension(:,:,:) :: frac_tex

      real(SP) , allocatable , dimension(:,:) :: coriol_s , dlat_s ,     &
                        & dlon_s , dmap_s , htgrid_s , lndout_s ,       &
                        & mask_s , dpth_s , snowam_s , texout_s ,       &
                        & xlat_s , xlon_s , xmap_s
      real(SP) , allocatable , dimension(:,:,:) :: frac_tex_s

      real(SP) , allocatable , dimension(:) :: sigma
      real(DP) :: xn

      contains

      subroutine allocate_grid(iy,jx,kz,ntex)
        implicit none
        integer , intent(in) :: iy,jx,kz,ntex
        integer :: ierr
        allocate(sigma(kz+1), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate sigma',ierr)
        call mall_mci(sigma,'maps')
        allocate(coriol(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate coriol',ierr)
        call mall_mci(coriol,'maps')
        allocate(xlat(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate xlat',ierr)
        call mall_mci(xlat,'maps')
        allocate(xlon(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate xlon',ierr)
        call mall_mci(xlon,'maps')
        allocate(dlat(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate dlat',ierr)
        call mall_mci(dlat,'maps')
        allocate(dlon(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate dlon',ierr)
        call mall_mci(dlon,'maps')
        allocate(xmap(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate xmap',ierr)
        call mall_mci(xmap,'maps')
        allocate(dmap(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate dmap',ierr)
        call mall_mci(dmap,'maps')
        allocate(htgrid(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate htgrid',ierr)
        call mall_mci(htgrid,'maps')
        allocate(dpth(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate dpth',ierr)
        call mall_mci(dpth,'maps')
        allocate(lndout(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate lndout',ierr)
        call mall_mci(lndout,'maps')
        allocate(mask(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate mask',ierr)
        call mall_mci(mask,'maps')
        allocate(snowam(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate snowam',ierr)
        call mall_mci(snowam,'maps')
        allocate(texout(iy,jx), stat=ierr)
        if (ierr /= 0) call die('allocate_grid','allocate texout',ierr)
        call mall_mci(texout,'maps')
        allocate(frac_tex(iy,jx,ntex), stat=ierr)
        if (ierr /= 0) call die('allocate_grid', &
                                'allocate frac_tex',ierr)
        call mall_mci(frac_tex,'maps')
      end subroutine allocate_grid

      subroutine allocate_subgrid(iysg,jxsg,ntex)
        implicit none
        integer , intent(in) :: iysg,jxsg,ntex
        integer :: ierr
        allocate(coriol_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate coriol_s',ierr)
        call mall_mci(coriol_s,'maps')
        allocate(xlat_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate xlat_s',ierr)
        call mall_mci(xlat_s,'maps')
        allocate(xlon_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate xlon_s',ierr)
        call mall_mci(xlon_s,'maps')
        allocate(dlat_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate dlat_s',ierr)
        call mall_mci(dlat_s,'maps')
        allocate(dlon_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate dlon_s',ierr)
        call mall_mci(dlon_s,'maps')
        allocate(xmap_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate xmap_s',ierr)
        call mall_mci(xmap_s,'maps')
        allocate(dmap_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate dmap_s',ierr)
        call mall_mci(dmap_s,'maps')
        allocate(htgrid_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate htgrid_s',ierr)
        call mall_mci(htgrid_s,'maps')
        allocate(dpth_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate dpth_s',ierr)
        call mall_mci(dpth_s,'maps')
        allocate(lndout_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate lndout_s',ierr)
        call mall_mci(lndout_s,'maps')
        allocate(mask_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate mask_s',ierr)
        call mall_mci(mask_s,'maps')
        allocate(snowam_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate snowam_s',ierr)
        call mall_mci(snowam_s,'maps')
        allocate(texout_s(iysg,jxsg), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate texout_s',ierr)
        call mall_mci(texout_s,'maps')
        allocate(frac_tex_s(iysg,jxsg,ntex), stat=ierr)
        if (ierr /= 0) call die('allocate_subgrid', &
                                'allocate frac_tex_s',ierr)
        call mall_mci(frac_tex_s,'maps')
      end subroutine allocate_subgrid

      subroutine free_grid
        implicit none
        call mall_mco(sigma,'maps')
        deallocate(sigma)
        call mall_mco(coriol,'maps')
        deallocate(coriol)
        call mall_mco(xlat,'maps')
        deallocate(xlat)
        call mall_mco(xlon,'maps')
        deallocate(xlon)
        call mall_mco(dlat,'maps')
        deallocate(dlat)
        call mall_mco(dlon,'maps')
        deallocate(dlon)
        call mall_mco(xmap,'maps')
        deallocate(xmap)
        call mall_mco(dmap,'maps')
        deallocate(dmap)
        call mall_mco(htgrid,'maps')
        deallocate(htgrid)
        call mall_mco(dpth,'maps')
        deallocate(dpth)
        call mall_mco(lndout,'maps')
        deallocate(lndout)
        call mall_mco(mask,'maps')
        deallocate(mask)
        call mall_mco(snowam,'maps')
        deallocate(snowam)
        call mall_mco(texout,'maps')
        deallocate(texout)
        call mall_mco(frac_tex,'maps')
        deallocate(frac_tex)
      end subroutine free_grid

      subroutine free_subgrid
        implicit none
        call mall_mco(coriol_s,'maps')
        deallocate(coriol_s)
        call mall_mco(xlat_s,'maps')
        deallocate(xlat_s)
        call mall_mco(xlon_s,'maps')
        deallocate(xlon_s)
        call mall_mco(dlat_s,'maps')
        deallocate(dlat_s)
        call mall_mco(dlon_s,'maps')
        deallocate(dlon_s)
        call mall_mco(xmap_s,'maps')
        deallocate(xmap_s)
        call mall_mco(dmap_s,'maps')
        deallocate(dmap_s)
        call mall_mco(htgrid_s,'maps')
        deallocate(htgrid_s)
        call mall_mco(dpth_s,'maps')
        deallocate(dpth_s)
        call mall_mco(lndout_s,'maps')
        deallocate(lndout_s)
        call mall_mco(mask_s,'maps')
        deallocate(mask_s)
        call mall_mco(snowam_s,'maps')
        deallocate(snowam_s)
        call mall_mco(texout_s,'maps')
        deallocate(texout_s)
        call mall_mco(frac_tex_s,'maps')
        deallocate(frac_tex_s)
      end subroutine free_subgrid

      end module mod_maps
