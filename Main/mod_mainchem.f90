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

      module mod_mainchem

      use mod_runparams

      implicit none
!
      real(8) , allocatable, dimension(:,:,:,:) :: chemsrc
      real(8) , allocatable, dimension(:,:,:,:) :: chia , chib
      real(8) , allocatable, dimension(:,:,:) :: srclp2
!
      real(8) , allocatable , dimension(:,:,:) :: ddsfc , dtrace ,      &
                                       & wdcvc , wdlsc , wxaq , wxsg

      contains 

        subroutine allocate_mod_mainchem(lmpi)
        implicit none
        logical , intent(in) :: lmpi

        if ( ichem == 1 ) then
          if (lmpi) then
            allocate(chemsrc(iy,jxp,mpy,ntr))
            allocate(chia(iy,kz,-1:jxp+2,ntr))
            allocate(chib(iy,kz,-1:jxp+2,ntr))
            allocate(srclp2(iy,jxp,ntr))
            allocate(ddsfc(iy,jxp,ntr)) 
            allocate(dtrace(iy,jxp,ntr))
            allocate(wdcvc(iy,jxp,ntr))
            allocate(wdlsc(iy,jxp,ntr))
            allocate(wxaq(iy,jxp,ntr))
            allocate(wxsg(iy,jxp,ntr))
          else
            allocate(chemsrc(iy,jx,mpy,ntr))
            allocate(chia(iy,kz,jx,ntr))
            allocate(chib(iy,kz,jx,ntr))
            allocate(srclp2(iy,jx,ntr))
            allocate(ddsfc(iy,jx,ntr))
            allocate(dtrace(iy,jx,ntr))
            allocate(wdcvc(iy,jx,ntr))
            allocate(wdlsc(iy,jx,ntr))
            allocate(wxaq(iy,jx,ntr))
            allocate(wxsg(iy,jx,ntr))
          end if
          chemsrc = d_zero
          chia = d_zero
          chib = d_zero
          srclp2 = d_zero
          ddsfc = d_zero
          dtrace = d_zero
          wdcvc = d_zero
          wdlsc = d_zero
          wxaq = d_zero
          wxsg = d_zero
        end if
       end subroutine allocate_mod_mainchem
!
      end module mod_mainchem
