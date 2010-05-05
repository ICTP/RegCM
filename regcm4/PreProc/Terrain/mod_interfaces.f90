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

      module mod_interfaces
      implicit none

      interface
        subroutine header(myid)
         integer :: myid
        end subroutine
      end interface

      interface
        subroutine output(nunitc,iy,jx,dsinm,clat,clon,plat,plon,       &
                      & iproj,htgrid,htsdgrid,lndout,xlat,xlon,dlat,    &
                      & dlon,xmap,dattyp,dmap,coriol,snowam,igrads,     &
                      & ibigend,kz,sigma,mask,ptop,nsg,truelatl,        &
                      & truelath,grdfac,filout,lsmtyp,sanda,sandb,claya,&
                      & clayb,frac_lnd,nveg,aertyp,texout,frac_tex,ntex,&
                      & lcoarse)
          character(7) :: aertyp
          real(4) :: clat , clon , dsinm , grdfac , plat , plon , ptop ,&
                   & truelath , truelatl
          character(5) :: dattyp
          character(50) :: filout
          integer :: ibigend , igrads , iy , jx , kz , ntex , nunitc ,  &
                   & nveg , nsg
          character(6) :: iproj
          character(4) :: lsmtyp
          real(4) , dimension(iy,jx) :: claya , clayb , coriol , dlat , &
                                  & dlon , dmap , htgrid , htsdgrid ,   &
                                  & lndout , mask , sanda , sandb ,     &
                                  & snowam , texout , xlat , xlon , xmap
          real(4) , dimension(iy,jx,nveg) :: frac_lnd
          real(4) , dimension(iy,jx,ntex) :: frac_tex
          real(4) , dimension(kz+1) :: sigma
          logical :: lcoarse
        end subroutine
      end interface

      interface
        subroutine setup(nunit,iy,jx,ntypec,iproj,ds,clat,clon,igrads,  &
                       & ibyte,filout,filctl)
          real(4) :: clat , clon , ds
          character(256) :: filctl , filout
          integer :: ibyte , igrads , iy , jx , ntypec , nunit
          character(6) :: iproj
        end subroutine setup
      end interface

      interface
        subroutine surf(xlat,xlon,lnduse,iy,jx,incr,dsgrid,lndout,land, &
                    & nrec,h2opct,lsmtyp,sanda,sandb,claya,clayb,       &
                    & frac_lnd,nveg,aertyp,intext,texout,frac_tex,ntex)
          character(7) :: aertyp
          real(4) :: dsgrid , h2opct
          integer :: incr , iy , jx , nrec , ntex , nveg
          character(4) :: lsmtyp
          real(4) , dimension(iy,jx) :: claya , clayb , lndout , sanda ,&
                                  & sandb , texout , xlat , xlon
          real(4) , dimension(iy,jx,nveg) :: frac_lnd
          real(4) , dimension(iy,jx,ntex) :: frac_tex
          integer , dimension(iy,jx) :: intext , lnduse
          real(4) , dimension(iy,jx,2) :: land
        end subroutine surf
      end interface

      interface
        subroutine xyobsll(iy,jx,iproj,clat,clon,plat,plon,truelath)
          real(4) :: clat , clon , plat , plon , truelath
          character(6) :: iproj
          integer :: iy , jx
        end subroutine xyobsll
      end interface

      end module mod_interfaces
