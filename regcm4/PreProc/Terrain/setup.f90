!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine setup(nunit,iy,jx,ntypec,iproj,ds,clat,clon,igrads,    &
                     & ibyte,filout,filctl)
      use mod_const
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds
      character(50) :: filctl , filout
      integer :: ibyte , igrads , iy , jx , ntypec , nunit
      character(6) :: iproj
      intent (in) clat , clon , ds , ibyte , igrads , iproj , iy , jx , &
                & ntypec , nunit
!
      rin = 1.5          ! 1.5 rad of influence-coarse mesh
 
      write (6,*) 'ntypec=' , ntypec
      write (6,*) 'iy=' , iy
      write (6,*) 'jx=' , jx
      write (6,*) 'ds=' , ds
      write (6,*) 'clat=' , clat
      write (6,*) 'clon=' , clon
      write (6,*) 'rin=' , rin
      write (6,*) 'iproj=' , iproj
!
      call fexist(filout)
      open (nunit,file=filout,status='unknown',form='unformatted',      &
          & access='direct',recl=iy*jx*ibyte)
      if ( igrads==1 ) then
        call fexist(filctl)
        open (31,file=filctl,status='unknown')
      end if
 
!
      dsinm = ds*1000.
!
      nnc = nint(60./float(ntypec))
      xnc = float(ntypec)/60.
      print * , '***** Terrain resolution (min): ' , xnc*60.
!
      end subroutine setup
