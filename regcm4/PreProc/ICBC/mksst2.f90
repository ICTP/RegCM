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

      subroutine mksst2(tsccm,sst1,sst2,topogm,xlandu,jx,iy,kdate)
      use mod_datewk
      use mod_datenum
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , kdate
      real(4) , dimension(jx,iy) :: sst1 , sst2 , topogm , tsccm ,      &
                                 &  xlandu
      intent (in) iy , jx , kdate , topogm , xlandu
      intent (out) tsccm
      intent (inout) sst1 , sst2
!
! Local variables
!
      integer :: i , j , k , kdate1 , kdate2 , ks , ks1 , ks2 , lat ,   &
               & lon , nday , nmo , nyear
      real(4) :: wt
!
      ks = 427 + 1045
      do k = 427 + 1045 , 1 , -1
        if ( wkday(k)<=kdate ) then
          ks = k
          exit
        end if
      end do
      kdate1 = wkday(ks)
!
      do k = 1 , 427 + 1045
        if ( wkday(k)>kdate ) then
          ks = k
          exit
        end if
      end do
      kdate2 = wkday(ks)
      call finddate_icbc(ks1,kdate1*100)
      call finddate_icbc(ks,kdate*100)
      call finddate_icbc(ks2,kdate2*100)
      wt = float(ks-ks1)/float(ks2-ks1)
 
!     ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
      do lon = 1 , jx
        do lat = 1 , iy
          sst1(lon,lat) = 0.
          sst2(lon,lat) = 0.
        end do
      end do
 
!     ******           READ IN RCM MONTHLY SST DATASET
 100  continue
      read (60,end=200) nday , nmo , nyear , ((sst1(i,j),j=1,iy),i=1,jx)
      if ( nyear*10000+nmo*100+nday/=kdate1 ) go to 100
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
 
!     ******           READ IN RCM MONTHLY SST DATASET
      read (60,end=200) nday , nmo , nyear , ((sst2(i,j),j=1,iy),i=1,jx)
!     PRINT *, 'READING RCM SST DATA:', NMO, NYEAR
 
      do i = 1 , jx
        do j = 1 , iy
          if ( (topogm(i,j)<=1.) .and.                                  &
             & (xlandu(i,j)>13.9 .and. xlandu(i,j)<15.1) .and.          &
             & (sst1(i,j)>-900.0 .and. sst2(i,j)>-900.0) ) tsccm(i,j)   &
             & = (1.-wt)*sst1(i,j) + wt*sst2(i,j)
        end do
      end do
 
      rewind (60)
 
      return
 200  continue
      print * , 'SST file is not the right one'
      stop 12
      end subroutine mksst2
