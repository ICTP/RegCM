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

      subroutine mksst3(tsccm,sst1,topogm,xlandu,jx,iy,kdate)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , kdate
      real , dimension(jx,iy) :: sst1 , topogm , tsccm , xlandu
      intent (in) iy , jx , kdate , topogm , xlandu
      intent (out) tsccm
      intent (inout) sst1
!
! Local variables
!
      integer :: i , j , mday , mhour , mmo , myear , nday , nhour ,    &
               & nmo , nyear
!
      nyear = kdate/1000000
      nmo = (kdate-nyear*1000000)/10000
      nday = (kdate-nyear*1000000-nmo*10000)/100
      nhour = kdate - nyear*1000000 - nmo*10000 - nday*100
!
!     ******           INITIALIZE SST1, SST2 (NEEDED FOR 82 JAN CASE)
 
!     ******           READ IN RCM 6 HOUR SST DATASET
 100  continue
      read (60,end=200) mhour , mday , mmo , myear ,                    &
                      & ((sst1(i,j),j=1,iy),i=1,jx)
      if ( nyear/=myear .or. nmo/=mmo .or. nday/=mday .or.              &
         & nhour/=mhour ) go to 100
 
      do i = 1 , jx
        do j = 1 , iy
          if ( (topogm(i,j)<=1.) .and.                                  &
             & (xlandu(i,j)>13.9 .and. xlandu(i,j)<15.1) .and.          &
             & (sst1(i,j)>-900.0) ) tsccm(i,j) = sst1(i,j)
        end do
      end do
 
      rewind (60)
 
      return
 200  continue
      print * , 'SST file is not the right one'
      stop 12
      end subroutine mksst3
