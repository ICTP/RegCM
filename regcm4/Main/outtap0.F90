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
 
      subroutine outtap0

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output to tape or disk for use c
!     in dataflow analyses.                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param2
      use mod_param3 , only : sigma
      use mod_main
      use mod_bats
      use mod_grads
      use mod_date , only : mdate0
      use mod_constants , only : rgti
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Local variables
!
      real(4) :: dtb , dtc , dto , dtr
      integer :: i , j , k
      real(4) , dimension(kzp1) :: sp1d
      real(4) , dimension(jxm2,iym2) :: fout
!
!-----output large-domain variables:
!
      open (20,file='output/OUT_HEAD',form='unformatted',recl=(iym2)    &
          & *(jxm2)*ibyte,access='direct')
      do k = 1 , kzp1
        sp1d(k) = sigma(kzp1-k+1)
      end do
      dto = tapfrq
      dtb = batfrq
      dtr = radisp
      dtc = chemfrq
      write (20,rec=1) mdate0 , ibltyp , icup , ipptls , iboudy , iy ,  &
                     & jx , kz , sp1d , dxsp , ptsp , clat , clon ,     &
                     & plat , plon , proj , dto , dtb , dtr , dtc ,     &
                     & iotyp, truelatl, truelath

      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = ht_io(i+1,j+1)*rgti
#else
          fout(j,i) = ht(i+1,j+1)*rgti
#endif
        end do
      end do
      write (20,rec=2) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = htsd_io(i+1,j+1)
#else
          fout(j,i) = htsd(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=3) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = veg2d_io(i+1,j+1)
#else
          fout(j,i) = veg2d(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=4) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = satbrt_io(i+1,j+1)
#else
          fout(j,i) = satbrt(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=5) fout
      do i = 1 , iym2
#ifdef MPP1
        do j = 1 , jxm2
          fout(j,i) = xlat_io(i+1,j+1)
        end do
#else
        do j = 1 , jxm2
          fout(j,i) = xlat(i+1,j+1)
        end do
#endif
      end do
      write (20,rec=6) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = xlong_io(i+1,j+1)
#else
          fout(j,i) = xlong(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=7) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = 1./msfx_io(i+1,j+1)
#else
          fout(j,i) = 1./msfx(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=8) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = 1./msfd_io(i+1,j+1)
#else
          fout(j,i) = 1./msfd(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=9) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = f_io(i+1,j+1)
#else
          fout(j,i) = f(i+1,j+1)
#endif
        end do
      end do
      write (20,rec=10) fout
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          if ( satbrt_io(i+1,j+1).gt.13.5 .and. satbrt_io(i+1,j+1)      &
             & .lt.15.5 ) then
            fout(j,i) = 0.
          else
            fout(j,i) = 2.
          end if
#else
          if ( satbrt(i+1,j+1).gt.13.5 .and. satbrt(i+1,j+1).lt.15.5 )  &
             & then
            fout(j,i) = 0.
          else
            fout(j,i) = 2.
          end if
#endif
        end do
      end do
      write (20,rec=11) fout
 
      close (20)
!
      end subroutine outtap0
