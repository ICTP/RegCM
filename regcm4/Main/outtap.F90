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
 
      subroutine outtap

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output to tape or disk for use c
!     in dataflow analyses.                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : ptop
      use mod_iunits
      use mod_main
      use mod_date
      use mod_bats
      use mod_cvaria
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Local variables
!
      real(4) , dimension(jxm2,iym2) :: fout
      real(8) , dimension(jx,iy) :: out1
      integer :: i , j , k , n , nn
      real(8) :: mmpd
!
 
!
!
!---------------------------------------------------------------------
!-----output RegCM-domain variables:
!
      if ( iotyp.eq.2 ) write (iutdat) idatex
!
!     ******  write one time of data to mm4 output file.
 
      do k = kz , 1 , -1
        do i = 1 , iy
          do j = 1 , jx
#ifdef MPP1
            out1(j,i) = ua_io(i,k,j)
#else
            out1(j,i) = ua(i,k,j)
#endif
          end do
        end do
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa_io(i+1,j+1)
#else
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa(i+1,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kz , 1 , -1
        do i = 1 , iy
          do j = 1 , jx
#ifdef MPP1
            out1(j,i) = va_io(i,k,j)
#else
            out1(j,i) = va(i,k,j)
#endif
          end do
        end do
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa_io(i+1,j+1)
#else
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa(i+1,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kz , 1 , -1
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = omega_io(i+1,k,j+1)
#else
            fout(j,i) = omega(i+1,k,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kz , 1 , -1
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = ta_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
            fout(j,i) = ta(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kz , 1 , -1
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = qva_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
            fout(j,i) = qva(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kz , 1 , -1
        do i = 1 , iym2
          do j = 1 , jxm2
#ifdef MPP1
            fout(j,i) = qca_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
            fout(j,i) = qca(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = (psa_io(i+1,j+1)+ptop)*10.
#else
          fout(j,i) = (psa(i+1,j+1)+ptop)*10.
#endif
        end do
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
      mmpd = 24./tapfrq
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = (rainc_io(i+1,j+1)+rainnc_io(i+1,j+1))*mmpd
#else
          fout(j,i) = (rainc(i+1,j+1)+rainnc(i+1,j+1))*mmpd
#endif
        end do
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
!     ****** write out the following surface fields:
!     1.  temp of lower soil layer (17)
!     2.  total soil water in mm h2o (13)
!     3.  accum infiltration (30)
      do i = 1 , iym2
        do j = 1 , jxm2
#ifdef MPP1
          fout(j,i) = tgb2d_io(1,i+1,j+1)
          do n = 2 , nnsg
            fout(j,i) = fout(j,i) + tgb2d_io(n,i+1,j+1)
          end do
#else
          fout(j,i) = tgb2d(1,i+1,j+1)
          do n = 2 , nnsg
            fout(j,i) = fout(j,i) + tgb2d(n,i+1,j+1)
          end do
#endif
          fout(j,i) = fout(j,i)/float(nnsg)
        end do
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
      do i = 1 , iym2
        do j = 1 , jxm2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
#ifdef MPP1
            if ( ocld2d_io(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + swt2d_io(n,i+1,j+1)
              nn = nn + 1
            end if
#else
            if ( ocld2d(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + swt2d(n,i+1,j+1)
              nn = nn + 1
            end if
#endif
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
      do i = 1 , iym2
        do j = 1 , jxm2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
#ifdef MPP1
            if ( ocld2d_io(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + rno2d_io(n,i+1,j+1)
              nn = nn + 1
            end if
#else
            if ( ocld2d(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + rno2d(n,i+1,j+1)
              nn = nn + 1
            end if
#endif
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
!     changes for accum infiltration
   
#ifdef MPP1
      rno2d_io  = 0.0
      rainc_io  = 0.0
      rainnc_io = 0.0
#else
      rno2d  = 0.0
      rainc  = 0.0
      rainnc = 0.0
#endif
 
      print * , 'OUT-history written date = ' , ldatez + xtime/1440.
 
      end subroutine outtap
