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
 
      subroutine outtap

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output to tape or disk for use c
!     in dataflow analyses.                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use regcm_param
      use param1
      use param2
      use param3
      use iunits
      use main
      use date
      use bats
      use cvaria
      implicit none
!
! Local variables
!
#ifdef MPP1
      real(4) , dimension(mjx-2,ix-2) :: fout
      real(8) , dimension(mjx,ix) :: out1
#else
      real(4) , dimension(jx-2,ix-2) :: fout
      real(8) , dimension(jx,ix) :: out1
#endif
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
 
      do k = kx , 1 , -1
        do i = 1 , ix
#ifdef MPP1
          do j = 1 , mjx
            out1(j,i) = ua_io(i,k,j)
          end do
#else
          do j = 1 , jx
            out1(j,i) = ua(i,k,j)
          end do
#endif
        end do
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa_io(i+1,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa(i+1,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kx , 1 , -1
        do i = 1 , ix
#ifdef MPP1
          do j = 1 , mjx
            out1(j,i) = va_io(i,k,j)
          end do
#else
          do j = 1 , jx
            out1(j,i) = va(i,k,j)
          end do
#endif
        end do
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa_io(i+1,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)+out1(j+1,i+2) &
                      & +out1(j+2,i+2))/psa(i+1,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = omega_io(i+1,k,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = omega(i+1,k,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = ta_io(i+1,k,j+1)/psa_io(i+1,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = ta(i+1,k,j+1)/psa(i+1,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = qva_io(i+1,k,j+1)/psa_io(i+1,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = qva(i+1,k,j+1)/psa(i+1,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do k = kx , 1 , -1
        do i = 1 , ix - 2
#ifdef MPP1
          do j = 1 , mjx - 2
            fout(j,i) = qca_io(i+1,k,j+1)/psa_io(i+1,j+1)
          end do
#else
          do j = 1 , jx - 2
            fout(j,i) = qca(i+1,k,j+1)/psa(i+1,j+1)
          end do
#endif
        end do
        if ( iotyp.eq.1 ) then
          nrcout = nrcout + 1
          write (iutdat,rec=nrcout) fout
        else if ( iotyp.eq.2 ) then
          write (iutdat) fout
        else
        end if
      end do
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fout(j,i) = (psa_io(i+1,j+1)+ptop)*10.
        end do
#else
        do j = 1 , jx - 2
          fout(j,i) = (psa(i+1,j+1)+ptop)*10.
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
      mmpd = 24./tapfrq
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fout(j,i) = (rainc_io(i+1,j+1)+rainnc_io(i+1,j+1))*mmpd
        end do
#else
        do j = 1 , jx - 2
          fout(j,i) = (rainc(i+1,j+1)+rainnc(i+1,j+1))*mmpd
        end do
#endif
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
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fout(j,i) = tgb2d_io(1,i+1,j+1)
          do n = 2 , nnsg
            fout(j,i) = fout(j,i) + tgb2d_io(n,i+1,j+1)
          end do
          fout(j,i) = fout(j,i)/float(nnsg)
        end do
#else
        do j = 1 , jx - 2
          fout(j,i) = tgb2d(1,i+1,j+1)
          do n = 2 , nnsg
            fout(j,i) = fout(j,i) + tgb2d(n,i+1,j+1)
          end do
          fout(j,i) = fout(j,i)/float(nnsg)
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
            if ( ocld2d_io(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + swt2d_io(n,i+1,j+1)
              nn = nn + 1
            end if
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
#else
        do j = 1 , jx - 2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
            if ( ocld2d(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + swt2d(n,i+1,j+1)
              nn = nn + 1
            end if
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
      do i = 1 , ix - 2
#ifdef MPP1
        do j = 1 , mjx - 2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
            if ( ocld2d_io(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + rno2d_io(n,i+1,j+1)
              nn = nn + 1
            end if
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
#else
        do j = 1 , jx - 2
          fout(j,i) = 0.0
          nn = 0
          do n = 1 , nnsg
            if ( ocld2d(n,i+1,j+1).ge.0.5 ) then
              fout(j,i) = fout(j,i) + rno2d(n,i+1,j+1)
              nn = nn + 1
            end if
          end do
          if ( nn.ge.max0(nnsg/2,1) ) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.E34
          end if
        end do
#endif
      end do
      if ( iotyp.eq.1 ) then
        nrcout = nrcout + 1
        write (iutdat,rec=nrcout) fout
      else if ( iotyp.eq.2 ) then
        write (iutdat) fout
      else
      end if
 
!     changes for accum infiltration
 
      do i = 1 , ix - 1
#ifdef MPP1
        do j = 1 , mjx - 1
          do n = 1 , nnsg
            rno2d_io(n,i,j) = 0.
          end do
          rainc_io(i,j) = 0.
          rainnc_io(i,j) = 0.
        end do
#else
        do j = 1 , jx - 1
          do n = 1 , nnsg
            rno2d(n,i,j) = 0.
          end do
          rainc(i,j) = 0.
          rainnc(i,j) = 0.
        end do
#endif
      end do
 
      print * , 'OUT-history written date = ' , ldatez + xtime/1440.
 
      end subroutine outtap
