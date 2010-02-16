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
 
      subroutine mapsmp(fld,fldr,iyy,jxx,ia,ib,iny,ja,jb,jnx,const,     &
                      & ichos,c40nam,time)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!                                                                     c
!     this subroutine prints a sample of a two-dimensional data field c
!     on the line printer with 5 significant digits.                  c
!                                                                     c
!     *** note *** the values of fld(i,j) should be limited within    c
!                  1.e30 --- 1.e-30. if the value outside this        c
!                  range is desired, the program should be changed    c
!                  accordingly (in do loop 20).                       c
!                                                                     c
!                                                                     c
!     fld    : a two-dimensional array to hold the data field to be   c
!              sampled and printed. fld could be a horizontal slice,  c
!              fld(i,j), or a vertical slice fld(k,i) or fld(k,j).    c
!
!     fldr   : reverse array of fld; i.e., fldr(j,i)=fld(i,j)
!                                                                     c
!     iyy    : the first dimension of fld.                            c
!              for the horizontal slice, iyy is the dimension in the  c
!                                        y direction.                 c
!              for the vertical slice, iyy is the dimension in the    c
!                                      z direction.                   c
!                                                                     c
!     jxx    : the second dimension of fld.                           c
!              for the horizontal slice, jxx is the dimension in the  c
!                                        x direction.                 c
!              for the vertical slice, jxx is the dimension in either c
!                                      the x or y direction.          c
!                                                                     c
!     ia     : initial sampling point in the first dimension.         c
!                                                                     c
!     ib     : final sampling point in the first dimension.           c
!                                                                     c
!     iny    : sampling interval in the first dimension.              c
!                                                                     c
!     ja     : initial sampling point in the second dimension.        c
!                                                                     c
!     jb     : final sampling point in the second dimension.          c
!                                                                     c
!     jnx    : sampling interval in the second dimension.             c
!                                                                     c
!     const  : constant used to be subtracted from fldr.              c
!                                                                     c
!     ichos > 0 : for horizontal array fld(y,x)                       c
!           < 0 : for vertical cross section fld(z,y) or fld(z,x)     c
!                                                                     c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      implicit none
!
! Dummy arguments
!
      character(40) :: c40nam
      real(8) :: const , time
      integer :: ia , ib , ichos , iny , iyy , ja , jb , jnx , jxx
      real(8) , dimension(iyy,jxx) :: fld
      real(8) , dimension(jxx,iyy) :: fldr
      intent (in) c40nam , const , fld , ia , ib , ichos , iny , iyy ,  &
                & ja , jb , jnx , jxx , time
      intent (inout) fldr
!
! Local variables
!
      real(8) :: fldl , fldmax , fldmin , fldu
      integer :: i , i1 , i2 , iexp , ir , it , iy , j , j1 , j2 , j2n ,&
               & j3 , jj , jl , jn , jn1 , jt , jtn , k1 , k2 , k3 ,    &
               & k4 , ksigt , n , n1
      character(24) :: ifmt1 , ifmt2
      integer , dimension(100) :: jm
!
      data ksigt/5/
!
      do i = 1 , iyy
        do j = 1 , jxx
          fldr(j,i) = fld(i,j)
        end do
      end do
!
      n = 6
      k1 = ksigt + 2
      k2 = 124/k1
      k3 = ksigt/2
      k4 = ksigt - k3
!
      do i = ia , ib , iny
        do j = ja , jb , jnx
          fldr(j,i) = fldr(j,i) - const
        end do
      end do
!
      fldmax = fldr(ja,ia)
      fldmin = fldr(ja,ia)
      fldu = 10.**ksigt
      fldl = 10.**(ksigt-1)
      do i = ia , ib , iny
        do j = ja , jb , jnx
          if ( dabs(fldr(j,i)).le.1.E30 .and. dabs(fldr(j,i))           &
             & .ge.1.E-30 ) then
            if ( dabs(fldr(j,i)).gt.fldmax ) fldmax = dabs(fldr(j,i))
            if ( dabs(fldr(j,i)).lt.fldmin ) fldmin = dabs(fldr(j,i))
          end if
        end do
      end do
!
      if ( fldmax.ne.fldmin ) then
        iexp = 0
        do n1 = 1 , 500
          if ( fldmax.ge.fldu ) then
            fldmax = fldmax/10.
            iexp = iexp - 1
          else if ( fldmax.lt.fldl ) then
            fldmax = fldmax*10.
            iexp = iexp + 1
          else
            exit
          end if
        end do
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i)*10.**iexp
          end do
        end do
        iy = ib - ia + 1
        jn = k2*jnx
        jn1 = jn - 1
        write (n,99001) c40nam , iexp , time
        do j1 = ja , jb , jn
          jl = min0(j1+jn1,jb)
          jt = jl - j1 + 1
          jtn = (jt-1)/jnx + 1
          j2n = 0
          do j2 = 1 , jt , jnx
            j2n = j2n + 1
            jm(j2n) = j1 + j2 - 1
          end do
          write (ifmt1,99002) jtn , k4 , k3
          write (n,ifmt1) (jm(jj),jj=1,j2n)
          write (ifmt2,99003) jtn , k1
!110      format(1x,i2,1x,i2)
          it = (iy-1)/iny
          ir = iy - it*iny
          do i2 = ia , ib , iny
            i1 = ib + ia - i2 - ir + 1
            if ( ichos.lt.0 ) i1 = i2
            write (n,ifmt2) i1 , (fldr(j3,i1),j3=j1,jl,jnx) , i1
          end do
          write (n,ifmt1) (jm(jj),jj=1,j2n)
        end do
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i)/(10.**iexp) + const
          end do
        end do
      else
        do i = ia , ib , iny
          do j = ja , jb , jnx
            fldr(j,i) = fldr(j,i) + const
          end do
        end do
        write (n,99004) c40nam , fldmax , time
      end if      !end if(fldmax.ne.fldmin)test
99001 format (////' this is a list of  ',a40,'  ,scaled by  1.e',i3,5x, &
             &'at time = ',f10.3)
99002 format ('(/4x,',i2,'(',i2,'x,i2,',i2,'x)/)')
99003 format ('(1x,i2,',i2,'f',i2,'.0,2x,i2)')
99004 format (/'   all of the values of ',a40,' are equal to ',e15.5,5x,&
             &'at time = ',f10.3)
!
      end subroutine mapsmp
