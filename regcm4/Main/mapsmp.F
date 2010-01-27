C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       subroutine mapsmp(fld,fldr,iyy,jxx,ia,ib,iny,ja,jb,jnx,const,
     1                  ichos,c40nam,time)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c                                                                     c
c     this subroutine prints a sample of a two-dimensional data field c
c     on the line printer with 5 significant digits.                  c
c                                                                     c
c     *** note *** the values of fld(i,j) should be limited within    c
c                  1.e30 --- 1.e-30. if the value outside this        c
c                  range is desired, the program should be changed    c
c                  accordingly (in do loop 20).                       c
c                                                                     c
c                                                                     c
c     fld    : a two-dimensional array to hold the data field to be   c
c              sampled and printed. fld could be a horizontal slice,  c
c              fld(i,j), or a vertical slice fld(k,i) or fld(k,j).    c
c
c     fldr   : reverse array of fld; i.e., fldr(j,i)=fld(i,j)
c                                                                     c
c     iyy    : the first dimension of fld.                            c
c              for the horizontal slice, iyy is the dimension in the  c
c                                        y direction.                 c
c              for the vertical slice, iyy is the dimension in the    c
c                                      z direction.                   c
c                                                                     c
c     jxx    : the second dimension of fld.                           c
c              for the horizontal slice, jxx is the dimension in the  c
c                                        x direction.                 c
c              for the vertical slice, jxx is the dimension in either c
c                                      the x or y direction.          c
c                                                                     c
c     ia     : initial sampling point in the first dimension.         c
c                                                                     c
c     ib     : final sampling point in the first dimension.           c
c                                                                     c
c     iny    : sampling interval in the first dimension.              c
c                                                                     c
c     ja     : initial sampling point in the second dimension.        c
c                                                                     c
c     jb     : final sampling point in the second dimension.          c
c                                                                     c
c     jnx    : sampling interval in the second dimension.             c
c                                                                     c
c     const  : constant used to be subtracted from fldr.              c
c                                                                     c
c     ichos > 0 : for horizontal array fld(y,x)                       c
c           < 0 : for vertical cross section fld(z,y) or fld(z,x)     c
c                                                                     c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      integer iyy,jxx,ia,ib,iny,ja,jb,jnx,ichos
      real(kind=8)  const,time
      real(kind=8)  fldr(jxx,iyy),fld(iyy,jxx)
      integer jm(100)
      character c40nam*40,ifmt1*24,ifmt2*24
      integer ksigt
      data ksigt/5/
      integer i,j,k1,k2,k3,k4,n,n1,jx,iy,jn,jn1,j1,jl,jt,jtn
      integer j2n,j2,iexp,jj,it,ir,i1,i2,j3
      real(kind=8)  fldmax,fldmin,fldu,fldl
c
      do i=1,iyy
      do j=1,jxx
        fldr(j,i) = fld(i,j)
      end do
      end do
c
      n=6
      k1=ksigt+2
      k2=124/k1
      k3=ksigt/2
      k4=ksigt-k3
c
      do 10 i=ia,ib,iny
      do 10 j=ja,jb,jnx
10     fldr(j,i)=fldr(j,i)-const
c
      fldmax=fldr(ja,ia)
      fldmin=fldr(ja,ia)
      fldu=10.**ksigt
      fldl=10.**(ksigt-1)
      do 20 i=ia,ib,iny
      do 20 j=ja,jb,jnx
       if(dabs(fldr(j,i)) .gt. 1.e30 .or. dabs(fldr(j,i)) .lt. 1.e-30)
     1  go to 20
       if (dabs(fldr(j,i)) .gt. fldmax) fldmax=dabs(fldr(j,i))
       if (dabs(fldr(j,i)) .lt. fldmin) fldmin=dabs(fldr(j,i))
20    continue
c
      if (fldmax .ne. fldmin) then
       iexp=0
       do 50 n1=1,500
        if(fldmax .ge. fldu) then
         fldmax=fldmax/10.
         iexp=iexp-1
        else if(fldmax .lt. fldl) then
         fldmax=fldmax*10.
         iexp=iexp+1
        else
         goto 60
        end if
50     continue
60     do 70 i=ia,ib,iny
       do 70 j=ja,jb,jnx
70      fldr(j,i)=fldr(j,i)*10.**iexp
       jx=jb-ja+1
       iy=ib-ia+1
       jn=k2*jnx
       jn1=jn-1
       write(n,80) c40nam,iexp,time
80     format(////' this is a list of  ',a40,'  ,scaled by  1.e',i3,5x,
     &       'at time = ',f10.3)
       do 130 j1=ja,jb,jn
        jl=min0(j1+jn1,jb)
        jt=jl-j1+1
        jtn=(jt-1)/jnx+1
        j2n=0
        do 90 j2=1,jt,jnx
         j2n=j2n+1
90       jm(j2n)=j1+j2-1
        write(ifmt1,100)jtn,k4,k3
100     format(5h(/4x,,i2,1h(,i2,5hx,i2,,i2,4hx)/))
        write(n,ifmt1) (jm(jj),jj=1,j2n)
        write(ifmt2,110)jtn,k1
110     format(7h(1x,i2,,i2,1hf,i2,9h.0,2x,i2))
c110    format(1x,i2,1x,i2)
        it=(iy-1)/iny
        ir=iy-it*iny
        do 120 i2=ia,ib,iny
         i1=ib+ia-i2-ir+1
         if (ichos .lt. 0) i1=i2
         write(n,ifmt2) i1,(fldr(j3,i1),j3=j1,jl,jnx),i1
120     continue
130     write(n,ifmt1) (jm(jj),jj=1,j2n)
       do 140 i=ia,ib,iny
       do 140 j=ja,jb,jnx
140     fldr(j,i)=fldr(j,i)/(10.**iexp)+const
      else
       do 160 i=ia,ib,iny
       do 160 j=ja,jb,jnx
160     fldr(j,i)=fldr(j,i)+const
       write(n,170) c40nam,fldmax,time
170    format(/'   all of the values of ',a40,' are equal to ',e15.5,5x,
     1       'at time = ',f10.3)
      end if      !end if(fldmax.ne.fldmin)test
c
c
      return
      end subroutine mapsmp
