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

      subroutine inidust

c  ***********************************************************
c  * description of 12- soil categories                  *****
c  *                                                     *****
c  * i         cat                     sizing            *****
c  * ------------------------------------------------    *****
c  * 1         sand                   coarse             *****
c  * 2         lomay sand             coarse             *****
c  * 3         sand lomay             coarse-medium      *****
c  * 4         silt loma              medium-fine        *****
c  * 5         silt                   medium             *****
c  * 6         loam                   fine               *****
c  * 7         sandy clay loam        coarse-medium-fine *****
c  * 8         silty clay loam        medium             *****
c  * 9         clay loam              medium-fine        *****
c  * 10        sandy clay             coarse-fine        *****
c  * 11        silty clay             medium-fine        *****
c  * 12        clay                   fine               *****
c  ***********************************************************

      implicit none

#     include "parameter.inc"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/trachem.cb"
#     include "include/dust.cb"
#     include "include/param1.cb"

      integer i,j,t,itex,n,ns,nm,nt
      
      real(kind=8)  srel(ix,nsoil,nats)
      real(kind=8)  deldp

      real(kind=8)  bcly(nats),bsnd(nats),bslt(nats)

      real(kind=8)  rhop

      data    bcly/0.00,0.10,0.10,0.15,0.15,0.15,
     &             0.20,0.20,0.30,0.35,0.40,0.50/
      data    bsnd/0.90,0.60,0.60,0.50,0.45,0.35,
     &             0.30,0.30,0.20,0.65,0.60,0.50/
      data    bslt/0.10,0.30,0.30,0.35,0.40,0.50,
     &             0.50,0.50,0.50,0.00,0.00,0.00/
             
      data    rhop/2650.000/
      real(kind=8)  pi, eps
      data    pi/3.1415926535897/,eps/1.0e-7/
      real(kind=8)  mmd(3,12),sigma(3,12),pcent(3,12)
      data    mmd/1000.0,   100.0,   10.0,
     &             690.0,   100.0,   10.0,
     &             520.0,   100.0,    5.0,
     &             520.0,   100.0,    5.0,
     &             520.0,    75.0,    2.5,
     &             520.0,    75.0,    2.5,
     &             210.0,    75.0,    2.5,
     &             210.0,    50.0,    2.5,
     &             125.0,    50.0,    1.0,
     &             100.0,    10.0,    1.0,
     &             100.0,    10.0,    0.5,
     &             100.0,    10.0,    0.5/

      data    sigma/1.6, 1.7, 1.8,
     &              1.6, 1.7, 1.8,
     &              1.6, 1.7, 1.8,
     &              1.6, 1.7, 1.8,
     &              1.6, 1.7, 1.8,
     &              1.6, 1.7, 1.8,
     &              1.7, 1.7, 1.8,
     &              1.7, 1.7, 1.8,
     &              1.7, 1.7, 1.8,
     &              1.8, 1.8, 1.8,
     &              1.8, 1.8, 1.8,
     &              1.8, 1.8, 1.8/

      data    pcent/0.90, 0.10, 0.00,
     &              0.60, 0.30, 0.10,
     &              0.60, 0.30, 0.10,
     &              0.50, 0.35, 0.15,
     &              0.45, 0.40, 0.15,
     &              0.35, 0.50, 0.15,
     &              0.30, 0.50, 0.20,
     &              0.30, 0.50, 0.20,
     &              0.20, 0.50, 0.30,
     &              0.65, 0.00, 0.35,
     &              0.60, 0.00, 0.40,
     &              0.50, 0.00, 0.50/
      real(kind=8)  ss(nsoil), stotal, xk, xl, xm, xn

#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do n=1,nats
      do i=1,ix
        clay2row2(i,n,j) = 0. 
        sand2row2(i,n,j) = 0.
        silt2row2(i,n,j) = 0.
      end do
      end do
      end do
      do n=1,nsoil
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        srel2d(i,j,n) =0.
      end do
      end do
      end do
#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif
      do i= 2, ilxm
        itex = nint( dustsotex(i,j)) 
        if (itex .ne. 0 .and. itex .le. 12) then 

cremember for the moment one texture type per grid cell !

          clay2row2(i,itex,j) = bcly(itex) * 100.0 
          sand2row2(i,itex,j) = bsnd(itex) * 100.0
          silt2row2(i,itex,j) = bslt(itex) * 100.0

        end if   

        do t=1,nats
          sandrow2(i,j) = sandrow2(i,j)+ sand2row2(i,t,j)
          clayrow2(i,j) = clayrow2(i,j)+ clay2row2(i,t,j)
        end do

      end do !i
      end do !j

      dp(1)=0.0001    !cm
      do ns=2,nsoil
        dp(ns)=dp(ns-1) * exp(0.0460517018598807)
        deldp  =dp(ns) - dp(ns-1)
      end do

#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif

        do n=1,nats
        do ns=1,nsoil
        do i=1,ix
          srel(i,ns,n) = 0.  
        end do
        end do
        end do
        do i=2,ilxm
          if (sandrow2(i,j).gt.0.0 .or. clayrow2(i,j).gt.0.0) then
            do nt=1, nats                     !soil types
              do ns=1,nsoil
                ss(ns) = 0.0
              end do
              stotal = 0.0

              if (sand2row2(i,nt,j).gt.0.) then
                do ns=1, nsoil             !soil size segregatoin no
                do nm=1, mode            !soil mode = 3
                  if ((pcent(nm,nt) .gt. eps ) .and.
     %              (sigma(nm,nt) .ne. 0.0)) then
                    xk = pcent(nm,nt)/(sqrt(2.0*pi)*log(sigma(nm,nt)))
                    xl = ((log(dp(ns))-log(mmd(nm,nt)*1.e-4))**2)/
     1                           (2.0*(log(sigma(nm,nt)))**2)
                    xm = xk*exp(-xl)
                  else
                    xm = 0.0
                  end if
                  xn = rhop*(2.0/3.0)*(dp(ns)/2.0)
                  deldp = 0.0460517018598807  ! dp(2)-dp(1) ss(nsoil)
                  ss(ns) = ss(ns)+(xm*deldp/xn)
                end do
                  stotal=stotal+ss(ns)
                end do

                do ns=1,nsoil
                  if (stotal .gt. 0.0) then
                    srel(i,ns,nt)=ss(ns)/stotal !srel(ix,nsoil,nats)
                  end if
                end do
              end if
            end do                    ! soil types
          end if
        end do
         
        do i=2,ilxm
          itex  =   nint( dustsotex(i,j)) 
          if (itex .ge. 1 .and. itex .le. nats) then
            do ns=1,nsoil
              srel2d(i,j,ns) = srel(i,ns,itex)  
            end do
          end if
        end do
            
      end do 

      return
      end  
