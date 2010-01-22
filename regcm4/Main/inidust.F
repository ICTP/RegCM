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
c  * Description of 12- soil categories                  *****
c  *                                                     *****
c  * I         cat                     sizing            *****
c  * ------------------------------------------------    *****
c  * 1         Sand                   coarse             *****
c  * 2         Lomay Sand             coarse             *****
c  * 3         Sand Lomay             coarse-medium      *****
c  * 4         Silt Loma              medium-fine        *****
c  * 5         Silt                   medium             *****
c  * 6         Loam                   fine               *****
c  * 7         Sandy Clay Loam        coarse-medium-fine *****
c  * 8         Silty Clay Loam        medium             *****
c  * 9         Clay Loam              medium-fine        *****
c  * 10        Sandy Clay             coarse-fine        *****
c  * 11        Silty Clay             medium-fine        *****
c  * 12        clay                   fine               *****
c  ***********************************************************

      implicit none

#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
#endif
#     include "include/trachem.cb"
#     include "include/dust.cb"
#     include "include/param1.cb"

      integer i,j,t,itex,n,ns,nm,nt
      
      real(kind=8)  SREL(ix,nsoil,nats)
      real(kind=8)  deldp

      real(kind=8)  BCLY(nats),BSND(nats),BSLT(nats)

      real(kind=8)  rhop

      data    BCLY/0.00,0.10,0.10,0.15,0.15,0.15,
     &             0.20,0.20,0.30,0.35,0.40,0.50/
      data    BSND/0.90,0.60,0.60,0.50,0.45,0.35,
     &             0.30,0.30,0.20,0.65,0.60,0.50/
      data    BSLT/0.10,0.30,0.30,0.35,0.40,0.50,
     &             0.50,0.50,0.50,0.00,0.00,0.00/
             
      data    rhop/2650.000/
      real(kind=8)  pi, eps
      data    pi/3.1415926535897/,eps/1.0e-7/
      real(kind=8)  mmd(3,12),sigma(3,12),pcent(3,12)
      DATA    MMD/1000.0,   100.0,   10.0,
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

      DATA    SIGMA/1.6, 1.7, 1.8,
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

      DATA    PCENT/0.90, 0.10, 0.00,
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
      enddo
      enddo
      enddo
      do n=1,nsoil
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        SREL2D(i,j,n) =0.
      enddo
      enddo
      enddo
#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif
      DO i= 2, ilxm
        itex = nint( dustsotex(i,j)) 
        if (itex .ne. 0 .and. itex .le. 12) then 

cremember for the moment one texture type per grid cell !

          clay2row2(i,itex,j) = BCLY(itex) * 100.0 
          sand2row2(i,itex,j) = BSND(itex) * 100.0
          silt2row2(i,itex,j) = BSLT(itex) * 100.0

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
        DELDP  =dp(ns) - dp(ns-1)
      end do

#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif

        do n=1,nats
        do ns=1,nsoil
        do i=1,ix
          SREL(i,ns,n) = 0.  
        enddo
        enddo
        enddo
        do i=2,ilxm
          IF (sandrow2(i,j).gt.0.0 .OR. clayrow2(i,j).gt.0.0) then
            DO NT=1, NATS                     !SOIL TYPES
              do NS=1,NSOIL
                SS(NS) = 0.0
              enddo
              STOTAL = 0.0

              IF(sand2row2(I,NT,j).GT.0.) THEN
                DO NS=1, NSOIL             !SOIL SIZE SEGREGATOIN NO
                DO NM=1, MODE            !SOIL MODE = 3
                  if((PCENT(NM,NT) .gt. eps ) .and.
     %              (SIGMA(NM,NT) .ne. 0.0)) then
                    XK = PCENT(NM,NT)/(SQRT(2.0*PI)*LOG(SIGMA(NM,NT)))
                    XL = ((LOG(DP(NS))-LOG(MMD(NM,NT)*1.E-4))**2)/
     1                           (2.0*(LOG(SIGMA(NM,NT)))**2)
                    XM = XK*EXP(-XL)
                  else
                    xm = 0.0
                  end if
                  XN = RHOP*(2.0/3.0)*(DP(NS)/2.0)
                  DELDP = 0.0460517018598807  ! DP(2)-DP(1) SS(nsoil)
                  SS(NS) = SS(NS)+(XM*DELDP/XN)
                END DO
                  STOTAL=STOTAL+SS(NS)
                END DO

                DO NS=1,NSOIL
                  IF (STOTAL .GT. 0.0) THEN
                    SREL(I,NS,NT)=SS(NS)/STOTAL !srel(ix,nsoil,nats)
                  END IF
                END DO
              END IF
            END DO                    ! soil types
          END IF
        END DO
         
        do i=2,ilxm
          itex  =   nint( dustsotex(i,j)) 
          if (itex .ge. 1 .and. itex .le. nats) then
            do ns=1,nsoil
              SREL2D(i,j,ns) = SREL(i,ns,itex)  
            enddo
          end if
        end do
            
      end do 

      return
      end  
