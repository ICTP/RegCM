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

      module mod_dust

      use mod_dynparam

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nsoil = 152
      integer , parameter :: nats = 12
      integer , parameter :: mode = 3
      integer , parameter :: jsoilm = 2 
      integer , parameter :: jfs =1 
      integer , parameter :: ust = 1
!
      real(8) , allocatable, dimension(:,:,:) :: clay2row2 , sand2row2 ,&
           & silt2row2
      real(8) ,allocatable,  dimension(:,:) :: clayrow2 , dustsotex ,   &
           & sandrow2
      real(8) ,allocatable,  dimension(:,:,:) :: srel2d
      real(8) , dimension(nsoil) :: dp

      contains 

        subroutine allocate_mod_dust 
        implicit none
#ifdef MPP1
        allocate(clay2row2(iy,nats,jxp))
        allocate(sand2row2(iy,nats,jxp))
        allocate(silt2row2(iy,nats,jxp))
        allocate(clayrow2(iy,jxp))
        allocate(dustsotex(iy,jxp))
        allocate(sandrow2(iy,jxp))
        allocate(srel2d(iy,jxp,nsoil))
#else        
        allocate(clay2row2(iy,nats,jx))
        allocate(sand2row2(iy,nats,jx))
        allocate(silt2row2(iy,nats,jx))
        allocate(clayrow2(iy,jx))
        allocate(dustsotex(iy,jx))
        allocate(sandrow2(iy,jx))
        allocate(srel2d(iy,jx,nsoil))
#endif 

        end subroutine allocate_mod_dust

      subroutine inidust
 
!  ***********************************************************
!  * description of 12- soil categories                  *****
!  *                                                     *****
!  * i         cat                     sizing            *****
!  * ------------------------------------------------    *****
!  * 1         sand                   coarse             *****
!  * 2         lomay sand             coarse             *****
!  * 3         sand lomay             coarse-medium      *****
!  * 4         silt loma              medium-fine        *****
!  * 5         silt                   medium             *****
!  * 6         loam                   fine               *****
!  * 7         sandy clay loam        coarse-medium-fine *****
!  * 8         silty clay loam        medium             *****
!  * 9         clay loam              medium-fine        *****
!  * 10        sandy clay             coarse-fine        *****
!  * 11        silty clay             medium-fine        *****
!  * 12        clay                   fine               *****
!  ***********************************************************
 
      use mod_constants , only : twopi
      implicit none
!
! Local variables
!
      real(8) , dimension(nats) :: bcly , bslt , bsnd
      real(8) :: deldp , eps , rhop , stotal , xk , xl , xm , xn
      integer :: i , itex , j , n , nm , ns , nt , t
      real(8) , dimension(3,12) :: mmd , pcent , sigma
      real(8) , dimension(iy,nsoil,nats) :: srel
      real(8) , dimension(nsoil) :: ss
!
      data bcly/0.00 , 0.10 , 0.10 , 0.15 , 0.15 , 0.15 , 0.20 , 0.20 , &
         & 0.30 , 0.35 , 0.40 , 0.50/
      data bsnd/0.90 , 0.60 , 0.60 , 0.50 , 0.45 , 0.35 , 0.30 , 0.30 , &
         & 0.20 , 0.65 , 0.60 , 0.50/
      data bslt/0.10 , 0.30 , 0.30 , 0.35 , 0.40 , 0.50 , 0.50 , 0.50 , &
         & 0.50 , 0.00 , 0.00 , 0.00/
 
      data rhop/2650.000/
      data eps/1.0E-7/
      data mmd/1000.0 , 100.0 , 10.0 , 690.0 , 100.0 , 10.0 , 520.0 ,   &
         & 100.0 , 5.0 , 520.0 , 100.0 , 5.0 , 520.0 , 75.0 , 2.5 ,     &
         & 520.0 , 75.0 , 2.5 , 210.0 , 75.0 , 2.5 , 210.0 , 50.0 ,     &
         & 2.5 , 125.0 , 50.0 , 1.0 , 100.0 , 10.0 , 1.0 , 100.0 ,      &
         & 10.0 , 0.5 , 100.0 , 10.0 , 0.5/
 
      data sigma/1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 ,  &
         & 1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 , 1.7 ,  &
         & 1.7 , 1.8 , 1.7 , 1.7 , 1.8 , 1.7 , 1.7 , 1.8 , 1.8 , 1.8 ,  &
         & 1.8 , 1.8 , 1.8 , 1.8 , 1.8 , 1.8 , 1.8/
 
      data pcent/0.90 , 0.10 , 0.00 , 0.60 , 0.30 , 0.10 , 0.60 , 0.30 ,&
         & 0.10 , 0.50 , 0.35 , 0.15 , 0.45 , 0.40 , 0.15 , 0.35 ,      &
         & 0.50 , 0.15 , 0.30 , 0.50 , 0.20 , 0.30 , 0.50 , 0.20 ,      &
         & 0.20 , 0.50 , 0.30 , 0.65 , 0.00 , 0.35 , 0.60 , 0.00 ,      &
         & 0.40 , 0.50 , 0.00 , 0.50/
 
      clay2row2 = 0.0
      sand2row2 = 0.0
      silt2row2 = 0.0
      srel2d    = 0.0

#ifdef MPP1
      do j = jbegin , jendm
        do i = 2 , iym2
          itex = nint(dustsotex(i,j))
          if ( itex.ge.1 .and. itex.le.nats ) then
!           remember for the moment one texture type per grid cell !
            clay2row2(i,itex,j) = bcly(itex)*100.0
            sand2row2(i,itex,j) = bsnd(itex)*100.0
            silt2row2(i,itex,j) = bslt(itex)*100.0
          end if
          do t = 1 , nats
            sandrow2(i,j) = sandrow2(i,j) + sand2row2(i,t,j)
            clayrow2(i,j) = clayrow2(i,j) + clay2row2(i,t,j)
          end do
        end do
      end do
#else
      do j = 2 , jxm2
        do i = 2 , iym2
          itex = nint(dustsotex(i,j))
          if ( itex.ge.1 .and. itex.le.nats ) then
!           remember for the moment one texture type per grid cell !
            clay2row2(i,itex,j) = bcly(itex)*100.0
            sand2row2(i,itex,j) = bsnd(itex)*100.0
            silt2row2(i,itex,j) = bslt(itex)*100.0
          end if
          do t = 1 , nats
            sandrow2(i,j) = sandrow2(i,j) + sand2row2(i,t,j)
            clayrow2(i,j) = clayrow2(i,j) + clay2row2(i,t,j)
          end do
        end do
      end do
#endif
 
      dp(1) = 0.0001  !cm
      do ns = 2 , nsoil
        dp(ns) = dp(ns-1)*exp(0.0460517018598807)
        deldp = dp(ns) - dp(ns-1)
      end do
 
#ifdef MPP1
      do j = jbegin , jendm
        do n = 1 , nats
          do ns = 1 , nsoil
            do i = 1 , iy
              srel(i,ns,n) = 0.
            end do
          end do
        end do
        do i = 2 , iym2
          if ( sandrow2(i,j).gt.0.0 .or. clayrow2(i,j).gt.0.0 ) then
            do nt = 1 , nats                  !soil types
              do ns = 1 , nsoil
                ss(ns) = 0.0
              end do
              stotal = 0.0
              if ( sand2row2(i,nt,j).gt.0. ) then
                do ns = 1 , nsoil          !soil size segregatoin no
                  do nm = 1 , mode       !soil mode = 3
                    if ( (pcent(nm,nt).gt.eps) .and.                    &
                       & (sigma(nm,nt).ne.0.0) ) then
                      xk = pcent(nm,nt)/(sqrt(twopi)*log(sigma(nm,nt)))
                      xl = ((log(dp(ns))-log(mmd(nm,nt)*1.E-4))**2)     &
                         & /(2.0*(log(sigma(nm,nt)))**2)
                      xm = xk*exp(-xl)
                    else
                      xm = 0.0
                    end if
                    xn = rhop*(2.0/3.0)*(dp(ns)/2.0)
                    deldp = 0.0460517018598807
                                              ! dp(2)-dp(1) ss(nsoil)
                    ss(ns) = ss(ns) + (xm*deldp/xn)
                  end do
                  stotal = stotal + ss(ns)
                end do
                do ns = 1 , nsoil
                  if ( stotal.gt.0.0 ) srel(i,ns,nt) = ss(ns)/stotal
                                                !srel(iy,nsoil,nats)
                end do
              end if
            end do                    ! soil types
          end if
        end do
        do i = 2 , iym2
          itex = nint(dustsotex(i,j))
          if ( itex.ge.1 .and. itex.le.nats ) then
            do ns = 1 , nsoil
              srel2d(i,j,ns) = srel(i,ns,itex)
            end do
          end if
        end do
      end do
#else 
      do j = 2 , jxm2
        do n = 1 , nats
          do ns = 1 , nsoil
            do i = 1 , iy
              srel(i,ns,n) = 0.
            end do
          end do
        end do
        do i = 2 , iym2
          if ( sandrow2(i,j).gt.0.0 .or. clayrow2(i,j).gt.0.0 ) then
            do nt = 1 , nats                  !soil types
              do ns = 1 , nsoil
                ss(ns) = 0.0
              end do
              stotal = 0.0
              if ( sand2row2(i,nt,j).gt.0. ) then
                do ns = 1 , nsoil          !soil size segregatoin no
                  do nm = 1 , mode       !soil mode = 3
                    if ( (pcent(nm,nt).gt.eps) .and.                    &
                       & (sigma(nm,nt).ne.0.0) ) then
                      xk = pcent(nm,nt)/(sqrt(twopi)*log(sigma(nm,nt)))
                      xl = ((log(dp(ns))-log(mmd(nm,nt)*1.E-4))**2)     &
                         & /(2.0*(log(sigma(nm,nt)))**2)
                      xm = xk*exp(-xl)
                    else
                      xm = 0.0
                    end if
                    xn = rhop*(2.0/3.0)*(dp(ns)/2.0)
                    deldp = 0.0460517018598807
                                              ! dp(2)-dp(1) ss(nsoil)
                    ss(ns) = ss(ns) + (xm*deldp/xn)
                  end do
                  stotal = stotal + ss(ns)
                end do
                do ns = 1 , nsoil
                  if ( stotal.gt.0.0 ) srel(i,ns,nt) = ss(ns)/stotal
                                                !srel(iy,nsoil,nats)
                end do
              end if
            end do                    ! soil types
          end if
        end do
        do i = 2 , iym2
          itex = nint(dustsotex(i,j))
          if ( itex.ge.1 .and. itex.le.nats ) then
            do ns = 1 , nsoil
              srel2d(i,j,ns) = srel(i,ns,itex)
            end do
          end if
        end do
      end do
#endif

      end subroutine inidust

      end module mod_dust
