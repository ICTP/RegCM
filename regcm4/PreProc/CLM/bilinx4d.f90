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

      subroutine bilinx4d(mti,loni,lati,nloni,nlati,mto,lono,lato,iy,jx,&
                        & nz,nt,xming,vmisdat)
!
!  PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A BIGGER
!  RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF GRID2.
!  A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON GRID4.THE
!  GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE TRAPPED POINT.
!  THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES IN BOTH X AND Y
!  DIRECTION OF THE TRAPPED GRID POINT AND USES THE INFORMATION
!  AS WEIGHTING FACTORS IN THE INTERPOLATION.
!  THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!  INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
!  THE CROSS POINTS IN THE RCM MODEL.
!
!  IN(NLONI,NLATI,NZ)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!  OUT(NLATO,NLONO,NZ) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL GRID.
!  LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!  LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!  P.........EAST-WEST WEIGHTING FACTOR.
!  Q.........NORTH-SOUTH WEIGHTING FACTOR.
!  IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!  IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID POINT.
 
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nlati , nloni , nt , nz
      real(4) :: vmisdat , xming
      real(4) , dimension(nloni,nlati,nz,nt) :: mti
      real(4) , dimension(nlati) :: lati
      real(4) , dimension(iy,jx) :: lato , lono
      real(4) , dimension(nloni) :: loni
      real(4) , dimension(iy,jx,nz,nt) :: mto
      intent (in) mti , iy , jx , lati , lato , loni , lono , nlati ,   &
                & nloni , nt , nz , vmisdat , xming
      intent (out) mto
!
! Local variables
!
      integer :: i , ip , ipp1 , j , jq , jqp1 , k , l
      real(4) :: lon360 , p , q , temp1 , temp2 , xind , yind
!
      do j = 1 , jx
        do i = 1 , iy
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon360 = lono(i,j)
          xind = (((lon360-loni(1))/(loni(nloni)-loni(1)))              &
               & *float(nloni-1)) + 1.
          ip = int(xind)
          ipp1 = min0(ip+1,nloni)
          p = xind - ip
 
          do l = 1 , nt
            do k = 1 , nz
 
              if ( mti(ip,jq,k,l)<=xming .and. mti(ipp1,jq,k,l)>xming ) &
                 & then
                temp1 = mti(ipp1,jq,k,l)
              else if ( mti(ip,jq,k,l)<=xming .and. mti(ipp1,jq,k,l)    &
                      & <=xming ) then
                temp1 = vmisdat
              else if ( mti(ip,jq,k,l)>xming .and. mti(ipp1,jq,k,l)     &
                      & <=xming ) then
                temp1 = mti(ip,jq,k,l)
              else
                temp1 = (1.0-p)*mti(ip,jq,k,l) + p*mti(ipp1,jq,k,l)
              end if
 
              if ( mti(ipp1,jqp1,k,l)<=xming .and. mti(ip,jqp1,k,l)     &
                 & >xming ) then
                temp2 = mti(ip,jqp1,k,l)
              else if ( mti(ipp1,jqp1,k,l)<=xming .and. mti(ip,jqp1,k,l)&
                      & <=xming ) then
                temp2 = vmisdat
              else if ( mti(ipp1,jqp1,k,l)>xming .and. mti(ip,jqp1,k,l) &
                      & <=xming ) then
                temp2 = mti(ipp1,jqp1,k,l)
              else
                temp2 = p*mti(ipp1,jqp1,k,l) + (1.0-p)*mti(ip,jqp1,k,l)
              end if
 
              if ( temp1<=xming .and. temp2>xming ) then
                mto(i,j,k,l) = temp2
              else if ( temp1<=xming .and. temp2<=xming ) then
                mto(i,j,k,l) = vmisdat
              else if ( temp1>xming .and. temp2<=xming ) then
                mto(i,j,k,l) = temp1
              else
                mto(i,j,k,l) = (1.-q)*temp1 + q*temp2
              end if
 
            end do
          end do
 
        end do
 
      end do
 
      end subroutine bilinx4d
