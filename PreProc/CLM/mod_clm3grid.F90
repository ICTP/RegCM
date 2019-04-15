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

module mod_clm3grid

  use mod_intkinds
  use mod_realkinds
  use mod_stdio

  contains

  subroutine clm3grid1(nlon,nlat,nlev,ntim, &
                       glon1,glon2,glat1,glat2,dlat,dlon, &
                       xlonmin,xlonmax,xlatmin,xlatmax,istart,icount)
  implicit none
!
  real(rk4) :: glat1 , glat2 , glon1 , glon2 , dlat , dlon
  real(rk4) :: xlatmax , xlatmin , xlonmax , xlonmin
  integer(ik4) :: nlat , nlev , nlon , ntim
  integer(ik4) , dimension(4) :: icount , istart
  intent (in) glat1 , glat2 , nlat , nlev , nlon , ntim , dlat , dlon
  intent (out) icount , istart
  intent (inout) xlatmax , xlatmin , xlonmax , xlonmin
!
  integer(ik4) :: ilatmax , ilatmin , ilonmax , ilonmin

  xlonmin = max(xlonmin-dlon,glon1)
  xlonmax = min(xlonmax+dlon,glon2)
  xlatmin = max(xlatmin-dlat,glat1)
  xlatmax = min(xlatmax+dlat,glat2)
  ilonmin = max(min(nint((glon2+xlonmin)/dlon)-1,nlon),1)
  ilonmax = max(min(nint((glon2+xlonmax)/dlon)+1,nlon),1)
!  ilatmin = max(min(nint((glat2+xlatmin+glat2+glat1)/dlat)-1,nlat),1)
!  ilatmax = max(min(nint((glat2+xlatmax+glat2+glat1)/dlat)+1,nlat),1)
  ilatmin = max(min(nint((glat2+xlatmin-nint(glat2)-nint(glat1))/dlat)-1,nlat),1)
  ilatmax = max(min(nint((glat2+xlatmax-nint(glat2)-nint(glat1))/dlat)+1,nlat),1)

  istart(1) = ilonmin
  icount(1) = ilonmax - ilonmin + 1
  istart(2) = ilatmin
  icount(2) = ilatmax - ilatmin + 1
  istart(3) = 1
  icount(3) = nlev
  istart(4) = 1
  icount(4) = ntim

  end subroutine clm3grid1

  subroutine clm3grid2(nlon,nlat,glon,glat,istart,icount,zlon,      &
                       zlat,zlev)

  implicit none
!
  integer(ik4) :: nlat , nlon
  real(rk4) , dimension(nlat) :: glat
  real(rk4) , dimension(nlon) :: glon
  integer(ik4) , dimension(4) :: icount , istart
  real(rk4) , dimension(icount(2)) :: zlat
  real(rk4) , dimension(icount(3)) :: zlev
  real(rk4) , dimension(icount(1)) :: zlon
  intent (in) glat , glon , icount , istart , nlat , nlon
  intent (out) zlat , zlev , zlon
!
  integer(ik4) :: i , j , k
!
  do i = 1 , icount(1)
    zlon(i) = glon(i+istart(1)-1)
  end do
  do j = 1 , icount(2)
    zlat(j) = glat(j+istart(2)-1)
  end do
  do k = 1 , icount(3)
    zlev(k) = real(icount(3) - k + 1)
  end do

  end subroutine clm3grid2

  subroutine bilinx4d(mti,lmsk,loni,lati,nloni,nlati,mto,lono,lato,jx,iy,&
                      nz,nt,xming,vmisdat)
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
  integer(ik4) :: iy , jx , nlati , nloni , nt , nz
  real(rk4) :: vmisdat , xming
  real(rk4) , dimension(nloni,nlati,nz,nt) :: mti
  real(rk4) , dimension(nloni,nlati) :: lmsk
  real(rk4) , dimension(nlati) :: lati
  real(rk4) , dimension(jx,iy) :: lato , lono
  real(rk4) , dimension(nloni) :: loni
  real(rk4) , dimension(jx,iy,nz,nt) :: mto
  intent (in) mti , iy , jx , lati , lato , loni , lono , nlati ,   &
              nloni , nt , nz , vmisdat , xming
  intent (out) mto
!
  integer(ik4) :: i , ip , ipp1 , j , jq , jqp1 , k , l
  real(rk4) :: lon360 , p , q , temp1 , temp2 , xind , yind
  logical :: gt1 , gt2
!
  do i = 1 , iy
    do j = 1 , jx

      yind = (((lato(j,i)-lati(1))/(lati(nlati)-lati(1))) * &
              float(nlati-1)) + 1.0
      jq = int(yind)
      jqp1 = min0(jq+1,nlati)
      q = yind - real(jq)

      lon360 = lono(j,i)
      xind = (((lon360-loni(1))/(loni(nloni)-loni(1))) * &
              float(nloni-1)) + 1.0
      ip = int(xind)
      ipp1 = min0(ip+1,nloni)
      p = xind - real(ip)

      do l = 1 , nt
        do k = 1 , nz
          gt1 = .false.
          gt2 = .false.
          temp1 = vmisdat
          temp2 = vmisdat
          if ( (lmsk(ip,jq)   < 0.5 .or.  mti(ip,jq,k,l)   <= xming)  .and. &
               (lmsk(ipp1,jq) > 0.5 .and. mti(ipp1,jq,k,l) >  xming) )  then
            temp1 = mti(ipp1,jq,k,l)
            gt1 = .true.
          else if ( (lmsk(ip,jq) > 0.5 .and. mti(ip,jq,k,l)   >  xming) .and. &
                  (lmsk(ipp1,jq) < 0.5 .or.  mti(ipp1,jq,k,l) <= xming) )  then
            temp1 = mti(ip,jq,k,l)
            gt1 = .true.
          else if ( (lmsk(ip,jq)   > 0.5 .and. mti(ip,jq,k,l)   >  xming) .and.&
                    (lmsk(ipp1,jq) > 0.5 .and. mti(ipp1,jq,k,l) >  xming) ) then
            temp1 = (1.0-p)*mti(ip,jq,k,l) + p*mti(ipp1,jq,k,l)
            gt1 = .true.
          end if
          if ( (lmsk(ipp1,jqp1) < 0.5 .or.  mti(ipp1,jqp1,k,l) <= xming) .and. &
               (lmsk(ip,  jqp1) > 0.5 .and. mti(ip,jqp1,k,l)   >  xming) )  then
            temp2 = mti(ip,jqp1,k,l)
            gt2 = .true.
          else if ((lmsk(ipp1,jqp1) > 0.5 .and. &
                    mti(ipp1,jqp1,k,l) > xming) .and. &
                   (lmsk(ip,jqp1)   < 0.5 .or.  &
                   mti(ip,jqp1,k,l)  <= xming) )  then
            temp2 = mti(ipp1,jqp1,k,l)
            gt2 = .true.
          else if ( (lmsk(ipp1,jqp1) > 0.5 .and. &
                     mti(ipp1,jqp1,k,l) > xming) .and. &
                    (lmsk(ip,jqp1)   > 0.5 .and. &
                     mti(ip,jqp1,k,l) > xming) ) then
            temp2 = (1.0-p)*mti(ip,jqp1,k,l) + p*mti(ipp1,jqp1,k,l)
            gt2 = .true.
          end if
          if ( .not. gt1 .and. .not. gt2 ) then
            mto(j,i,k,l) = vmisdat
          else if ( .not. gt1 ) then
            mto(j,i,k,l) = temp2
          else if ( .not. gt2 ) then
            mto(j,i,k,l) = temp1
          else
            mto(j,i,k,l) = (1.0-q)*temp1 + q*temp2
          end if
        end do
      end do
    end do
  end do

  end subroutine bilinx4d

  subroutine maskme(landmask,vals,vmisdat,nlon,nlat,nlev,ntim)
  implicit none
!
  integer(ik4) :: nlat , nlev , nlon , ntim
  real(rk4) :: vmisdat
  real(rk4) , dimension(nlon,nlat) :: landmask
  real(rk4) , dimension(nlon,nlat,nlev,ntim) :: vals
  intent (in) landmask , nlat , nlev , nlon , ntim , vmisdat
  intent (inout) vals
!
  integer(ik4) :: i , j , k , l
!
  do l = 1 , ntim
    do k = 1 , nlev
      do j = 1 , nlat
        do i = 1 , nlon
          if ( landmask(i,j)>0.5 ) then
            vals(i,j,k,l) = vals(i,j,k,l)
          else
            vals(i,j,k,l) = vmisdat
          end if
        end do
      end do
    end do
  end do

  end subroutine maskme

end module mod_clm3grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
