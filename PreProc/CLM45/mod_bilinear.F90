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

module mod_bilinear

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_message

  interface bilinear
    module procedure bilinear2d
    module procedure bilinear2d_3d_in
    module procedure bilinear2d_4d_in
  end interface bilinear

  contains
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
  !
  subroutine bilinear2d(mti,lmsk,loni,lati,mto,lono,lato,xming,vmisdat)
    implicit none
    real(rk4) , intent(in) , dimension(:,:) :: mti
    real(rk4) , intent(in) , dimension(:,:) :: lmsk
    real(rk4) , intent(in) , dimension(:) :: lati , loni
    real(rk4) , intent(out) , dimension(:,:) :: mto
    real(rk4) , intent(in) , dimension(:,:) :: lato , lono
    real(rk4) , intent(in) :: vmisdat , xming

    integer(ik4) :: ni , nj , nloni , nlati
    integer(ik4) :: i , ip , ipp1 , j , jq , jqp1
    real(rk4) :: lon360 , p , q , temp1 , temp2 , xind , yind
    logical :: gt1 , gt2

    nj = size(mto,1)
    ni = size(mto,2)
    nloni = size(loni)
    nlati = size(lati)

    ! Check dimensions
    if ( nj /= size(lato,1) .or. nj /= size(lono,1) .or. &
         ni /= size(lato,2) .or. ni /= size(lono,2) .or. &
         nloni /= size(mti,1) .or. nloni /= size(lmsk,1) .or. &
         nlati /= size(mti,2) .or. nlati /= size(lmsk,2) ) then
      write(stderr,*) 'Dimension error in bilinear2d'
      call die(__FILE__,'Now stopping',__LINE__)
    end if

    do i = 1 , ni
      do j = 1 , nj
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
        gt1 = .false.
        gt2 = .false.
        temp1 = vmisdat
        temp2 = vmisdat
        if ( (lmsk(ip,jq)   < 0.5 .or.  mti(ip,jq)   <= xming)  .and. &
             (lmsk(ipp1,jq) > 0.5 .and. mti(ipp1,jq) >  xming) )  then
          temp1 = mti(ipp1,jq)
          gt1 = .true.
        else if ( (lmsk(ip,jq)   > 0.5 .and. mti(ip,jq)   >  xming) .and. &
                  (lmsk(ipp1,jq) < 0.5 .or.  mti(ipp1,jq) <= xming) )  then
          temp1 = mti(ip,jq)
          gt1 = .true.
        else if ( (lmsk(ip,jq)   > 0.5 .and. mti(ip,jq)   >  xming) .and. &
                  (lmsk(ipp1,jq) > 0.5 .and. mti(ipp1,jq) >  xming) ) then
          temp1 = (1.0-p)*mti(ip,jq) + p*mti(ipp1,jq)
          gt1 = .true.
        end if
        if ( (lmsk(ipp1,jqp1) < 0.5 .or.  mti(ipp1,jqp1) <= xming) .and. &
             (lmsk(ip,  jqp1) > 0.5 .and. mti(ip,jqp1)   >  xming) )  then
          temp2 = mti(ip,jqp1)
          gt2 = .true.
        else if ( (lmsk(ipp1,jqp1) > 0.5 .and. mti(ipp1,jqp1) > xming) .and. &
                  (lmsk(ip,jqp1)   < 0.5 .or.  mti(ip,jqp1)  <= xming) )  then
          temp2 = mti(ipp1,jqp1)
          gt2 = .true.
        else if ( (lmsk(ipp1,jqp1) > 0.5 .and. mti(ipp1,jqp1) > xming) .and. &
                  (lmsk(ip,jqp1)   > 0.5 .and. mti(ip,jqp1)   > xming) ) then
          temp2 = (1.0-p)*mti(ip,jqp1) + p*mti(ipp1,jqp1)
          gt2 = .true.
        end if
        if ( .not. gt1 .and. .not. gt2 ) then
          mto(j,i) = vmisdat
        else if ( .not. gt1 ) then
          mto(j,i) = temp2
        else if ( .not. gt2 ) then
          mto(j,i) = temp1
        else
          mto(j,i) = (1.0-q)*temp1 + q*temp2
        end if
      end do
    end do
  end subroutine bilinear2d

  subroutine bilinear2d_3d_in(mti,lmsk,loni,lati,mto,lono,lato,xming,vmisdat)
    implicit none
    real(rk4) , intent(in) , dimension(:,:,:) :: mti
    real(rk4) , intent(in) , dimension(:,:) :: lmsk
    real(rk4) , intent(in) , dimension(:) :: lati , loni
    real(rk4) , intent(out) , dimension(:,:,:) :: mto
    real(rk4) , intent(in) , dimension(:,:) :: lato , lono
    real(rk4) , intent(in) :: vmisdat , xming
    integer(ik4) :: k , nk
    nk = size(mti,3)
    if ( size(mto,3) /= nk ) then
      write(stderr,*) 'Dimension error in bilinear2d_3d_in'
      call die(__FILE__,'Now stopping',__LINE__)
    end if
    do k = 1 , nk
      call bilinear2d(mti(:,:,k),lmsk,loni,lati,mto(:,:,k), &
                      lono,lato,xming,vmisdat)
    end do
  end subroutine bilinear2d_3d_in

  subroutine bilinear2d_4d_in(mti,lmsk,loni,lati,mto,lono,lato,xming,vmisdat)
    implicit none
    real(rk4) , intent(in) , dimension(:,:,:,:) :: mti
    real(rk4) , intent(in) , dimension(:,:) :: lmsk
    real(rk4) , intent(in) , dimension(:) :: lati , loni
    real(rk4) , intent(out) , dimension(:,:,:,:) :: mto
    real(rk4) , intent(in) , dimension(:,:) :: lato , lono
    real(rk4) , intent(in) :: vmisdat , xming
    integer(ik4) :: k , l , nk , nl
    nk = size(mti,3)
    nl = size(mti,4)
    if ( size(mto,3) /= nk .or. size(mto,4) /= nl ) then
      write(stderr,*) 'Dimension error in bilinear2d_4d_in'
      call die(__FILE__,'Now stopping',__LINE__)
    end if
    do l = 1 , nl
      do k = 1 , nk
        call bilinear2d(mti(:,:,k,l),lmsk,loni,lati,mto(:,:,k,l), &
                        lono,lato,xming,vmisdat)
      end do
    end do
  end subroutine bilinear2d_4d_in

end module mod_bilinear
