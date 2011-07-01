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

module mod_fvgcm

  use mod_dynparam
  use mod_memutil
  use m_realkinds
  use m_die
  use m_zeit
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil

  private

  integer , parameter :: nlev = 18 , nlat = 181 , nlon = 288

  real(sp) , dimension(nlev+1) :: ak , bk
  real(sp) , dimension(nlev) :: pplev , sigma1 , sigmar
  real(sp) , dimension(nlat) :: vlat
  real(sp) , dimension(nlon) :: vlon

  real(sp) , target , dimension(nlon,nlat,nlev*4+1) :: bb
  real(sp) , target , dimension(nlon,nlat,nlev*3) :: b2
  real(sp) , target , dimension(nlon,nlat,nlev*2) :: d2
  real(sp) , pointer , dimension(:,:,:) :: b3
  real(sp) , pointer , dimension(:,:,:) :: d3

  real(sp) , dimension(nlon,nlat) :: zs2
  real(sp) , dimension(nlon,nlat,nlev) :: pp3d , z1

  real(sp) , pointer , dimension(:,:) :: ps2
  real(sp) , pointer , dimension(:,:,:) :: q2 , t2 , u2 , v2
  real(sp) , pointer , dimension(:,:,:) :: tp , qp , hp
  real(sp) , pointer , dimension(:,:,:) :: up , vp
  real(sp) , pointer , dimension(:,:,:) :: t3 , q3 , h3
  real(sp) , pointer , dimension(:,:,:) :: u3 , v3

  public :: getfvgcm , headerfv

  contains

  subroutine getfvgcm(idate)
  implicit none
!
  type(rcm_time_and_date) , intent(in) :: idate
!
  character(3) , dimension(12) :: chmon
  character(20) :: finm , fips
  character(5) :: fn_a2 , fn_rf , pn_a2 , pn_rf
  integer :: i , i2 , ii , j , j2 , k , mrec , nrec , numx , numy
  integer(2) , dimension(288,181) :: itmp
  real(dp) :: offset , xscale
  real(sp) , dimension(288,181) :: temp
  logical :: there
  character(4) , dimension(30) :: yr_a2 , yr_rf
!
  data fn_rf/'FV_RF'/ , fn_a2/'FV_A2'/
  data pn_rf/'PS_RF'/ , pn_a2/'PS_A2'/
  data yr_rf/'1961' , '1962' , '1963' , '1964' , '1965' , '1966' , &
             '1967' , '1968' , '1969' , '1970' , '1971' , '1972' , &
             '1973' , '1974' , '1975' , '1976' , '1977' , '1978' , &
             '1979' , '1980' , '1981' , '1982' , '1983' , '1984' , &
             '1985' , '1986' , '1987' , '1988' , '1989' , '1990'/
  data yr_a2/'2071' , '2072' , '2073' , '2074' , '2075' , '2076' ,  &
             '2077' , '2078' , '2079' , '2080' , '2081' , '2082' ,  &
             '2083' , '2084' , '2085' , '2086' , '2087' , '2088' ,  &
             '2089' , '2090' , '2091' , '2092' , '2093' , '2094' ,  &
             '2095' , '2096' , '2097' , '2098' , '2099' , '2100'/
  data chmon/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , 'JUL' , &
             'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/
!
  call zeit_ci('getfvgcm')
  if ( idate == globidate1 ) then
    numx = nint((lon1-lon0)/1.25) + 1
    numy = nint(lat1-lat0) + 1
    inquire (file=trim(inpglob)//'/FVGCM/HT_SRF',exist=there)
    if ( .not.there ) then
      write (stderr,*) trim(inpglob)//'/FVGCM/HT_SRF is not present'
      call die('getfvgcm')
    end if
    open (61,file=trim(inpglob)//'/FVGCM/HT_SRF',form='unformatted', &
          recl=numx*numy*ibyte,access='direct')
    read (61,rec=1) ((temp(i,j),i=1,numx),j=1,numy)
    do j = nint(lat0) , nint(lat1)
      do i = nint(lon0/1.25) , nint(lon1/1.25)
        ii = i + 1
        if ( ii <= 0 ) ii = ii + 288
        if ( ii > 288 ) ii = ii - 288
        i2 = i - nint(lon0/1.25) + 1
        j2 = j - nint(lat0) + 1
        zs2(ii,j+91) = temp(i2,j2)/9.80616
      end do
    end do
    close (61)
  end if
 
  if ( idate%day /= 1 .or. idate%hour /= 0 ) then
    if ( ssttyp == 'FV_RF' ) then
      finm = 'RF/'//yr_rf(idate%year-1960)//'/'//fn_rf//yr_rf(idate%year-1960) &
             //chmon(idate%month)
      fips = 'RF/'//yr_rf(idate%year-1960)//'/'//pn_rf//yr_rf(idate%year-1960) &
             //chmon(idate%month)
    else if ( ssttyp == 'FV_A2' ) then
      finm = 'A2/'//yr_a2(idate%year-2070)//'/'//fn_a2//yr_a2(idate%year-2070) &
             //chmon(idate%month)
      fips = 'A2/'//yr_a2(idate%year-2070)//'/'//pn_a2//yr_a2(idate%year-2070) &
             //chmon(idate%month)
    else
      write (stderr,*) 'Unknown sstyp. Supported FV_RF and FV_A2'
      call die('getfvgcm')
    end if
  else if ( idate%month /= 1 ) then
    if ( ssttyp == 'FV_RF' ) then
      finm = 'RF/'//yr_rf(idate%year-1960)//'/'//fn_rf//yr_rf(idate%year-1960) &
             //chmon(idate%month-1)
      fips = 'RF/'//yr_rf(idate%year-1960)//'/'//pn_rf//yr_rf(idate%year-1960) &
             //chmon(idate%month-1)
    else if ( ssttyp == 'FV_A2' ) then
      finm = 'A2/'//yr_a2(idate%year-2070)//'/'//fn_a2//yr_a2(idate%year-2070) &
             //chmon(idate%month-1)
      fips = 'A2/'//yr_a2(idate%year-2070)//'/'//pn_a2//yr_a2(idate%year-2070) &
             //chmon(idate%month-1)
    else
      write (stderr,*) 'Unknown sstyp. Supported FV_RF and FV_A2'
      call die('getfvgcm')
    end if
  else if ( ssttyp == 'FV_RF' ) then
    if ( idate%year == 1961 ) then
      write (stderr,*) 'Fields on 00z01jan1961 is not saved'
      write (stderr,*) 'Please run from 00z02jan1961'
      call die('getfvgcm')
    end if
    finm = 'RF/'//yr_rf(idate%year-1961)//'/'//fn_rf//yr_rf(idate%year-1961)  &
           //chmon(12)
    fips = 'RF/'//yr_rf(idate%year-1961)//'/'//pn_rf//yr_rf(idate%year-1961)  &
           //chmon(12)
  else if ( ssttyp == 'FV_A2' ) then
    if ( idate%year == 2071 ) then
      write (stderr,*) 'Fields on 00z01jan2071 is not saved'
      write (stderr,*) 'Please run from 00z02jan2071'
      call die('getfvgcm')
    end if
    finm = 'A2/'//yr_a2(idate%year-2071)//'/'//fn_a2//yr_a2(idate%year-2071)  &
           //chmon(12)
    fips = 'A2/'//yr_a2(idate%year-2071)//'/'//pn_a2//yr_a2(idate%year-2071)  &
           //chmon(12)
  else
    write (stderr,*) 'Unknown sstyp. Supported FV_RF and FV_A2'
    call die('getfvgcm')
  end if
  numx = nint((lon1-lon0)/1.25) + 1
  numy = nint(lat1-lat0) + 1
  do k = 1 , nlev*4 + 1
    do j = 1 , nlat
      do i = 1 , nlon
        bb(i,j,k) = -9999.
      end do
    end do
  end do
  inquire (file=trim(inpglob)//'/FVGCM/'//finm,exist=there)
  if ( .not.there ) then
    write (stderr,*) trim(inpglob)//'/FVGCM/'//finm//' is not available'
    write (stderr,*) 'please copy FVGCM output under ', &
                trim(inpglob)//'/FVGCM/'
    call die('getfvgcm')
  end if
  open (63,file=trim(inpglob)//'/FVGCM/'//finm,form='unformatted',  &
        recl=(numx*numy*2+16)/4*ibyte,access='direct')
  open (62,file=trim(inpglob)//'/FVGCM/'//fips,form='unformatted',  &
        recl=numx*numy*ibyte,access='direct')
  if ( idate%day /= 1 .or. idate%hour /= 0 ) then
    nrec = ((idate%day-1)*4+idate%hour/6-1)*(nlev*4)
    mrec = (idate%day-1)*4 + idate%hour/6 - 1
  else if ( idate%month == 1 .or. idate%month == 2 .or. &
            idate%month == 4 .or. idate%month == 6 .or. &
            idate%month == 8 .or. idate%month == 9 .or. &
            idate%month == 11 ) then
    nrec = (31*4-1)*(nlev*4)
    mrec = 31*4 - 1
  else if ( idate%month == 5 .or. idate%month == 7 .or. &
            idate%month == 10 .or. idate%month == 12 ) then
    nrec = (30*4-1)*(nlev*4)
    mrec = 30*4 - 1
  else if ( mod(idate%year,4) == 0 .and. idate%year /= 2100 ) then
    nrec = (29*4-1)*(nlev*4)
    mrec = 29*4 - 1
  else
    nrec = (28*4-1)*(nlev*4)
    mrec = 28*4 - 1
  end if
  mrec = mrec + 1
  read (62,rec=mrec) ((temp(i,j),i=1,numx),j=1,numy)
  do j = nint(lat0) , nint(lat1)
    do i = nint(lon0/1.25) , nint(lon1/1.25)
      ii = i + 1
      if ( ii <= 0 ) ii = ii + 288
      if ( ii > 288 ) ii = ii - 288
      i2 = i - nint(lon0/1.25) + 1
      j2 = j - nint(lat0) + 1
      ps2(ii,j+91) = temp(i2,j2)*0.01
    end do
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) offset , xscale , ((itmp(i,j),i=1,numx),j=1,numy)
    do j = nint(lat0) , nint(lat1)
      do i = nint(lon0/1.25) , nint(lon1/1.25)
        ii = i + 1
        if ( ii <= 0 ) ii = ii + 288
        if ( ii > 288 ) ii = ii - 288
        i2 = i - nint(lon0/1.25) + 1
        j2 = j - nint(lat0) + 1
        u2(ii,j+91,k) = real(dble(itmp(i2,j2))*xscale + offset)
      end do
    end do
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) offset , xscale , ((itmp(i,j),i=1,numx),j=1,numy)
    do j = nint(lat0) , nint(lat1)
      do i = nint(lon0/1.25) , nint(lon1/1.25)
        ii = i + 1
        if ( ii <= 0 ) ii = ii + 288
        if ( ii > 288 ) ii = ii - 288
        i2 = i - nint(lon0/1.25) + 1
        j2 = j - nint(lat0) + 1
        v2(ii,j+91,k) = real(dble(itmp(i2,j2))*xscale + offset)
      end do
    end do
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) offset , xscale , ((itmp(i,j),i=1,numx),j=1,numy)
    do j = nint(lat0) , nint(lat1)
      do i = nint(lon0/1.25) , nint(lon1/1.25)
        ii = i + 1
        if ( ii <= 0 ) ii = ii + 288
        if ( ii > 288 ) ii = ii - 288
        i2 = i - nint(lon0/1.25) + 1
        j2 = j - nint(lat0) + 1
        t2(ii,j+91,k) = real(dble(itmp(i2,j2))*xscale + offset)
      end do
    end do
  end do
  do k = 1 , nlev
    nrec = nrec + 1
    read (63,rec=nrec) offset , xscale , ((itmp(i,j),i=1,numx),j=1,numy)
    do j = nint(lat0) , nint(lat1)
      do i = nint(lon0/1.25) , nint(lon1/1.25)
        ii = i + 1
        if ( ii <= 0 ) ii = ii + 288
        if ( ii > 288 ) ii = ii - 288
        i2 = i - nint(lon0/1.25) + 1
        j2 = j - nint(lat0) + 1
        q2(ii,j+91,k) = real(dble(itmp(i2,j2))*xscale + offset)
      end do
    end do
  end do
  close (63)
  close (62)
  write (stdout,*) 'READ IN fields at DATE:' , idate
  do k = 1 , nlev
    do j = 1 , nlat
      do i = 1 , nlon
        if ( ps2(i,j) > -9995. ) then
          pp3d(i,j,k) = ps2(i,j)*0.5*(bk(k)+bk(k+1))    &
                        + 0.5*(ak(k)+ak(k+1))
        else
          pp3d(i,j,k) = -9999.0
        end if
      end do
    end do
  end do
!
!     to calculate Heights on sigma surfaces.
  call htsig(t2,z1,pp3d,ps2,zs2,nlon,nlat,nlev)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
  call height(hp,z1,t2,ps2,pp3d,zs2,nlon,nlat,nlev,pplev,nlev)
!     2. For Zonal and Meridional Winds
  call intlin(up,u2,ps2,pp3d,nlon,nlat,nlev,pplev,nlev)
  call intlin(vp,v2,ps2,pp3d,nlon,nlat,nlev,pplev,nlev)
!     3. For Temperatures
  call intlog(tp,t2,ps2,pp3d,nlon,nlat,nlev,pplev,nlev)
!     4. For Moisture qva & qca
  call humid1fv(t2,q2,pp3d,nlon,nlat,nlev)
  call intlin(qp,q2,ps2,pp3d,nlon,nlat,nlev,pplev,nlev)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
  call bilinx2(b3,b2,xlon,xlat,vlon,vlat,nlon,nlat,jx,iy,nlev*3)
  call bilinx2(d3,d2,dlon,dlat,vlon,vlat,nlon,nlat,jx,iy,nlev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
  call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,nlev,plon,plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
  call top2btm(t3,jx,iy,nlev)
  call top2btm(q3,jx,iy,nlev)
  call top2btm(h3,jx,iy,nlev)
  call top2btm(u3,jx,iy,nlev)
  call top2btm(v3,jx,iy,nlev)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,nlev)
 
  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,nlev)
 
  call readsst(ts4,idate)

!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,nlev)
!
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev)
 
  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,nlev)
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4     DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  call zeit_co('getfvgcm')
  end subroutine getfvgcm
!
!
  subroutine headerfv

  implicit none
!
  integer :: i , j , k , kr
!
  do j = 1 , nlat
    vlat(j) = float(j-1) - 90.0
  end do
!
  do i = 1 , nlon
    vlon(i) = float(i-1)*1.25
  end do
 
  pplev(1) = 30.
  pplev(2) = 50.
  pplev(3) = 70.
  pplev(4) = 100.
  pplev(5) = 150.
  pplev(6) = 200.
  pplev(7) = 250.
  pplev(8) = 300.
  pplev(9) = 350.
  pplev(10) = 420.
  pplev(11) = 500.
  pplev(12) = 600.
  pplev(13) = 700.
  pplev(14) = 780.
  pplev(15) = 850.
  pplev(16) = 920.
  pplev(17) = 960.
  pplev(18) = 1000.
 
  do k = 1 , nlev
    sigmar(k) = pplev(k)*0.001
  end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , nlev
    kr = nlev - k + 1
    sigma1(k) = sigmar(kr)
  end do
!
  ak(1) = 2.9170
  ak(2) = 7.9292
  ak(3) = 21.5539
  ak(4) = 49.1834
  ak(5) = 83.1425
  ak(6) = 79.9308
  ak(7) = 75.7738
  ak(8) = 70.5752
  ak(9) = 64.2963
  ak(10) = 56.9838
  ak(11) = 48.7913
  ak(12) = 39.9895
  ak(13) = 30.9631
  ak(14) = 22.1902
  ak(15) = 14.2039
  ak(16) = 7.5413
  ak(17) = 2.6838
  ak(18) = 0.0
  ak(19) = 0.0
 
  bk(1) = 0.0
  bk(2) = 0.0
  bk(3) = 0.0
  bk(4) = 0.0
  bk(5) = 0.0
  bk(6) = 0.0380541
  bk(7) = 0.0873088
  bk(8) = 0.1489307
  bk(9) = 0.2232996
  bk(10) = 0.3099406
  bk(11) = 0.4070096
  bk(12) = 0.5112977
  bk(13) = 0.6182465
  bk(14) = 0.7221927
  bk(15) = 0.8168173
  bk(16) = 0.8957590
  bk(17) = 0.9533137
  bk(18) = 0.9851222
  bk(19) = 1.0
 
  call getmem3d(b3,1,jx,1,iy,1,nlev*3,'mod_fvgcm:b3')
  call getmem3d(d3,1,jx,1,iy,1,nlev*2,'mod_fvgcm:d3')

!     Set up pointers

  ps2 => bb(:,:,1)
  t2 => bb(:,:,2:nlev+1)
  q2 => bb(:,:,nlev+2:2*nlev+1)
  u2 => bb(:,:,2*nlev+2:3*nlev+1)
  v2 => bb(:,:,3*nlev+2:4*nlev+1)
  tp => b2(:,:,1:nlev)
  qp => b2(:,:,nlev+1:2*nlev)
  hp => b2(:,:,2*nlev+1:3*nlev)
  up => d2(:,:,1:nlev)
  vp => d2(:,:,nlev+1:2*nlev)
  t3 => b3(:,:,1:nlev)
  q3 => b3(:,:,nlev+1:2*nlev)
  h3 => b3(:,:,2*nlev+1:3*nlev)
  u3 => d3(:,:,1:nlev)
  v3 => d3(:,:,nlev+1:2*nlev)

  end subroutine headerfv

end module mod_fvgcm
