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

    use mod_intkinds
    use mod_realkinds
    use mod_dynparam
    use mod_memutil
    use mod_message
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

    integer(ik4) , parameter :: nlev = 18
    integer(ik4) :: numx , numy

    real(rk8) , dimension(nlev+1) :: ak , bk
    real(rk8) , dimension(nlev) :: pplev , sigma1 , sigmar

    real(rk8) , pointer , dimension(:) :: vlat
    real(rk8) , pointer , dimension(:) :: vlon

    real(rk8) , pointer , dimension(:,:,:) :: bb
    real(rk8) , pointer , dimension(:,:,:) :: b2
    real(rk8) , pointer , dimension(:,:,:) :: d2
    real(rk8) , pointer , dimension(:,:,:) :: b3
    real(rk8) , pointer , dimension(:,:,:) :: d3
    real(rk4) , pointer , dimension(:,:) :: temp
    integer(2) , pointer , dimension(:,:) :: itmp

    real(rk8) , pointer , dimension(:,:) :: zs2
    real(rk8) , pointer , dimension(:,:,:) :: pp3d , z1

    real(rk8) , pointer , dimension(:,:) :: ps2
    real(rk8) , pointer , dimension(:,:,:) :: q2 , t2 , u2 , v2
    real(rk8) , pointer , dimension(:,:,:) :: tp , qp , hp
    real(rk8) , pointer , dimension(:,:,:) :: up , vp
    real(rk8) , pointer , dimension(:,:,:) :: t3 , q3 , h3
    real(rk8) , pointer , dimension(:,:,:) :: u3 , v3

    public :: getfvgcm , headerfv

    contains

    subroutine getfvgcm(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      character(len=3) , dimension(12) :: chmon
      character(len=20) :: finm , fips
      character(len=5) :: fn_a2 , fn_rf , pn_a2 , pn_rf
      integer(ik4) :: i , j , k , mrec , nrec
      real(rk8) :: offset , xscale
      logical :: there
      character(len=4) , dimension(30) :: yr_a2 , yr_rf
      integer(ik4) :: year , month , day , hour
      integer :: ilenrec

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

      call split_idate(idate,year,month,day,hour)

    if ( idate == globidate1 ) then
      inquire (file=trim(inpglob)//'/FVGCM/HT_SRF',exist=there)
      if ( .not.there ) then
        write (stderr,*) trim(inpglob)//'/FVGCM/HT_SRF is not present'
        call die('getfvgcm')
      end if
      inquire(iolength=ilenrec) temp
      open (61,file=trim(inpglob)//'/FVGCM/HT_SRF',form='unformatted', &
            recl=ilenrec,access='direct',action='read',status='old')
      read (61,rec=1) temp
      zs2 = temp*regrav
      close (61)
    end if

    if ( day /= 1 .or. hour /= 0 ) then
      if ( year < 2001 ) then
        finm = 'RF/'//yr_rf(year-1960)//'/'//fn_rf//yr_rf(year-1960) &
               //chmon(month)
        fips = 'RF/'//yr_rf(year-1960)//'/'//pn_rf//yr_rf(year-1960) &
               //chmon(month)
      else
        finm = dattyp(4:5)//'/'//yr_a2(year-2070)//'/'//&
                fn_a2//yr_a2(year-2070)//chmon(month)
        fips = dattyp(4:5)//'/'//yr_a2(year-2070)//'/'//&
                pn_a2//yr_a2(year-2070)//chmon(month)
      end if
    else if ( month /= 1 ) then
      if ( year < 2001 ) then
        finm = 'RF/'//yr_rf(year-1960)//'/'//fn_rf//yr_rf(year-1960) &
               //chmon(month-1)
        fips = 'RF/'//yr_rf(year-1960)//'/'//pn_rf//yr_rf(year-1960) &
               //chmon(month-1)
      else
        finm = dattyp(4:5)//'/'//yr_a2(year-2070)//'/'// &
                fn_a2//yr_a2(year-2070)//chmon(month-1)
        fips = dattyp(4:5)//'/'//yr_a2(year-2070)//'/'// &
                pn_a2//yr_a2(year-2070)//chmon(month-1)
      end if
    else
      if ( year == 1961 ) then
        write (stderr,*) 'Fields on 00z01jan1961 is not saved'
        write (stderr,*) 'Please run from 00z02jan1961'
        call die('getfvgcm')
      else if ( year == 2071 ) then
        write (stderr,*) 'Fields on 00z01jan2071 is not saved'
        write (stderr,*) 'Please run from 00z02jan2071'
        call die('getfvgcm')
      else if ( year < 2001 ) then
        finm = 'RF/'//yr_rf(year-1961)//'/'//fn_rf//yr_rf(year-1961)  &
               //chmon(12)
        fips = 'RF/'//yr_rf(year-1961)//'/'//pn_rf//yr_rf(year-1961)  &
               //chmon(12)
      else
        finm = dattyp(4:5)//'/'//yr_a2(year-2071)//'/'// &
                fn_a2//yr_a2(year-2071)//chmon(12)
        fips = dattyp(4:5)//'/'//yr_a2(year-2071)//'/'// &
                pn_a2//yr_a2(year-2071)//chmon(12)
      end if
    end if
    do k = 1 , nlev*4 + 1
      do j = 1 , numy
        do i = 1 , numx
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
    inquire(iolength=ilenrec)  offset , xscale , itmp
    open (63,file=trim(inpglob)//'/FVGCM/'//finm,form='unformatted',  &
          recl=ilenrec,access='direct',action='read',status='old')
    inquire(iolength=ilenrec) temp
    open (62,file=trim(inpglob)//'/FVGCM/'//fips,form='unformatted',  &
          recl=ilenrec,access='direct',action='read',status='old')
    if ( day /= 1 .or. hour /= 0 ) then
      nrec = ((day-1)*4+hour/6-1)*(nlev*4)
      mrec = (day-1)*4 + hour/6 - 1
    else if ( month == 1 .or. month == 2 .or. &
              month == 4 .or. month == 6 .or. &
              month == 8 .or. month == 9 .or. &
              month == 11 ) then
      nrec = (31*4-1)*(nlev*4)
      mrec = 31*4 - 1
    else if ( month == 5 .or. month == 7 .or. &
              month == 10 .or. month == 12 ) then
      nrec = (30*4-1)*(nlev*4)
      mrec = 30*4 - 1
    else if ( mod(year,4) == 0 .and. year /= 2100 ) then
      nrec = (29*4-1)*(nlev*4)
      mrec = 29*4 - 1
    else
      nrec = (28*4-1)*(nlev*4)
      mrec = 28*4 - 1
    end if
    mrec = mrec + 1
    read (62,rec=mrec) temp
    ps2 = temp*0.01
    do k = 1 , nlev
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      u2(:,:,k) = itmp*xscale + offset
    end do
    do k = 1 , nlev
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      v2(:,:,k) = itmp*xscale + offset
    end do
    do k = 1 , nlev
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      t2(:,:,k) = itmp*xscale + offset
    end do
    do k = 1 , nlev
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      q2(:,:,k) = itmp*xscale + offset
    end do
    close (63)
    close (62)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    do k = 1 , nlev
      do j = 1 , numy
        do i = 1 , numx
          if ( ps2(i,j) > -9995. ) then
            pp3d(i,j,k) = ps2(i,j)*0.5*(bk(k)+bk(k+1))    &
                          + 0.5*(ak(k)+ak(k+1))
          else
            pp3d(i,j,k) = -9999.0
          end if
        end do
      end do
    end do

    call htsig(t2,z1,pp3d,ps2,zs2,numx,numy,nlev)

    call height(hp,z1,t2,ps2,pp3d,zs2,numx,numy,nlev,pplev,nlev)

    call intlin(up,u2,ps2,pp3d,numx,numy,nlev,pplev,nlev)
    call intlin(vp,v2,ps2,pp3d,numx,numy,nlev,pplev,nlev)

    call intlog(tp,t2,ps2,pp3d,numx,numy,nlev,pplev,nlev)

    call humid1fv(t2,q2,pp3d,numx,numy,nlev)
    call intlin(qp,q2,ps2,pp3d,numx,numy,nlev,pplev,nlev)

    call bilinx2(b3,b2,xlon,xlat,vlon,vlat,numx,numy,jx,iy,nlev*3)
    call bilinx2(d3,d2,dlon,dlat,vlon,vlat,numx,numy,jx,iy,nlev*2)

    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,nlev,plon,plat,iproj)

    call top2btm(t3,jx,iy,nlev)
    call top2btm(q3,jx,iy,nlev)
    call top2btm(h3,jx,iy,nlev)
    call top2btm(u3,jx,iy,nlev)
    call top2btm(v3,jx,iy,nlev)

    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,nlev)

    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    if(i_band == 1) then
      call p1p2_band(b3pd,ps4,jx,iy)
    else
      call p1p2(b3pd,ps4,jx,iy)
    endif

    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,nlev)

    call readsst(ts4,idate)

    call intv1(u4,u3,b3pd,sigmah,sigmar,ptop,jx,iy,kz,nlev)
    call intv1(v4,v3,b3pd,sigmah,sigmar,ptop,jx,iy,kz,nlev)

    call intv2(t4,t3,ps4,sigmah,sigmar,ptop,jx,iy,kz,nlev)

    call intv1(q4,q3,ps4,sigmah,sigmar,ptop,jx,iy,kz,nlev)
    call humid2(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
    call hydrost(h4,t4,topogm,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine getfvgcm

  subroutine headerfv
    implicit none
    integer(ik4) :: i , j , k , kr
    numx = nint((lon1-lon0)/1.25) + 1
    numy = nint(lat1-lat0) + 1

    write (stdout,*) 'Reading a ',numx,'x',numy,' point grid'
    call getmem1d(vlat,1,numy,'fvgc:vlat')
    call getmem1d(vlon,1,numx,'fvgc:vlon')

    do j = 1 , numy
      vlat(j) = lat0 + dble(j-1)*1.0D0
    end do
    do i = 1 , numx
      vlon(i) = lon0 + dble(i-1)*1.25D0
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

    do k = 1 , nlev
      kr = nlev - k + 1
      sigma1(k) = sigmar(kr)
    end do

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

    call getmem3d(bb,1,numx,1,numy,1,nlev*4+1,'fvgc:bb')
    call getmem3d(b2,1,numx,1,numy,1,nlev*3,'fvgc:b2')
    call getmem3d(d2,1,numx,1,numy,1,nlev*2,'fvgc:d2')
    call getmem3d(pp3d,1,numx,1,numy,1,nlev,'fvgc:pp3d')
    call getmem3d(z1,1,numx,1,numy,1,nlev,'fvgc:z1')
    call getmem2d(zs2,1,numx,1,numy,'fvgc:zs2')
    call getmem2d(temp,1,numx,1,numy,'fvgc:temp')
    call getmem2d(itmp,1,numx,1,numy,'fvgc:itmp')

    call getmem3d(b3,1,jx,1,iy,1,nlev*3,'mod_fvgcm:b3')
    call getmem3d(d3,1,jx,1,iy,1,nlev*2,'mod_fvgcm:d3')

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
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
