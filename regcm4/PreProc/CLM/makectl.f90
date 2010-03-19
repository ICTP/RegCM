      subroutine makectl(iunctl,outfil,nx,ny,nzmax,nfld,dsm,clat,clon,  &
                       & plat,plon,iproj,ibigend,truelatl,truelath,     &
                       & grdfac,vmisdat,vnam,lnam,nlev,ntim,varmin,     &
                       & varmax)
      implicit none
!
! Dummy arguments
!
      real :: clat , clon , dsm , grdfac , plat , plon , truelath ,     &
            & truelatl , vmisdat
      integer :: ibigend , iunctl , nfld , nx , ny , nzmax
      character(6) :: iproj
      character(30) :: outfil
      character(30) , dimension(nfld) :: lnam
      integer , dimension(nfld) :: nlev , ntim
      real , dimension(3) :: varmax , varmin
      character(15) , dimension(nfld) :: vnam
      intent (in) clat , clon , dsm , ibigend , iproj , iunctl , lnam , &
                & nfld , nlev , ntim , nx , ny , nzmax , outfil , plat ,&
                & plon , truelath , truelatl , varmax , varmin ,        &
                & vmisdat , vnam
!
! Local variables
!
      character(2) , dimension(12) :: amonth
      real :: ddeg , dlat , dlon
      integer :: i , ifld , ilat1 , ilon1 , itim , nlat , nlon , nvars
      character(15) :: vnam_ctl
!
      data (amonth(i),i=1,12)/'01' , '02' , '03' , '04' , '05' , '06' , &
           &'07' , '08' , '09' , '10' , '11' , '12'/
 
      ddeg = (dsm/1000.)*(45./atan(1.))/6371229.
 
      write (iunctl,99001) outfil
 
      write (iunctl,99002)
 
      if ( ibigend==1 ) then
        write (iunctl,99003)
      else
        write (iunctl,99004)
      end if
 
      write (iunctl,99005) vmisdat
 
      if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
 
        if ( iproj=='LAMCON' ) then
          write (iunctl,99006) nx , ny , clon , clat , nx/2 , ny/2 ,    &
                             & truelatl , truelath , clon , dsm , dsm
        else if ( iproj=='ROTMER' ) then
          write (iunctl,99007) nx , ny , plon , plat , ddeg ,           &
                             & ddeg*0.95238
        else
        end if
 
        ilon1 = nint(varmin(1))
        nlon = 2*nx
        dlon = (varmax(1)-varmin(1))/nlon
        write (iunctl,99008) nlon , ilon1 , dlon
 
        ilat1 = nint(varmin(2))
        nlat = 2*ny
        dlat = (varmax(2)-varmin(2))/nlat
        write (iunctl,99009) nlon , ilon1 , dlon
 
      else
 
        dlon = (varmax(1)-varmin(1))/nx
        write (iunctl,99010) nx , varmin(1) , dlon
 
        dlat = (varmax(2)-varmin(2))/ny
        write (iunctl,99011) ny , varmin(2) , dlat
 
      end if
 
      write (iunctl,99012) nzmax
 
      write (iunctl,99013)
 
      nvars = 1
      do ifld = 1 , nfld
        nvars = nvars + ntim(ifld)
      end do
      write (iunctl,99014) nvars
 
      write (iunctl,99015)
      do ifld = 1 , nfld
        do itim = 1 , ntim(ifld)
!         print*,ifld,nfld,itim,ntim(nfld)
          vnam_ctl = trim(vnam(ifld))//amonth(itim)
          write (iunctl,99016) vnam_ctl , nlev(ifld) , lnam(ifld)
        end do
      end do
      write (iunctl,99017)
99001 format ('dset ^',a30)
99002 format ('title CLM3/RegCM information')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef ',e10.4)
99006 format ('pdef ',2(i4,x),'lcc ',2(f9.3,x),2(i4,x),5(f9.3,x))
99007 format ('pdef ',2(i4,x),'eta.u',2(x,f9.3),2(x,f9.3))
99008 format ('xdef ',i3,' linear ',i5,x,f10.5)
99009 format ('ydef ',i3,' linear ',i5,x,f10.5)
99010 format ('xdef ',i3,' linear',2(x,f10.5))
99011 format ('ydef ',i3,' linear',2(x,f10.5))
99012 format ('zdef ',i3,' linear 1 1')
99013 format ('tdef   1 linear 00z15Jan1900 1mo')
99014 format ('vars ',i3)
99015 format (' HEAD                  1 0 header information')
99016 format (x,a15,x,i2,' 0 ',a30)
99017 format ('endvars')
 
      end subroutine makectl
