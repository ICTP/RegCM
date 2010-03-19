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

      subroutine makectl(iunctl,outfil,nx,ny,nzmax,nfld,dsm,clat,clon,  &
                       & plat,plon,iproj,ibigend,truelatl,truelath,     &
                       & grdfac,vmisdat,vnam,lnam,nlev,ntim,varmin,     &
                       & varmax)
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , dsm , grdfac , plat , plon , truelath ,  &
            & truelatl , vmisdat
      integer :: ibigend , iunctl , nfld , nx , ny , nzmax
      character(6) :: iproj
      character(256) :: outfil
      character(64) , dimension(nfld) :: vnam
      character(64) , dimension(nfld) :: lnam
      integer , dimension(nfld) :: nlev , ntim
      real(4) , dimension(3) :: varmax , varmin
      intent (in) clat , clon , dsm , ibigend , iproj , iunctl , lnam , &
                & nfld , nlev , ntim , nx , ny , nzmax , outfil , plat ,&
                & plon , truelath , truelatl , varmax , varmin ,        &
                & vmisdat , vnam
!
! Local variables
!
      character(2) , dimension(12) :: amonth
      real(4) :: ddeg , dlat , dlon , lon1 , lat1
      integer :: i , ifld , itim , nlat , nlon , nvars
      character(64) :: vnam_ctl
      integer :: imi
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
          write (iunctl,99006) nx , ny , clat , clon , nx/2 , ny/2 ,    &
                             & truelatl , truelath , clon , dsm , dsm
        else
          write (iunctl,99007) nx , ny , plon , plat , ddeg ,           &
                             & ddeg*0.95238
        end if
 
        lon1 = varmin(1)
        dlon = dsm*0.001/111./2.
        nlon = 0.5 + nint(abs((varmax(1)-varmin(1))/dlon))
        write (iunctl,99008) nlon , lon1 , dlon
 
        lat1 = varmin(2)
        dlat = dsm*0.001/111./2.
        nlat = 2 + nint(abs(varmax(2)-varmin(2))/dlat)
        write (iunctl,99009) nlat , lat1 , dlat
 
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
          imi = index(vnam(ifld),'MONTHLY')
          if (imi == 0) then
            vnam_ctl = vnam(ifld)
          else
            vnam_ctl = vnam(ifld)(imi+8:)
          end if
          vnam_ctl = trim(vnam_ctl)//trim(amonth(itim))
          write (iunctl,99016) vnam_ctl,nlev(ifld),lnam(ifld)
        end do
      end do
      write (iunctl,99017)
99001 format ('dset ^',a30)
99002 format ('title CLM3/RegCM information')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef ',f8.1)
99006 format ('pdef ',2(i4,x),'lcc ',2(f9.3,x),2(i4,x),5(f9.3,x))
99007 format ('pdef ',2(i4,x),'eta.u',2(x,f9.3),2(x,f9.3))
99008 format ('xdef ',i3,' linear ',f10.5,x,f10.5)
99009 format ('ydef ',i3,' linear ',f10.5,x,f10.5)
99010 format ('xdef ',i3,' linear',2(x,f10.5))
99011 format ('ydef ',i3,' linear',2(x,f10.5))
99012 format ('zdef ',i3,' linear 1 1')
99013 format ('tdef   1 linear 00z15Jan1900 1mo')
99014 format ('vars ',i3)
99015 format (' HEAD                  1 0 header information')
99016 format (a15,x,i2,' 0 ',a40)
99017 format ('endvars')
 
      end subroutine makectl
