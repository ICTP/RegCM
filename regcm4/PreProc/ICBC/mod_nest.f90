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

      module mod_nest
      use mod_regcm_param , only : ix , jx , kx , ibyte , dattyp
      use mod_preproc_param

      implicit none

      private

      integer , parameter :: np = 15

      integer :: nrec

      real , target , dimension(jx,ix,np*4) :: b3
      real , target , dimension(jx,ix,np*2) :: d3
      real , dimension(jx,ix) :: b3pd

      real , allocatable , target , dimension(:,:,:) :: b2
      real , allocatable , target , dimension(:,:,:) :: d2

      real , allocatable , dimension(:,:,:) :: c , q , t
      real , allocatable , dimension(:,:,:) :: u , v
      real , allocatable , dimension(:,:) :: ps
      real , allocatable , dimension(:,:) :: ht_in , xlat_in , xlon_in

      real , pointer , dimension(:,:,:) :: c3 , h3 , q3 , t3
      real , pointer , dimension(:,:,:) :: u3 , v3
 
      real , pointer , dimension(:,:,:) :: cp , hp , qp , tp
      real , pointer , dimension(:,:,:) :: up , vp

      real , dimension(jx,ix,kx) :: z1

      real , dimension(np) :: plev , sigmar
      real , allocatable , dimension(:) :: sigf
      real , allocatable , dimension(:) :: sig

      character(6) :: iproj_in

      integer :: il , jl , kl , iotyp_in , idate0

      real :: clat_in , clon_in , plat_in , plon_in , ptop_in

      public :: get_nest , headnest

      contains

      subroutine get_nest(idate,ncr)
      use mod_mxncom
      use mod_grid
      use mod_datenum
      use mod_write
      implicit none
!
! Dummy arguments
!
      integer :: idate , ncr
      intent (in) ncr
!
! Local variables
!
      character(14) :: fillin
      integer :: i , idatek , j , k , mn0 , mn1 , nd0 , nd1 , nh0 ,     &
               & nh1 , nmop , ny0 , ny1 , nyrp
      real , dimension(jx,ix) :: d1xa , d1xb , d1xc , d1xd , d1xt
      integer , dimension(jx,ix) :: i1dl , i1dr , i1ul , i1ur , j1dl ,  &
                                  & j1dr , j1ul , j1ur
      real , dimension(jx,ix) :: d2xa , d2xb , d2xc , d2xd , d2xt
      integer , dimension(jx,ix) :: i2dl , i2dr , i2ul , i2ur , j2dl ,  &
                                  & j2dr , j2ul , j2ur
      logical :: there
      real :: wt
!
      if (.not. allocated(b2)) then
        write (*,*) 'Called get_nest before headnest !'
        stop
      end if
!
      if ( idate==idate0 ) then
        write (fillin,99001) idate
        inquire (file='../DATA/RegCM/'//fillin,exist=there)
        if ( .not.there ) then
          write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
          write (*,*) 'please copy (or link)' , fillin
          stop
        end if
        if ( iotyp_in==1 ) then
          open (55,file='../DATA/RegCM/'//fillin,form='unformatted',    &
              & recl=il*jl*ibyte,access='direct')
          nrec = 0
        else if ( iotyp_in==2 ) then
          open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
          rewind (55)
        else
        end if
      else if ( idate==idate1 ) then
        ny0 = idate0/1000000
        mn0 = mod(idate0/10000,100)
        nd0 = mod(idate0/100,100)
        nh0 = mod(idate0,100)
 
        ny1 = idate1/1000000
        mn1 = mod(idate1/10000,100)
        nd1 = mod(idate1/100,100)
        nh1 = mod(idate1,100)
 
        if ( ny0==ny1 .and. mn0==mn1 ) then
          write (fillin,99001) idate0
          inquire (file='../DATA/RegCM/'//fillin,exist=there)
          if ( .not.there ) then
            write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
            write (*,*) 'please copy (or link)' , fillin
            stop
          end if
          if ( iotyp_in==1 ) then
            open (55,file='../DATA/RegCM/'//fillin,form='unformatted',  &
                & recl=il*jl*ibyte,access='direct')
            nrec = ((nd1-nd0)*4+(nh1-nh0)/6)*(kl*6+5)
          else if ( iotyp_in==2 ) then
            open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
            rewind (55)
          else
          end if
        else if ( nd1==1 .and. nh1==0 ) then
          if ( (ny1-ny0)*12+(mn1-mn0)==1 ) then
            write (fillin,99001) idate0
            inquire (file='../DATA/RegCM/'//fillin,exist=there)
            if ( .not.there ) then
              write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
              write (*,*) 'please copy (or link)' , fillin
              stop
            end if
            if ( iotyp_in==1 ) then
              open (55,file='../DATA/RegCM/'//fillin,form='unformatted',&
                  & recl=il*jl*ibyte,access='direct')
              if ( mn0==1 .or. mn0==3 .or. mn0==5 .or. mn0==7 .or.      &
                 & mn0==8 .or. mn0==10 .or. mn0==12 ) then
                nrec = (124-(nd0-1)*4+nh0/6)*(kl*6+5)
              else if ( mn0==4 .or. mn0==6 .or. mn0==9 .or. mn0==11 )   &
                      & then
                nrec = (120-(nd0-1)*4+nh0/6)*(kl*6+5)
              else
                nrec = 112 - (nd0-1)*4 + nh0/6
                if ( mod(ny0,4)==0 ) nrec = nrec + 4
                if ( mod(ny0,100)==0 ) nrec = nrec - 4
                if ( mod(ny0,400)==0 ) nrec = nrec + 4
                nrec = nrec*(kl*6+5)
              end if
            else if ( iotyp_in==2 ) then
              open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
              rewind (55)
            else
            end if
          else
            if ( mn1>1 ) then
              write (fillin,99001) ny1*1000000 + (mn1-1)*10000 + 100
            else
              write (fillin,99001) (ny1-1)*1000000 + 120100
            end if
            inquire (file='../DATA/RegCM/'//fillin,exist=there)
            if ( .not.there ) then
              write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
              write (*,*) 'please copy (or link)' , fillin
              stop
            end if
            if ( iotyp_in==1 ) then
              open (55,file='../DATA/RegCM/'//fillin,form='unformatted',&
                  & recl=il*jl*ibyte,access='direct')
              if ( mn0==1 .or. mn0==3 .or. mn0==5 .or. mn0==7 .or.      &
                 & mn0==8 .or. mn0==10 .or. mn0==12 ) then
                nrec = 123*(kl*6+5)
              else if ( mn0==4 .or. mn0==6 .or. mn0==9 .or. mn0==11 )   &
                      & then
                nrec = 119*(kl*6+5)
              else
                nrec = 111
                if ( mod(ny0,4)==0 ) nrec = nrec + 4
                if ( mod(ny0,100)==0 ) nrec = nrec - 4
                if ( mod(ny0,400)==0 ) nrec = nrec + 4
                nrec = nrec*(kl*6+5)
              end if
            else if ( iotyp_in==2 ) then
              open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
              rewind (55)
            else
            end if
          end if
        else
          write (fillin,99001) ny1*1000000 + mn1*10000 + 100
          inquire (file='../DATA/RegCM/'//fillin,exist=there)
          if ( .not.there ) then
            write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
            write (*,*) 'please copy (or link)' , fillin
            stop
          end if
          if ( iotyp_in==1 ) then
            open (55,file='../DATA/RegCM/'//fillin,form='unformatted',  &
                & recl=il*jl*ibyte,access='direct')
            nrec = ((nd1-1)*4+nh1/6-1)*(kl*6+5)
          else if ( iotyp_in==2 ) then
            open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
            rewind (55)
          else
          end if
        end if
      else
      end if
!     WRITE(*,*) 'Open ATM file:', '../DATA/RegCM/'//fillin
 
      if ( iotyp_in==1 ) then
        if ( idate/=idate1 .and. mod(idate,10000)==100 .and. ncr==1 )   &
           & nrec = nrec - (kl*6+5)
        idatek = idate
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((u(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((v(i,j,k),i=1,jl),j=1,il)
        end do
        nrec = nrec + kl         ! skip omega
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((t(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((q(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((c(i,j,k),i=1,jl),j=1,il)
        end do
        nrec = nrec + 1
        read (55,rec=nrec) ((ps(i,j),i=1,jl),j=1,il)
        nrec = nrec + 4
      else if ( iotyp_in==2 ) then
        if ( idate/=idate1 .and. mod(idate,10000)==100 .and. ncr==1 )   &
           & rewind (55)
 50     continue
        read (55) idatek
        if ( idatek/=idate ) then
          do k = 1 , kl*6 + 5
            read (55)
          end do
!         WRITE(*,*) 'READ IN fields at DATE:',idateK
          go to 50
        end if
!       idate=idateK
 
!       print*,' IDATE = ',idate
        do k = kl , 1 , -1
          read (55) ((u(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          read (55) ((v(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          read (55)
        end do
        do k = kl , 1 , -1
          read (55) ((t(i,j,k),i=1,jl),j=1,il)
        end do
        do k = kl , 1 , -1
          read (55) ((q(i,j,k),i=1,jl),j=1,il)
        end do
        do k = 1 , kl
          read (55) ((c(i,j,k),i=1,jl),j=1,il)
        end do
        read (55) ((ps(i,j),i=1,jl),j=1,il)
        do k = 1 , 4
          read (55)
        end do
      else
      end if
      write (*,*) 'READ IN fields at DATE:' , idatek , ' from ' , fillin
      if ( idate/=idate1 .and. mod(idate,10000)==100 .and. ncr==1 ) then
        write (fillin,99001) idate
        inquire (file='../DATA/RegCM/'//fillin,exist=there)
        if ( .not.there ) then
          write (*,*) '../DATA/RegCM/'//fillin , ' is not available'
          write (*,*) 'please copy (or link)' , fillin
          stop
        end if
        if ( iotyp_in==1 ) then
          open (55,file='../DATA/RegCM/'//fillin,form='unformatted',    &
              & recl=il*jl*ibyte,access='direct')
          nrec = 0
        else if ( iotyp_in==2 ) then
          open (55,file='../DATA/RegCM/'//fillin,form='unformatted')
          rewind (55)
        else
        end if
!       WRITE(*,*) 'Open ATM file:', fillin
      end if
!
!     to calculate Heights on sigma surfaces.
      call htsig_o(t,z1,ps,ht_in,sig,ptop_in,jl,il,kl)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
      call height_o(hp,z1,t,ps,ht_in,sig,ptop_in,jl,il,kl,    &
                  & plev,np)
!     2. For Zonal and Meridional Winds
      call intlin_o(up,u,ps,sig,ptop_in,jl,il,kl,plev,np)
      call intlin_o(vp,v,ps,sig,ptop_in,jl,il,kl,plev,np)
!     3. For Temperatures
      call intlog_o(tp,t,ps,sig,ptop_in,jl,il,kl,plev,np)
!     4. For Moisture qva & qca
      call humid1_o(t,q,ps,sig,ptop_in,jl,il,kl)
      call intlin_o(qp,q,ps,sig,ptop_in,jl,il,kl,plev,np)
      call intlog_o(cp,c,ps,sig,ptop_in,jl,il,kl,plev,np)
      call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in,grdfac,jl,    &
                  & il,np,plon_in,plat_in,iproj_in)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,ix,i1ur,i1ul,    &
                  & i1dr,i1dl,j1ur,j1ul,j1dr,j1dl,d1xt,d1xa,d1xb,d1xc,  &
                  & d1xd,jl,il,np)
      call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,ix,i2ur,i2ul,    &
                  & i2dr,i2dl,j2ur,j2ul,j2dr,j2dl,d2xt,d2xa,d2xb,d2xc,  &
                  & d2xd,jl,il,np)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,ix,np,plon,plat,  &
                & iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
      call top2btm(t3,jx,ix,np)
      call top2btm(q3,jx,ix,np)
      call top2btm(c3,jx,ix,np)
      call top2btm(h3,jx,ix,np)
      call top2btm(u3,jx,ix,np)
      call top2btm(v3,jx,ix,np)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,ix,np)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,ix)
      call p1p2(b3pd,ps4,jx,ix)
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,ix,np)
 
      if ( ssttyp=='EH5RF' .or. ssttyp=='EH5A2' .or.                    &
         & ssttyp=='EH5B1' .or. ssttyp=='EHA1B' ) then
        call mksst3(ts4,sst1,topogm,xlandu,jx,ix,idate)
      else if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!       F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!       PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp=='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,ix,nyrp, &
               &      nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,ix,nyrp,nmop,wt)
        end if
      else
        if ( ssttyp=='OI2WK' ) then
          call mksst2a(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,ix,    &
               &       idate/100)
        else
          call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,ix,idate/100)
        end if
      end if
 
!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,ix,kx,np)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,ix,kx,np)
!
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,ix,kx,np)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,ix,kx,np)
      call humid2fv(t4,q4,ps4,ptop,sigma2,jx,ix,kx)
      call intv1(c4,c3,ps4,sigma2,sigmar,ptop,jx,ix,kx,np)
!
!     F4     DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,ix,kx)
!
!     G      WRITE AN INITIAL FILE FOR THE RegCM
      call writef2(ptop,idate)
!
      deallocate(sigf)
      deallocate(sig)
      deallocate(b2)
      deallocate(d2)
      deallocate(c)
      deallocate(q)
      deallocate(t)
      deallocate(u)
      deallocate(v)
      deallocate(ps)
      deallocate(xlat_in)
      deallocate(xlon_in)
      deallocate(ht_in)

99001 format ('ATM.',i10)
!
      end subroutine get_nest
!
!
!
      subroutine headnest
      use mod_mxncom
      use mod_grid
      implicit none
!
! Local variables
!
      real :: d2r , dtb , dtc , dto , dtr , dxsp , ptsp , xsign ,       &
            & truelat1 , truelat2
      integer :: ibltyp , iboudy , icup , ipptls , k
      integer :: ias
      logical :: there
!
      plev(1) = 50.
      plev(2) = 70.
      plev(3) = 100.
      plev(4) = 150.
      plev(5) = 200.
      plev(6) = 250.
      plev(7) = 300.
      plev(8) = 400.
      plev(9) = 500.
      plev(10) = 600.
      plev(11) = 700.
      plev(12) = 775.
      plev(13) = 850.
      plev(14) = 925.
      plev(15) = 1000.
 
      do k = 1 , np
        sigmar(k) = plev(k)*0.001
      end do
 
      inquire (file='../DATA/RegCM/OUT_HEAD',exist=there)
      if ( .not.there ) then
        write (*,*) '../DATA/RegCM/OUT_HEAD is not available'
        write (*,*) 'please copy (or link) the previous output OUT_HEAD'
        stop
      end if

      open (49,file='../DATA/RegCM/OUT_HEAD',form='unformatted',        &
           &access='direct',recl=24**ibyte)
      read (49,rec=1) idate0 , ibltyp , icup , ipptls , iboudy , il ,   &
                    & jl , kl
      close (49)

!     Reserve space for I/O

      allocate(sigf(kl+1), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: sigf'
      allocate(sig(kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: sig'
      allocate(b2(il,jl,np*4), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: b2'
      allocate(d2(il,jl,np*2), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: d2'
      allocate(c(il,jl,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: c'
      allocate(q(il,jl,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: q'
      allocate(t(il,jl,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: t'
      allocate(u(il,jl,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: u'
      allocate(v(il,jl,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: v'
      allocate(ps(il,jl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: ps'
      allocate(xlat_in(il,jl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: xlat_in'
      allocate(xlon_in(il,jl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: xlon_in'
      allocate(ht_in(il,jl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: ht_in'

      open (49,file='../DATA/RegCM/OUT_HEAD',form='unformatted',        &
           &access='direct',recl=il*jl*ibyte)
      read (49,rec=1) idate0 , ibltyp , icup , ipptls , iboudy , il ,   &
                    & jl , kl , (sigf(k),k=kl+1,1,-1) , dxsp ,          &
                    & ptsp , clat_in , clon_in , plat_in , plon_in ,    &
                    & iproj_in , dto , dtb , dtr , dtc ,   &
                    & iotyp_in , truelat1 , truelat2
      ptop_in = ptsp*10.
      if ( iproj_in=='LAMCON' ) then
        d2r = atan(1.)*4./180.
        if ( clat_in<0. ) then
          xsign = -1.       ! SOUTH HEMESPHERE
        else
          xsign = 1.        ! NORTH HEMESPHERE
        end if
        if ( abs(truelat1-truelat2)>1.E-1 ) then
          grdfac = (alog10(cos(truelat1*d2r))                           &
                   & -alog10(cos(truelat2*d2r)))                        &
                   & /(alog10(tan((45.0-xsign*truelat1/2.0)*d2r))       &
                   & -alog10(tan((45.0-xsign*truelat2/2.0)*d2r)))
        else
          grdfac = xsign*sin(truelat1*d2r)
        end if
      else if ( iproj_in=='POLSTR' ) then
        grdfac = 1.0
      else if ( iproj_in=='NORMER' ) then
        grdfac = 0.0
      else
        grdfac = 0.0
      end if
      read (49,rec=2) ht_in
      read (49,rec=6) xlat_in
      read (49,rec=7) xlon_in
      close (49)
 
      do k = 1 , kl
        sig(k) = 0.5*(sigf(k)+sigf(k+1))
      end do
 
      imxmn = 0
      lcross = 0
      ldot = 0

!     Set up pointers
 
      tp => b2(:,:,1:np)
      qp => b2(:,:,np+1:2*np)
      cp => b2(:,:,2*np+1:3*np)
      hp => b2(:,:,3*np+1:4*np)
      up => d2(:,:,1:np)
      vp => d2(:,:,np+1:2*np)
      t3 => b3(:,:,1:np)
      q3 => b3(:,:,np+1:2*np)
      c3 => b3(:,:,2*np+1:3*np)
      h3 => b3(:,:,3*np+1:4*np)
      u3 => d3(:,:,1:np)
      v3 => d3(:,:,np+1:2*np)

      end subroutine headnest

      end module mod_nest
