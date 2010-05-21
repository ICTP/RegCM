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

      module mod_nest
      use mod_dynparam

      implicit none

      private

      integer , parameter :: np = 15

      integer :: nrec

      real(4) , allocatable , target , dimension(:,:,:) :: b3
      real(4) , allocatable , target , dimension(:,:,:) :: d3
      real(4) , allocatable , dimension(:,:) :: b3pd
      real(4) , allocatable , dimension(:,:,:) :: z1

      real(4) , allocatable , target , dimension(:,:,:) :: b2
      real(4) , allocatable , target , dimension(:,:,:) :: d2

      real(4) , allocatable , dimension(:,:,:) :: c , q , t
      real(4) , allocatable , dimension(:,:,:) :: u , v
      real(4) , allocatable , dimension(:,:) :: ps
      real(4) , allocatable , dimension(:,:) :: ht_in
      real(4) , allocatable , dimension(:,:) :: xlat_in , xlon_in

      real(4) , pointer , dimension(:,:,:) :: c3 , h3 , q3 , t3
      real(4) , pointer , dimension(:,:,:) :: u3 , v3
 
      real(4) , pointer , dimension(:,:,:) :: cp , hp , qp , tp
      real(4) , pointer , dimension(:,:,:) :: up , vp

      real(4) , dimension(np) :: plev , sigmar
      real(4) , allocatable , dimension(:) :: sigf
      real(4) , allocatable , dimension(:) :: sig

      character(6) :: iproj_in

      integer :: iy_in , jx_in , kl , iotyp_in , idate0

      real(4) :: clat_in , clon_in , plat_in , plon_in , ptop_in

      public :: get_nest , headnest

      contains

      subroutine get_nest(idate,ncr)
      use mod_grid
      use mod_date
      use mod_write
      use mod_interp , only : cressmcr , cressmdt
      use mod_vertint
      use mod_hgt
      use mod_humid
      use mod_mksst
      use mod_uvrot
      use mod_vectutil
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
      character(256) :: inpfile
      integer :: i , idatek , j , k , mn0 , mn1 , nd0 , nd1 , nh0 ,     &
               & nh1 , nmop , ny0 , ny1 , nyrp
      logical :: there
      real(4) :: wt
!
      if (.not. allocated(b2)) then
        write (*,*) 'Called get_nest before headnest !'
        stop
      end if
!
      if ( idate==idate0 ) then
        write (fillin,99001) idate
        inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
        inquire (file=inpfile,exist=there)
        if ( .not.there ) then
          write (*,*) trim(inpfile), ' is not available'
          write (*,*) 'please copy (or link)' , trim(inpfile)
          stop
        end if
        if ( iotyp_in==1 ) then
          open (55,file=trim(inpfile),form='unformatted',               &
              & recl=iy_in*jx_in*ibyte,access='direct')
          nrec = 0
        else if ( iotyp_in==2 ) then
          open (55,file=trim(inpfile),form='unformatted')
          rewind (55)
        else
        end if
      else if ( idate==globidate1 ) then
        ny0 = idate0/1000000
        mn0 = mod(idate0/10000,100)
        nd0 = mod(idate0/100,100)
        nh0 = mod(idate0,100)
 
        ny1 = globidate1/1000000
        mn1 = mod(globidate1/10000,100)
        nd1 = mod(globidate1/100,100)
        nh1 = mod(globidate1,100)
 
        if ( ny0==ny1 .and. mn0==mn1 ) then
          write (fillin,99001) idate0
          inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
          inquire (file=trim(inpfile),exist=there)
          if ( .not.there ) then
            write (*,*) trim(inpfile), ' is not available'
            write (*,*) 'please copy (or link)' , trim(inpfile)
            stop
          end if
          if ( iotyp_in==1 ) then
            open (55,file=trim(inpfile),form='unformatted',             &
                & recl=iy_in*jx_in*ibyte,access='direct')
            nrec = ((nd1-nd0)*4+(nh1-nh0)/6)*(kl*6+5)
          else if ( iotyp_in==2 ) then
            open (55,file=trim(inpfile),form='unformatted')
            rewind (55)
          else
          end if
        else if ( nd1==1 .and. nh1==0 ) then
          if ( (ny1-ny0)*12+(mn1-mn0)==1 ) then
            write (fillin,99001) idate0
            inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
            inquire (file=trim(inpfile),exist=there)
            if ( .not.there ) then
              write (*,*) trim(inpfile), ' is not available'
              write (*,*) 'please copy (or link)' , trim(inpfile)
              stop
            end if
            if ( iotyp_in==1 ) then
              open (55,file=trim(inpfile),form='unformatted',           &
                 &  recl=iy_in*jx_in*ibyte,access='direct')
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
              open (55,file=trim(inpfile),form='unformatted')
              rewind (55)
            else
            end if
          else
            if ( mn1>1 ) then
              write (fillin,99001) ny1*1000000 + (mn1-1)*10000 + 100
            else
              write (fillin,99001) (ny1-1)*1000000 + 120100
            end if
            inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
            inquire (file=trim(inpfile),exist=there)
            if ( .not.there ) then
              write (*,*) trim(inpfile), ' is not available'
              write (*,*) 'please copy (or link)' , trim(inpfile)
              stop
            end if
            if ( iotyp_in==1 ) then
              open (55,file=trim(inpfile),form='unformatted',           &
                  & recl=iy_in*jx_in*ibyte,access='direct')
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
              open (55,file=trim(inpfile),form='unformatted')
              rewind (55)
            else
            end if
          end if
        else
          write (fillin,99001) ny1*1000000 + mn1*10000 + 100
          inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
          inquire (file=trim(inpfile),exist=there)
          if ( .not.there ) then
            write (*,*) trim(inpfile), ' is not available'
            write (*,*) 'please copy (or link)' , trim(inpfile)
            stop
          end if
          if ( iotyp_in==1 ) then
            open (55,file=trim(inpfile),form='unformatted',             &
                & recl=iy_in*jx_in*ibyte,access='direct')
            nrec = ((nd1-1)*4+nh1/6-1)*(kl*6+5)
          else if ( iotyp_in==2 ) then
            open (55,file=trim(inpfile),form='unformatted')
            rewind (55)
          else
          end if
        end if
      else
      end if

      ! write (6,*) 'Open ATM file: ', trim(inpfile)
 
      if ( iotyp_in==1 ) then
        if ( idate/=globidate1 .and. mod(idate,10000)==100 .and.        &
           & ncr==1 ) nrec = nrec - (kl*6+5)
        idatek = idate
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((u(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((v(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        nrec = nrec + kl         ! skip omega
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((t(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((q(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          nrec = nrec + 1
          read (55,rec=nrec) ((c(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        nrec = nrec + 1
        read (55,rec=nrec) ((ps(i,j),i=1,jx_in),j=1,iy_in)
        nrec = nrec + 4
      else if ( iotyp_in==2 ) then
        if ( idate/=globidate1 .and. mod(idate,10000)==100 .and.        &
           & ncr==1 ) rewind (55)
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
          read (55) ((u(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          read (55) ((v(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          read (55)
        end do
        do k = kl , 1 , -1
          read (55) ((t(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = kl , 1 , -1
          read (55) ((q(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        do k = 1 , kl
          read (55) ((c(i,j,k),i=1,jx_in),j=1,iy_in)
        end do
        read (55) ((ps(i,j),i=1,jx_in),j=1,iy_in)
        do k = 1 , 4
          read (55)
        end do
      else
      end if
      write (*,*) 'READ IN fields at DATE:' , idatek , ' from ' , fillin

      if ( idate/=globidate1 .and. mod(idate,10000)==100 .and. ncr==1 ) then
        write (fillin,99001) idate
        inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep//fillin
        inquire (file=trim(inpfile),exist=there)
        if ( .not.there ) then
          write (*,*) trim(inpfile), ' is not available'
          write (*,*) 'please copy (or link)' , trim(inpfile)
          stop
        end if
        if ( iotyp_in==1 ) then
          open (55,file=trim(inpfile),form='unformatted',               &
              & recl=iy_in*jx_in*ibyte,access='direct')
          nrec = 0
        else if ( iotyp_in==2 ) then
          open (55,file=trim(inpfile),form='unformatted')
          rewind (55)
        else
        end if
!       WRITE(*,*) 'Open ATM file:', trim(inpfile)
      end if
!
!     to calculate Heights on sigma surfaces.
      call htsig_o(t,z1,ps,ht_in,sig,ptop_in,jx_in,iy_in,kl)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
      call height_o(hp,z1,t,ps,ht_in,sig,ptop_in,jx_in,iy_in,kl,    &
                  & plev,np)
!     2. For Zonal and Meridional Winds
      call intlin_o(up,u,ps,sig,ptop_in,jx_in,iy_in,kl,plev,np)
      call intlin_o(vp,v,ps,sig,ptop_in,jx_in,iy_in,kl,plev,np)
!     3. For Temperatures
      call intlog_o(tp,t,ps,sig,ptop_in,jx_in,iy_in,kl,plev,np)
!     4. For Moisture qva & qca
      call humid1_o(t,q,ps,sig,ptop_in,jx_in,iy_in,kl)
      call intlin_o(qp,q,ps,sig,ptop_in,jx_in,iy_in,kl,plev,np)
      call intlog_o(cp,c,ps,sig,ptop_in,jx_in,iy_in,kl,plev,np)
      call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in,grdfac,       &
             &      jx_in,iy_in,np,plon_in,plat_in,iproj_in)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np)
      call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,iy,              &
                  & jx_in,iy_in,np)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,np,plon,plat,  &
                & iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
      call top2btm(t3,jx,iy,np)
      call top2btm(q3,jx,iy,np)
      call top2btm(c3,jx,iy,np)
      call top2btm(h3,jx,iy,np)
      call top2btm(u3,jx,iy,np)
      call top2btm(v3,jx,iy,np)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,np)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call p1p2(b3pd,ps4,jx,iy)
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,np)
 
      if ( ssttyp=='EH5RF' .or. ssttyp=='EH5A2' .or.                    &
         & ssttyp=='EH5B1' .or. ssttyp=='EHA1B' ) then
        call mksst3(ts4,sst1,topogm,xlandu,jx,iy,idate)
      else if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!       F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!       PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp=='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,nyrp, &
               &      nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
        end if
      else
        if ( ssttyp=='OI2WK' ) then
          call mksst2a(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,    &
               &       idate/100)
        else
          call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,iy,idate/100)
        end if
      end if
 
!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
!
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
      call humid2fv(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
      call intv1(c4,c3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
!
!     F4     DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
!     G      WRITE AN INITIAL FILE FOR THE RegCM
      call writef2(ptop,idate)
!
99001 format ('ATM.',i10)
!
      end subroutine get_nest
!
!
!
      subroutine headnest
      use mod_grid
      use mod_interp, only : imxmn , lcross , ldot
      use mod_constants , only : degrad
      implicit none
!
! Local variables
!
      real(4) :: dtb , dtc , dto , dtr , dxsp , ptsp , xsign ,          & 
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
 
      inquire (file=trim(inpglob)//'/RegCM/OUT_HEAD',exist=there)
      if ( .not.there ) then
        write (*,*) trim(inpglob)//'/RegCM/OUT_HEAD is not available'
        write (*,*) 'please copy (or link) the previous output OUT_HEAD'
        stop
      end if

      open (49,file=trim(inpglob)//'/RegCM/OUT_HEAD',                   &
          & form='unformatted',access='direct',recl=24**ibyte)
      read (49,rec=1) idate0 , ibltyp , icup , ipptls , iboudy , iy_in ,&
                    & jx_in , kl
      close (49)

      iy_in = iy_in - 2
      jx_in = jx_in - 2

!     Reserve space for I/O

      allocate(sigf(kl+1), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: sigf'
      allocate(sig(kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: sig'
      allocate(b2(jx_in,iy_in,np*4), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: b2'
      allocate(d2(jx_in,iy_in,np*2), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: d2'
      allocate(c(jx_in,iy_in,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: c'
      allocate(q(jx_in,iy_in,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: q'
      allocate(t(jx_in,iy_in,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: t'
      allocate(u(jx_in,iy_in,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: u'
      allocate(v(jx_in,iy_in,kl), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: v'
      allocate(ps(jx_in,iy_in), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: ps'
      allocate(xlat_in(jx_in,iy_in), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: xlat_in'
      allocate(xlon_in(jx_in,iy_in), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: xlon_in'
      allocate(ht_in(jx_in,iy_in), stat=ias)
      if (ias /= 0) stop 'Allocation Error in headnest: ht_in'

      open (49,file=trim(inpglob)//'/RegCM/OUT_HEAD',form='unformatted',&
          & access='direct',recl=iy_in*jx_in*ibyte)
      read (49,rec=1) idate0 , ibltyp , icup , ipptls , iboudy , iy_in ,&
                    & jx_in , kl , (sigf(k),k=kl+1,1,-1) , dxsp ,       &
                    & ptsp , clat_in , clon_in , plat_in , plon_in ,    &
                    & iproj_in , dto , dtb , dtr , dtc ,   &
                    & iotyp_in , truelat1 , truelat2
      ptop_in = ptsp*10.
      iy_in = iy_in - 2
      jx_in = jx_in - 2
      if ( iproj_in=='LAMCON' ) then
        if ( clat_in<0. ) then
          xsign = -1.       ! SOUTH HEMESPHERE
        else
          xsign = 1.        ! NORTH HEMESPHERE
        end if
        if ( abs(truelat1-truelat2)>1.E-1 ) then
          grdfac = (log10(cos(truelat1*degrad))                         &
                   & -log10(cos(truelat2*degrad)))                      &
                   & /(log10(tan((45.0-xsign*truelat1/2.0)*degrad))     &
                   & -log10(tan((45.0-xsign*truelat2/2.0)*degrad)))
        else
          grdfac = xsign*sin(truelat1*degrad)
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

      if (allocated(b3)) deallocate(b3)
      if (allocated(d3)) deallocate(d3)
      if (allocated(b3pd)) deallocate(b3pd)
      if (allocated(z1)) deallocate(z1)
      allocate(b3(iy_in,jx_in,np*4))
      allocate(d3(iy_in,jx_in,np*2))
      allocate(b3pd(iy_in,jx_in))
      allocate(z1(iy_in,jx_in,kl))

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
