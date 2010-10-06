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
!
! DUST module
!
      use mod_constants
      use mod_dynparam
      use mod_ncio
      use mod_trachem
      use mod_message
#ifdef MPP1
      use mod_mppio
#endif

      implicit none
!
      integer , parameter :: nsoil = 152
      integer , parameter :: mode = 3
      integer , parameter :: jsoilm = 1 
      integer , parameter :: jfs =1 
      integer , parameter :: ust = 1
!
! lognormal alfaro parameters
! define the aerosol distribution at the emission and the corresponding 
! weighting factors in fuction of bin sizes
!
      integer, parameter :: isize = 12
      real(8), parameter :: d1 =1.5
      real(8), parameter :: d2 =6.7
      real(8), parameter :: d3 =14.2
      real(8), parameter :: sigma1=1.7 
      real(8), parameter :: sigma2 =1.2
      real(8), parameter :: sigma3 =1.5      

      ! corresponding binding energies/
      real(8),parameter  :: e1 = 3.61
      real(8),parameter  :: e2 = 3.52
      real(8),parameter  :: e3 = 3.46
!        
!     Basic dust aerosol density (ACE-2 ) in kg/m3
!        
      real(8) , parameter :: rhop = 2650.0

      real(8), dimension (2,isize)  ::  aerosize       
      real(8) , allocatable, dimension(:) ::  frac1,frac2,frac3

! soil variable, srel 2 d corresponds to the soil aggregate size distribution
! in each texture type.

      real(8) , allocatable, dimension(:,:,:) :: clay2row2 , sand2row2 ,&
           & silt2row2
      real(8) , allocatable,  dimension(:,:) :: clayrow2 ,   &
           & sandrow2
      real(8) , allocatable,  dimension(:,:,:,:) :: srel2d
      real(8) , allocatable , dimension(:,:,:) ::dustsotex
      real(8) , dimension(nsoil) :: dp_array  ! Name of variable changed ! SC. 06.10.2010
!
! Initialise sub bin aerosol distribution 
!
      data  aerosize/1.0E-08 , 2.0E-08 , 2.0E-08 , 4.0E-08 , 4.0E-08 ,   &
         & 8.0E-08 , 8.0E-08 , 1.6E-07 , 1.6E-07 , 3.2E-07 , 3.2E-07 ,  &
         & 6.4E-07 , 6.4E-07 , 1.28E-06 , 1.28E-06 , 2.56E-06 ,         &
         & 2.56E-06 , 5.12E-06 , 5.12E-06 , 10.4E-06 , 10.24E-06 ,      &
         & 20.48E-06 , 20.48E-06 , 40.6E-06/
!
      contains 
!
      subroutine allocate_mod_dust 
        implicit none

        allocate(frac1(isize))
        allocate(frac2(isize))
        allocate(frac3(isize)) 
        frac1 = 0.0D0
        frac2 = 0.0D0
        frac3 = 0.0D0
#ifdef MPP1
        allocate(clay2row2(iy,nats,jxp))
        allocate(sand2row2(iy,nats,jxp))
        allocate(silt2row2(iy,nats,jxp))
        allocate(clayrow2(iy,jxp))
        allocate(dustsotex(iy,jxp,nats))
        allocate(sandrow2(iy,jxp))
        allocate(srel2d(iy,jxp,nsoil,nats))
#else        
        allocate(clay2row2(iy,nats,jx))
        allocate(sand2row2(iy,nats,jx))
        allocate(silt2row2(iy,nats,jx))
        allocate(clayrow2(iy,jx))
        allocate(dustsotex(iy,jx,nats))
        allocate(sandrow2(iy,jx))
        allocate(srel2d(iy,jx,nsoil,nats))
#endif 
        clay2row2 = 0.0D0
        sand2row2 = 0.0D0
        silt2row2 = 0.0D0
        clayrow2 = 0.0D0
        dustsotex = 0.0D0
        sandrow2 = 0.0D0
        srel2d = 0.0D0

      end subroutine allocate_mod_dust
!
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
!
      subroutine inidust
! 
#ifdef MPP1
#ifdef IBM
      include 'mpif.h'
#else
      use mpi
#endif
#endif
      implicit none
!
#ifdef MPP1
      integer :: ierr
#endif
!
      real(8) , dimension(nats) :: bcly , bslt , bsnd
      real(8) :: deldp , eps , rhop , stotal , xk , xl , xm , xn
      integer :: i , j , n , nm , ns , nt , itr
      real(8) , dimension(3,12) :: mmd , pcent , sigma
      real(8) , dimension(iy,nsoil,nats) :: srel
      real(8) , dimension(nsoil) :: ss
      logical :: rd_tex 
      character(5) :: aerctl
      real(8) :: alogdi , amean1 , amean2 , amean3 , asigma1 ,          &
               & asigma2 , asigma3 , rwi , totv1 , totv2 , totv3
!
! Fab update 
! change type 1 and 2 and 3 to Laurent et al., 2008, marticorena et al., 1997
!  soil size parameter 
!
      data bcly/0.00 , 0.10 , 0.10 , 0.15 , 0.15 , 0.15 , 0.20 , 0.20 , &
         & 0.30 , 0.35 , 0.40 , 0.50/
      data bsnd/0.90 , 0.60 , 0.60 , 0.50 , 0.45 , 0.35 , 0.30 , 0.30 , &
         & 0.20 , 0.65 , 0.60 , 0.50/
      data bslt/0.10 , 0.30 , 0.30 , 0.35 , 0.40 , 0.50 , 0.50 , 0.50 , &
         & 0.50 , 0.00 , 0.00 , 0.00/
! 
      data rhop/2650.000/
      data eps/1.0E-7/
!
      data mmd/690.0 , 210.0 , 10.0 , 690.0 , 210.0 , 0.0 , 690.0 ,     &
         & 210.0 , 0.0 , 520.0 , 100.0 , 5.0 , 520.0 , 75.0 , 2.5 ,     &
         & 520.0 , 75.0 , 2.5 , 210.0 , 75.0 , 2.5 , 210.0 , 50.0 ,     &
         & 2.5 , 125.0 , 50.0 , 1.0 , 100.0 , 10.0 , 1.0 , 100.0 ,      &
         & 10.0 , 0.5 , 100.0 , 10.0 , 0.5/
! 
      data sigma/1.6 , 1.8 , 1.8 , 1.6 , 1.8 , 1.8 , 1.6 , 1.8 , 1.8 ,  &
         & 1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 , 1.6 , 1.7 , 1.8 , 1.7 ,  &
         & 1.7 , 1.8 , 1.7 , 1.7 , 1.8 , 1.7 , 1.7 , 1.8 , 1.8 , 1.8 ,  &
         & 1.8 , 1.8 , 1.8 , 1.8 , 1.8 , 1.8 , 1.8/
! 
       data pcent/1.0 , 0.00 , 0.00 , 0.90 , 0.10 , 0.00 , 0.80 , 0.20 ,&
         & 0.00 , 0.50 , 0.35 , 0.15 , 0.45 , 0.40 , 0.15 , 0.35 ,      &
         & 0.50 , 0.15 , 0.30 , 0.50 , 0.20 , 0.30 , 0.50 , 0.20 ,      &
         & 0.20 , 0.50 , 0.30 , 0.65 , 0.00 , 0.35 , 0.60 , 0.00 ,      &
         & 0.40 , 0.50 , 0.00 , 0.50/
!
      rd_tex = .false.
#ifdef MPP1
      if ( myid.eq.0 ) then
        do itr = 1 , ntr
          aerctl = chtrname(itr)
          if ( aerctl(1:4).eq.'DUST' ) then
            rd_tex = .true.
            exit
          end if
        end do
      end if
    call mpi_bcast(rd_tex,1,mpi_logical,0,mpi_comm_world,ierr)
#else
      do itr = 1 , ntr
        aerctl = chtrname(itr)
        if ( aerctl(1:4).eq.'DUST' ) then
          rd_tex = .true.
          exit
        end if
      end do
#endif

#ifdef MPP1
#ifdef CLM
!      if ( myid.eq.0 ) then
!        if ( rd_tex ) then
!          call clm_getsoitex()
!          do j = 1 , jx
!            do i = 1 , iy
!              dustsotex_io(i,j) = clm_soitex(i,j)
!            end do
!          end do
!        end if
!      end if
!      if ( allocated(clm_soitex) ) deallocate(clm_soitex)
      if ( myid.eq.0 ) then
        if ( rd_tex ) then
          call read_texture(nats,dustsotex_io)
        end if
      end if
#else
      if ( myid.eq.0 ) then
        if ( rd_tex ) then
          call read_texture(nats,dustsotex_io)
        end if
      end if
#endif
      if (myid.eq.0 ) then
        do j=1,jx
          do n=1,nats
            do i=1,iy
              src_1(i,n,j)=dustsotex_io(i,j,n)
            end do
          end do
        end do
      end if

      call mpi_scatter(src_1,iy*jxp*nats,mpi_real8, &
                     & src1,iy*jxp*nats,mpi_real8,0,  &
                     & mpi_comm_world,ierr)


        do j=1,jendl
          do n=1,nats
            do i=1,iy
              dustsotex(i,j,n)=src1(i,n,j)
            end do
        end do
      end do

#else
      if ( rd_tex ) then
        call read_texture(nats,dustsotex)
      end if
#endif
    
! end texture file reading

      clay2row2 = 0.0
      clayrow2  = 0.0
      sand2row2 = 0.0
      silt2row2 = 0.0
      srel2d    = 0.0

#ifdef MPP1
      do j = jbegin , jendm
#else

#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm2
#endif
#endif
       do i = 2 , iym2
         do nt=1,nats
            clay2row2(i,nt,j) = bcly(nt)*100.0
            sand2row2(i,nt,j) = bsnd(nt)*100.0
            silt2row2(i,nt,j) = bslt(nt)*100.0

            sandrow2(i,j) = sandrow2(i,j) + dustsotex(i,j,nt)* sand2row2(i,nt,j)
            clayrow2(i,j) = clayrow2(i,j) + dustsotex(i,j,nt)*clay2row2(i,nt,j)

          end do
        end do

      end do ! end j loop


      dp_array(1) = 0.0001  !cm
      do ns = 2 , nsoil
        dp_array(ns) = dp_array(ns-1)*exp(0.0460517018598807)
        deldp = dp_array(ns) - dp_array(ns-1)
      end do
 
#ifdef MPP1
      do j = jbegin , jendm
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm2
#endif
#endif
        srel(:,:,:) = 0.
        do i = 2 , iym2
          do nt = 1 , nats
            ss(:) =0.0
            stotal = 0.0
            if ( sand2row2(i,nt,j).gt.0. ) then
              do ns = 1 , nsoil          !soil size segregatoin no
                do nm = 1 , mode       !soil mode = 3
                  if ( (pcent(nm,nt).gt.eps) .and.                    &
                       & (sigma(nm,nt).ne.0.0) ) then
                    xk = pcent(nm,nt)/(sqrt(twopi)*log(sigma(nm,nt)))
                    xl = ((log(dp_array(ns))-log(mmd(nm,nt)*1.E-4))**2)     &
                         & /(2.0*(log(sigma(nm,nt)))**2)
                    xm = xk*exp(-xl)
                  else
                    xm = 0.0
                  end if
                  xn = rhop*(2.0/3.0)*(dp_array(ns)/2.0)
                  deldp = 0.0460517018598807
                                              ! dp_array(2)-dp_array(1) ss(nsoil)
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
        end do

        do nt=1,nats
          do i = 2 , iym2
            do ns = 1 , nsoil
              srel2d(i,j,ns,nt) = srel(i,ns,nt)
            end do
          end do
        end do

      end do  ! end J loop
!
! Finally calculate the emission stribution weights in function of Alfaro
! lognormal distribution parameters 
!
      totv1 = 0.0
      totv2 = 0.0
      totv3 = 0.0

      do n = 1 , isize

        rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.E6
        alogdi = log10(rwi)
        amean1 = log10(d1)
        amean2 = log10(d2)
        amean3 = log10(d3)

        asigma1 = log10(sigma1)
        asigma2 = log10(sigma2)
        asigma3 = log10(sigma3)

        frac1(n) = exp(-(alogdi-amean1)**2./(2*asigma1**2))
        frac2(n) = exp(-(alogdi-amean2)**2./(2*asigma2**2))
        frac3(n) = exp(-(alogdi-amean3)**2./(2*asigma3**2))

        totv1 = totv1 + frac1(n)
        totv2 = totv2 + frac2(n)
        totv3 = totv3 + frac3(n)

      end do

      do n = 1 , isize
        frac1(n) = frac1(n)/totv1
        frac2(n) = frac2(n)/totv2
        frac3(n) = frac3(n)/totv3
        if ( frac1(n).lt.1.E-9 ) frac1(n) = 0.0
        if ( frac2(n).lt.1.E-9 ) frac2(n) = 0.0
        if ( frac3(n).lt.1.E-9 ) frac3(n) = 0.0
      end do

      end subroutine inidust
! 
      function cvmgt(val1,val2,cond)
      implicit none
!
      logical :: cond
      real(8) :: val1 , val2
      real(8) :: cvmgt
      intent (in) cond , val1 , val2
!
      if ( cond ) then
        cvmgt = val1
      else
        cvmgt = val2
      end if
!
      end function cvmgt
!
!     *****************************************************************
!     * calculate of ustar01(d) using iversen and white (1982)     ****
!     * for smoth surface:                                         ****
!     * coded by :                                                 ****
!     * ashraf s. zakey, 2003                                      ****
!     * dum    : particle diameter [um]                            ****
!     * ustar0 : threshold frication velocity [m/s]                ****
!     *****************************************************************
!
      function ustart01(rhop,dum,rhair)
      implicit none
!
      real(8) , parameter :: a2 = 0.129 , c1 = 0.006 , c2 = 1.928 ,     &
                           & c3 = 0.0858 , c4 = -0.0617 , c5 = 2.5 ,    &
                           & y1 = 1331.647 , y2 = 1.561228 ,            &
                           & y3 = 0.38194
!
      real(8) :: dum , rhair , rhop
      real(8) :: ustart01
      intent (in) dum , rhair , rhop
!
      real(8) :: dm , rep , term , term1 , term2
! 
      dm = dum  !* 1.0e-4      ! cm
      rep = y1*(dm**y2) + y3
      term1 = sqrt(1.0+(c1/(rhop*gti*0.1*(dm**c5))))
      term2 = sqrt(rhop*gti*100.0*dm/rhair)
      term = term1*term2
      ustart01 = cvmgt(a2*term*(1.0-c3*exp(c4*(rep-10.0))),             &
               & a2*term/sqrt(c2*(rep**0.092)-1.0),rep.gt.10.0)
! 
      end function ustart01
!
!     *****************************************************************
!     *                                                            ****
!     * modified by a.s.zakey, nov.2003                            ****
!     * y. shao, 13 june 2000                                      ****
!     * calculate ustar0(d) using shao and lu (2000) for uncovered ****
!     * dry surface                                                ****
!     * dum:    particle diameter                   [um]           ****
!     * ustar0: threshold friction velocity       [cm/s]           ****
!     *****************************************************************
! 
      function ustart0(rhop,dum,rhoa)
      implicit none
!
      real(8) , parameter :: agamma = 3.0E-4 , f = 0.0123
!
      real(8) :: dum , rhoa , rhop
      real(8) :: ustart0
      intent (in) dum , rhoa , rhop
!
      real(8) :: dm , sigma
! 
      sigma = rhop/rhoa
      dm = dum*1.0E-2
      ustart0 = f*(sigma*gti*dm+agamma/(rhoa*dm))
      ustart0 = sqrt(ustart0)
      ustart0 = ustart0*100.0
!
      end function ustart0
! 
!  **********************************************************
!  *  dust emission scheme                             ******
!  *                                                   ******
!  * this scheme based on marticorena and bergametti,  ******
!  * 1995; gong et al.,(2003); alfaro et al.,(1997)    ******
!  *                                                   ******
!  * the modification coded by:                        ******
!  * ashraf s. zakey                                   ******
!  **********************************************************
! 
      subroutine sfflux(ilg,il1,il2,jloop,luc,ivegcov,vegfrac, &
                      & ustarnd,z0,soilw,surfwd,roarow,trsize,rsfrow)
! 
      implicit none
!
      integer :: il1 , il2 , ilg , jloop , luc
      integer , dimension(ilg) ::  ivegcov
      real(8) , dimension(ilg) :: roarow , soilw , surfwd , vegfrac ,   &
                       &          z0 , ustarnd
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(nbin,2) :: trsize
      intent (in) il1 , il2 ,  ivegcov , jloop , roarow ,     &
                & soilw , surfwd , vegfrac , z0 , ustarnd
      intent (out) rsfrow
!
      integer :: i , ieff , ieffmax , n , ns
      real(8) , dimension(ilg) :: xclayrow , xroarow , xsoilw ,         &
                         & xsurfwd , xvegfrac , xz0 , xustarnd
      real(8) , dimension(ilg,20) :: xfland
      real(8) , dimension(ilg,nbin) :: xrsfrow
      real(8) , dimension(ilg,nats) :: xsand2row,xftex
      real(8) , dimension(ilg,nsoil,nats) :: xsrel2d
! 
      rsfrow = 0.0
!     effective emitter cell ( depending on ivegcov)
      xvegfrac = 0.
      xftex = 0.
      xsoilw = 0.
      xsurfwd = 0.
      xz0 = 0.
      xclayrow = 0.
      xroarow = 0.
      xsrel2d = 0.
      xsand2row = 0.
      xustarnd=0.
      xfland = 0.0
      xrsfrow = 0.
 
      ieff = 0
      ieffmax = 0
      do i = il1 , il2
        if (ivegcov(i).eq.8 .or. ivegcov(i).eq.11) then   
          ieff = ieff + 1
          xvegfrac(ieff) = vegfrac(i)
          xsoilw(ieff) = soilw(i)
          xsurfwd(ieff) = surfwd(i)
          xz0(ieff) = z0(i)
          xroarow(ieff) = roarow(i)
          xustarnd(ieff) = ustarnd(i) 
          xclayrow(ieff) = clayrow2(i,jloop)
          do n = 1 , nats
            xftex(ieff,n) = dustsotex(i,jloop,n)
            xsand2row(ieff,n) = sand2row2(i,n,jloop)
          do  ns = 1 , nsoil
             xsrel2d(ieff,ns,n) = srel2d(i,jloop,ns,n)
          end do
          end do
          do n = 1 ,luc 
            xfland(ieff,n) = 1
          end do
        end if
      end do
 
      ieffmax = ieff

!     if ( ieffmax>0. ) print *, maxval(xustarnd)

      if ( ieffmax.gt.0 ) call dust_module(1,ieffmax,ilg,trsize,xsoilw, &
         & xvegfrac,xsurfwd,xftex,xclayrow,xroarow,xz0,xsrel2d,         &
         & xustarnd,xrsfrow)
        
!     if ( ieffmax>0. ) print *, 'FLUX : ' , maxval(xrsfrow)
!     put back the dust flux on the right grid
 
      ieff = 0
 
      do i = il1 , il2
        if  (ivegcov(i).eq.8 .or. ivegcov(i).eq.11) then     
          ieff = ieff + 1
          do n = 1 , nbin
            rsfrow(i,n) = xrsfrow(ieff,n)
          end do
        end if
      end do
 
      end subroutine sfflux
! 
      subroutine dust_module(il1,il2,ilg,trsize,soilw,vegfrac,surfwd,ftex,&
                       & clayrow,roarow,z0,srel,ustarnd,rsfrow)
 
      implicit none
!
      integer :: il1 , il2 , ilg
      real(8) , dimension(ilg) :: clayrow , roarow , soilw , surfwd ,   &
                                & vegfrac , z0 , ustarnd
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(ilg,nats) :: ftex
      real(8) , dimension(ilg,nsoil,nats) :: srel
      real(8) , dimension(nbin,2) :: trsize
      intent (in) clayrow , soilw , surfwd , z0 , ustarnd , ftex
!
      real(8) , dimension(ilg) :: alamda , hc , rc , srl , wprim
      real(8) :: arc1 , arc2 , br , cly1 , cly2 , sigr , tempd ,        &
               & umin , ustarns , uth , utmin , x , xz , ym , z0s
      integer :: i
      real(8) , dimension(ilg) :: ustar
      real(8) , dimension(ilg,nsoil) :: utheff
!
      data umin/15./
      data xz/0.25/ , br/202.0/ , ym/0.16/ , sigr/1.45/
      data z0s/3.E-3/ , x/10./

      do i = il1 , il2
 
        srl(i) = z0(i)*100.0
        rc(i) = 1.0
 
        if ( jfs.eq.0 ) then
 
!         * raupach et al. (1993)                                     

          if ( vegfrac(i).lt.1.0 ) then
            alamda(i) = xz*(log(1.0-vegfrac(i)))*(-1.0)
            arc1 = sigr*ym*alamda(i)
            arc2 = br*ym*alamda(i)
            if ( arc1.le.1.0 .and. arc2.le.1.0 ) rc(i)                  &
               & = (sqrt(1.0-arc1)*sqrt(1.0+arc2))
          end if
 
        else if ( jfs.eq.1 ) then
!
!         Marticorena et al., 1997: correction factor for non erodible elements
!  
          rc(i) = 1 - (dlog(0.5E-2/z0s)/(dlog(0.35*(x/z0s)**0.8)))
 
        end if
 
!       threshold velocity correction for soil humidity hc
 
        if ( jsoilm.eq.0 ) then
 
          if ( soilw(i).lt.0.0 ) then
            write (aline,*) 'hc, rc =' , soilw(i) , ' less than zero'
            call say
            call fatal(__FILE__,__LINE__,'NEGATIVE SOILW')
          else if ( soilw(i).lt.0.03 ) then
            hc(i) = exp(22.7*soilw(i))
          else if ( soilw(i).ge.0.03 ) then
            hc(i) = exp(95.3*soilw(i)-2.029)
          else
            hc(i) = 1.0
          end if
 
        else if ( jsoilm.eq.1 ) then
 
          cly1 = clayrow(i)
          cly2 = cly1*cly1
        
          wprim(i) = 0.0014*cly2 + 0.17*cly1
          tempd=  dmax1(0.00001d0,soilw(i)*100.0 -wprim(i))
!          print*,'humidity',i,cly1,soilw(i)*100,wprim(i),tempd
          if ( soilw(i)*100.gt.wprim(i) ) then
            hc(i) = sqrt(1.0+1.21*tempd**0.68)
!          print*,'hc',i,hc(i)
          else
            hc(i) = 1.0
          end if
 
!       no soil humidity correction facor if jsoilm > 1
        else
          hc(i)=1.0
        end if
 
!       * total correction factor for both hc and rc                
        rc(i) = rc(i)/hc(i)
 
!       *     computation of the wind friction velocity              
!       ***** *     accounting for the increase of the roughness length
!       ***** *     due to the saltation layer (gillette etal. jgr 103,
!       ***** *     no. d6, p6203-6209, 1998                           
!       *****
!       ustarns = (vonkar*100.*surfwd(i))/(log(1000./srl(i)))

        ustarns = ustarnd(i)*100 !cm.s-1
        utmin = (umin/(100.*vonkar*rc(i)))*log(1000./srl(i))
 
        if ( surfwd(i).ge.utmin ) then
          ustar(i) = ustarns + 0.3*(surfwd(i)-utmin)*(surfwd(i)-utmin)
        else
          ustar(i) = ustarns
        end if
 
      end do       ! end i loop
 
      call uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop)
 
      call emission(ilg,il1,il2,rhop,ftex,uth,roarow,rc,utheff,     &
                  & ustar,srel,rsfrow,trsize,vegfrac)
 
      end subroutine dust_module
! 
      subroutine uthefft(il1,il2,ilg,ust,nsoil,roarow,utheff,rhop)
      implicit none
!
      integer :: il1 , il2 , ilg , nsoil , ust
      real(8) :: rhop
      real(8) , dimension(ilg) :: roarow
      real(8) , dimension(ilg,nsoil) :: utheff
      intent (in) il1 , il2 , ilg , nsoil , ust
      intent (out) utheff
!
      integer :: i , j
!
      do i = 1 , nsoil
        do j = il1 , il2
          if ( ust.eq.0 ) utheff(j,i) = ustart0(rhop,dp_array(i),roarow(j))
          if ( ust.eq.1 ) utheff(j,i) = ustart01(rhop,dp_array(i),roarow(j))
        end do
      end do
 
      end subroutine uthefft
!
      subroutine emission(ilg,il1,il2,rhop,ftex,uth,roarow,rc,      &
                        & utheff,ustar,srel,rsfrow,trsize,vegfrac)
 
      implicit none
!
      integer :: il1 , il2 , ilg
      real(8) :: rhop , uth
      real(8) , dimension(ilg) :: rc ,ustar, roarow , vegfrac
      real(8) , dimension(ilg,nbin) :: rsfrow
      real(8) , dimension(ilg,nats) :: ftex
      real(8) , dimension(ilg,nsoil,nats) :: srel
      real(8) , dimension(nbin,2) :: trsize
      real(8) , dimension(ilg,nsoil) :: utheff
      intent (in)  il1 , il2 , ilg , rc , rhop , roarow , srel ,  &
                & trsize , ustar , utheff , vegfrac, ftex
      intent (inout) rsfrow , uth
!
!     real(8) :: beffect
      real(8) :: beta , const,                                & 
               & p1, p2 , p3, rwi, dec, ec , fdp1 , fdp2
      real(8) , dimension(ilg,nats) :: fsoil , fsoil1 , fsoil2 , fsoil3
      integer :: i , k , n , nt , ns
      real(8) , dimension(ilg,isize,nats) :: rsfrowsub
      real(8), dimension(ilg,nbin,nats):: rsfrowt
!
!     Put const consistent with soil parameters and Laurent et al., 08
      data const/1./, beta/16300./ 
 
      p1 = 0.0
      p2 = 0.0
      p3 = 0.0
      fsoil(:,:) = 0.
      fsoil1(:,:) = 0.
      fsoil2(:,:) = 0.
      fsoil3(:,:) = 0.
 
      do nt = 1 , nats
        do ns = 1 , nsoil
          do i = il1 , il2
 
            if ( rc(i).gt.0.0 .and. ustar(i).ne.0. ) then
              uth = utheff(i,ns)/(rc(i)*ustar(i))
 
              if ( uth.le.1.0 ) then
 
                fdp1 = ustar(i)**3*(1.0-uth*uth)
                fdp2 = (1.0+uth)*const*(1.E-5)*roarow(i)*rgti
 
                if ( fdp2.le.0.0 ) fdp2 = 0.
 
! FAB: with subgrid soil texture, the aggregation of vertical fluxes per texture type
! at the grid cell level is done in fine.  
!               fsoil(k) = srel(k,j,i)*fdp1*fdp2*aeffect*beffect
! FAB 
                fsoil(i,nt) = srel(i,ns,nt)*fdp1*fdp2 
 
!               size-distributed kinetic energy flux(per texture type)
                dec = fsoil(i,nt)*beta
!               individual kinetic energy for an aggregate of size dp (
!               g cm2 s-2) cf alfaro (dp) is in cm
                ec = (mathpi/12)*rhop*1E-3*(dp_array(ns)**3.0)*                &
                    & (20*ustar(i))**2.0
 
                if ( ec.gt.e1 ) then
                  p1 = (ec-e1)/(ec-e3)
                  p2 = (1-p1)*(ec-e2)/(ec-e3)
                  p3 = 1 - p1 - p2
                else if ( ec.gt.e2 .and. ec.le.e1 ) then
                  p1 = 0.
                  p2 = (ec-e2)/(ec-e3)
                  p3 = 1 - p2
                else if ( ec.gt.e3 .and. ec.le.e2 ) then
                  p1 = 0.
                  p2 = 0.
                  p3 = 1.
                else if ( ec.le.e3 ) then
                  p1 = 0.
                  p2 = 0.
                  p3 = 0.
                end if
 
                fsoil1(i,nt) = fsoil1(i,nt) + 1.E-2*p1*(dec/e1)* &
                               (mathpi/6.)*rhop*((d1*1.E-04)**3.)
                fsoil2(i,nt) = fsoil2(i,nt) + 1.E-2*p2*(dec/e2)* &
                               (mathpi/6.)*rhop*((d2*1.E-04)**3.)
                fsoil3(i,nt) = fsoil3(i,nt) + 1.E-2*p3*(dec/e3)* &
                               (mathpi/6.)*rhop*((d3*1.E-04)**3.)
              end if
            end if
          end do
        end do
      end do
!
! calculate fluxes for each of transport bins
!
      rsfrowt(:,:,:) = 0.
      do nt = 1 , nats
        do n = 1 , isize
          do i = il1 , il2
!         discretisation of the modal emission in isize emission sub bin
            rsfrowsub(i,n,nt) = fsoil1(i,nt)*frac1(n) + &
                                fsoil2(i,nt)*frac2(n) + &
                                fsoil3(i,nt)*frac3(n)
!         and in tranport bins (nbin)
            rwi = (aerosize(1,n)+aerosize(2,n))/2.0*1.E6

            do k = 1 , nbin
              if ( rwi.ge.trsize(k,1) .and. rwi.lt.trsize(k,2) )          &
                 & rsfrowt(i,k,nt) = rsfrowt(i,k,nt) + rsfrowsub(i,n,nt)
            end do
          end do
        end do
      end do 

! Finally, aggregation of the dust flux at the grid cell level =
! weighted sum over soil texture
! weighting by grid cell veg fraction, snow fraction, 
! EBL = erodibility factor, to be introduced ) 

!       f = 0.0D0
!       aeffect = (1-f)*(1-vegfrac(k))
!       beffect = 0.01*fland(k,i)*sand2row(k,i)
! Fab : fland is equal to 1 with bats
! the fraction of sand ( coarse particles) is intrinsically contained in dsrel soil aggregate distribution
! there is no need to multipky by sand2row.

       rsfrow(:,:) = 0.
       do k = 1 , nbin
         do nt = 1 , nats
           do i = il1 , il2
              rsfrow(i,k) =  rsfrow(i,k) + rsfrowt(i,k,nt)*ftex(i,nt)  &
                          &  * (1 - vegfrac(i)) 
! * EBL(i)
!                         &  * (1-snowfrac)     
           end do
         end do
       end do
!
      end subroutine emission
!
      end module mod_dust
