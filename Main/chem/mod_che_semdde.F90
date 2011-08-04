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
 
module mod_che_semdde
!
! Chemical and aerosol surface emission and dry deposition
!
  use m_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_main
  use mod_che_trac
  use mod_che_dust
  use mod_che_aerosol
  use mod_message
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_indices
!
  private
!
  public :: chsrfem , chdrydep
!
! Dynamic Viscosity Parameters
!
  real(dp) , parameter :: a1 = 145.8D0
  real(dp) , parameter :: a2 = 1.5D0
  real(dp) , parameter :: a3 = 110.4D0
!
! Molecular Free Path calculation parameters
!
  real(dp) , parameter :: c1 = 6.54D-8
  real(dp) , parameter :: c2 = 1.818D-5
  real(dp) , parameter :: c3 = 1.013D5
  real(dp) , parameter :: c4 = 293.15D0
!
! Cunningham slip correction factor parameters
!
  real(dp) , parameter :: aa1 = 1.257D0
  real(dp) , parameter :: aa2 = 0.4D0
  real(dp) , parameter :: aa3 = 1.1D0
!
  contains
!
! SURFACE EMISSION
!
  subroutine chsrfem

    use netcdf
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
!
    integer :: i , j , k , l , m , n
    integer , parameter :: iutopt = 12
    logical :: there
    integer :: itr , ierr
!
!   fisrt activate dust initialization
!
    write (aline, *) 'Calling inidust'
    call say(myid)
    call inidust

!   read the monthly aerosol emission files

    if ( myid == 0 ) then
      chemsrc_io = d_zero
      if (aertyp(4:5) /= '00') then
        call read_aerosol(chtrname,chemsrc_io)
      end if
      do j = 1 , jx
        do itr = 1 , ntr
          do m = 1 , mpy
            do i = 1 , iy
              src_0(i,m,itr,j) = chemsrc_io(i,j,m,itr)
            end do
          end do
        end do
      end do
    end if
!
    call mpi_scatter(src_0,iy*mpy*ntr*jxp,mpi_real8,         &
                     src0, iy*mpy*ntr*jxp,mpi_real8,         &
                     0,mpi_comm_world,ierr)
    do j = 1 , jendl
      do itr = 1 , ntr
        do m = 1 , mpy
          do i = 1 , iy
            chemsrc(i,j,m,itr) = src0(i,m,itr,j)
          end do
        end do
      end do
    end do

!   sulfates sources

    do m = 1 , mpy
      do j = 1 , jendl
        do i = 1 , iy
          if ( iso4 > 0 ) chemsrc(i,j,m,iso4)                        &
               = 0.02D0*chemsrc(i,j,m,iso2)
          if ( iso2 > 0 ) chemsrc(i,j,m,iso2)                        &
               = 0.98D0*chemsrc(i,j,m,iso2)
   
!         partition hydrophilic hydrophonic ( cooke et al.1999)
!         BC
          if ( ibchb > 0 .and. ibchl > 0 ) then
            chemsrc(i,j,m,ibchl) = 0.2D0*chemsrc(i,j,m,ibchb)
            chemsrc(i,j,m,ibchb) = 0.8D0*chemsrc(i,j,m,ibchb)
          end if
!         OC
          if ( iochb > 0 .and. iochl > 0 ) then
            chemsrc(i,j,m,iochl) = chemsrc(i,j,m,iochb)*d_half
            chemsrc(i,j,m,iochb) = chemsrc(i,j,m,iochb)*d_half
          end if
        end do
      end do
    end do
   
!   Optical properties / internal mixing

    if ( mixtype == 2 ) then
      if ( myid == 0 ) then
        inquire (file='optdat.bin',exist=there)
        if ( .not.there ) then
          write (*,*) 'For mixtype=2, optdat.bin is required'
          write (*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
          call fatal(__FILE__,__LINE__,'optdat.bin is required')
        end if
        open (iutopt,file='optdat.bin',form='unformatted',            &
              recl=4*19*11*11*11*11*ibyte,access='direct')
        read (iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            k=1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            k=1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                            1,11),l=1,11),m=1,11),n=1,11)
        close (iutopt)
      end if
      call mpi_bcast(dextmix,4*19*11*11*11*11,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(dssamix,4*19*11*11*11*11,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(dgmix,4*19*11*11*11*11,mpi_real,0,mpi_comm_world,ierr)

!     Check !
   
      do k = 1 , 11
        do l = 1 , 11
          do m = 1 , 11
            do n = 1 , 11
   
              if ( k+l+m+n == 14 ) then
                do i = 1 , 4
                  do j = 1 , 19
   
                    if ( (dextmix(i,j,k,l,m,n) < d_zero) .or.   &
                         (dextmix(i,j,k,l,m,n) > 20.0D0) ) then
                      write (aline,*) 'problem in dextmix ' ,   &
                                      dextmix(i,j,k,l,m,n)
                      call say(myid)
                      call fatal(__FILE__,__LINE__,'DETMIX ERROR')
                    end if
   
                    if ( (dssamix(i,j,k,l,m,n) < d_zero) .or.   &
                         (dssamix(i,j,k,l,m,n) > d_one) ) then
                      write (aline,*) 'problem in dssamix ' ,   &
                                      dssamix(i,j,k,l,m,n)
                      call say(myid)
                      call fatal(__FILE__,__LINE__,'DSSAMIX ERROR')
                    end if
   
                    if ( (dgmix(i,j,k,l,m,n) < d_zero) .or.     &
                         (dgmix(i,j,k,l,m,n) > d_one) ) then
                      write (aline,*) 'problem in dgmix ' ,     &
                                      dgmix(i,j,k,l,m,n)
                      call say(myid)
                      call fatal(__FILE__,__LINE__,'DGMIX ERROR')
                    end if
   
                  end do
                end do
              end if
   
            end do
          end do
        end do
      end do
   
      if ( myid == 0 ) write (*,*) '! OPDATA CHECKED !'

    end if

  end subroutine chsrfem
!
!**************************************************************
!  dry depostion scheme
!  this scheme based on :
! - zhang et al,(2001) : a size-segregated particle
!   dry deposition scheme for an atmospheric aerosol
!   module, atmos. env. 35, 549-560
!
! - giorgi, f. (1986): a particle dry deposition
!      parameterization scheme for use in tracer
!      transport models. jgr,91, 9794-9806
!
!  input:
!  =====
! - throw : temperature in k
! - roarow : air density
! - shj    : local mid-layer sigma value
! - pressg : grid row of surface pressure [pa]
! - temp2  : temperature at 10m. (deg k)
! - sutemp : surface temperature (deg k)
! - srad   : solar irradiance at the ground(w/m**2)
! - rh10   : relative humidity of air at 10m.
! -vegcover: vegetation cover
! - wind10 : wind at 10 m
! - xrow   : dust concentration [kg/kg]
! - ustar  : u*
! - ilev   : number of model level
!
! output
! ======
! - rtdry  : dry deposition tendency
!
!**************************************************************
!
  subroutine chdrydep(ilg,il1,il2,ilev,luc,nbin,ivegcov,throw,      &
                      roarow,shj,pressg,temp2,sutemp,srad,rh10,     &
                      wind10,zeff,trsize,pdepv)
!
    implicit none
!
    integer :: il1 , il2 , ilev , ilg , luc , nbin
    integer , dimension(ilg) :: ivegcov
    real(dp) , dimension(ilg,ilev,nbin) :: pdepv
    real(dp) , dimension(ilg) :: pressg , rh10 , srad , sutemp ,       &
                                temp2 , wind10 , zeff
    real(dp) , dimension(ilg,ilev) :: roarow , throw
    real(dp) , dimension(ilev) :: shj
    real(dp) , dimension(nbin,2) :: trsize
    intent (in) il1 , il2 , ilev , ilg , ivegcov , luc , nbin ,       &
                pressg , rh10 , roarow , shj , srad , sutemp , temp2 ,&
                throw , trsize , wind10 , zeff
    intent (inout) pdepv
!
    real(dp) :: amfp , amob , asq , ch , cm , cun , dtemp , dthv , eb ,&
               eim , ein , es , fh , fm , frx1 , kui , logratio ,     &
               mol , pre , prii , priiv , psit , psiu , ptemp2 , qs , &
               r1 , ratioz , aa , rib , st , tbar , thstar , tsv ,    &
               tsw , ustarsq , utstar , vp , vptemp , wvpm , ww , x , &
               y , z , z0water , zdl , zl
    real(dp) , dimension(ilg,ilev) :: amu
    real(dp) , dimension(ilg) :: anu , schm , zz0
    real(dp) , dimension(ilg,ilev,isize) :: cfac , pdepvsub , pdiff ,  &
           rhsize , taurel
    integer :: i , j , jc , k , kcov , l , lev , n , tot
    real(dp) , dimension(ilg,luc) :: ra , ustar , vegcover
    real(dp) , dimension(ilg,luc,isize) :: rs
   
    real(dp) , parameter :: z10 = d_10
    real(dp) , dimension(isize) :: avesize
!
    i = 0
    pdepvsub = d_zero

    do n = 1 , isize
      avesize(n) = (aerosize(1,n)+aerosize(2,n))*d_half
    end do
   
!======================================================================
!   ********************************************************
!   *   aerosize - dry radius                           ****
!   *   rhop  - density for each aerosol type           ****
!   ********************************************************
    do n = 1 , isize
      do l = 1 , ilev
        do i = il1 , il2
   
!         ********************************************************
!         *  aerosol gravitational settling velocity          ****
!         *  and diffusion coefficient                        ****
!         *                                                   ****
!         * air's dynamic viscosity                           ****
!         ********************************************************
   
          amu(i,l) = a1*1.0D-8*throw(i,l)**a2/(throw(i,l)+a3)
   
!         . . . . mid layer pressure in [pascal].
          pre = pressg(i)*shj(l)
!         ********************************************************
!         * mean molecular free path.                         ****
!         *     k.v. beard [1976], j atm. sci., 33            ****
!         ********************************************************
   
          amfp = c1*(amu(i,l)/c2)*(c3/pre)*(throw(i,l)/c4)**(d_half)
          prii = d_two/9.0D0*egrav/amu(i,l)
          priiv = prii*(rhop-roarow(i,l))
   
!         ********************************************************
!         * cunningham slip correction factor and             ****
!         * relaxation time = vg/grav.                        ****
!         ********************************************************
   
          cfac(i,l,n) = d_one + amfp/avesize(n) &
                        *(aa1+aa2*dexp(-aa3*avesize(n)/amfp))
          taurel(i,l,n) = dmax1(priiv*avesize(n)**d_two* &
                                cfac(i,l,n)*regrav,d_zero)
   
!         ********************************************************
!         * stokes friction                                  *****
!         pdepvsub(i,l,n) ' sellting dep. velocity = '
!         ********************************************************
   
          pdepvsub(i,l,n) = taurel(i,l,n)*egrav
        end do
      end do
    end do
!======================================================================
   
!======================================================================
!   ****************************************************
!   * ra : is the aerodynamic resistance above the  ****
!   *      canopy and it is function in u* and      ****
!   *      z0: roughness length and the stability   ****
!   *      function                                 ****
!   * mol  - monin obukhov length (m) - calculated  ****
!   *           for each land use category          ****
!   * ptemp2 -potential temperature at z2  (deg. k) ****
!   * temp2 - temperature at 10m. (deg k)           ****
!   * z10   - 10 m.                                 ****
!   * sutemp -surface temperature (deg k)           ****
!   * srad   -solar irradiance at the ground(w/m**2)****
!   * rh10  - relative humidity of air at 10m.      ****
!   *           (0.0-1.0)                           ****
!   * stdpmb - sea level pressure (mb)              ****
!   ****************************************************
    do j = 1 , luc
      do i = il1 , il2
        ww = dmax1(wind10(i),d_one)
        zz0(i) = zeff(i)
   
!       *********************************************************
!       *  potential temperature at z2  (deg. k)                 
!       *****
   
!       *********************************************************
        ptemp2 = temp2(i) + z10*0.0098D0
   
!       *********************************************************
!       * for calculations over water compute values of critical 
!       ***** * profile variables: l and ustar                       
!       ***** * begin for water                   
!       *****
   
!       *********************************************************
        if ( ivegcov(i) == 0 ) then
   
!       ********************************************************
!       *  vp  - vapour pressure at z2                          
!       ***** *  wvpm- water vapour mixing ratio at  z2            
!       ***** *  vptemp- virtual potential temperature at z2 (deg. k)
!       ********************************************************
          es = 6.108D0*dexp(17.27D0*(temp2(i)-tzero)/(temp2(i)-35.86D0))
          vp = rh10(i)*es
          wvpm = ep2*vp/(stdpmb-vp)
          vptemp = ptemp2*(d_one+0.61D0*wvpm)
   
!         ********************************************************
!         *  assume rh10 at water surface is 100%                 
!         ***** *   vp = es(tsw-tzero) !sat. vap press at surface   
!         ***** *   saturated vapour pressure at surface             
!         ***** *   saturated mixing ratio at surface                
!         ***** *   tsv - virtual potential temperature at surface   
!         ***** *           (deg. k)                                 
!         *****
   
!         ********************************************************
          tsw = sutemp(i)
          vp = 6.108D0*dexp(17.27D0*(tsw-tzero)/(tsw-35.86D0))
          qs = ep2*vp/(stdpmb-vp)
          tsv = tsw*(d_one+0.61D0*qs)
          z0water = 1.0D-4
!         ******************************************************
!         * scalet  :  not required if  z2 = 10m                  
!         *****
!         ******************************************************
          dthv = (vptemp-tsv)
!         ******************************************************
!         * calculate drag coefficient cun with neutral condition 
!         ***** * assumption  garratt (1977)                         
!         *****
   
!         ******************************************************
          cun = 7.5D-4 + 6.7D-5*ww
          mol = 9999.0D0
   
          if ( dabs(dthv) > 1.0D-6 ) then
            mol = vptemp*cun**1.5D0*ww**d_two/(5.096D-3*dthv)
          end if
          if ( mol > d_zero .and. mol < d_five ) mol =  d_five
          if ( mol > -d_five .and. mol < d_zero ) mol = -d_five
          zdl = z10/mol
!
          if ( zdl < d_zero ) then
   
!           **************************************************
!           * wind speed                     
!           **************************************************
            x = (d_one-15.0D0*zdl)**d_rfour
            psiu = d_two*dlog((d_one+x)*d_half) +   & 
                         dlog((d_one+x*x)*d_half) - &
                   d_two*datan(x) + halfpi
   
!           ****************************************************
!           * pot temp                        
!           ****************************************************
            y = dsqrt(d_one-9.0D0*zdl)
            psit = d_two*0.74D0*dlog((d_one+y)*d_half)
          else
            psiu = -4.7D0*zdl
            psit = psiu
          end if
          z0water = 0.000002D0*ww**2.5D0
!
          ustar(i,j) = vonkar*ww/(dlog(z10/z0water)-psiu)
          thstar = vonkar*(ptemp2-sutemp(i))/(0.74D0*dlog(z10/z0water)-psit)
          zz0(i) = z0water
!
        else
   
!         ******************************************************
!         * compute ustar and l for land use categories other than 
!         **** * water use louis method. !pkk 7/16/85, find bulk     
!         **** * richardson number.                                  
!         ****
   
!         ******************************************************
          rib = egrav*z10*(ptemp2-sutemp(i))/(sutemp(i)*ww**d_two)
   
!         *******************************************************
!         * ensure that conditions over land are never stable when 
!         ***** * there is incoming solar radiatiom                  
!         *****
!         *******************************************************
          if ( srad(i) > d_zero .and. rib > d_zero ) rib = 1.0D-15
!
          dtemp = ptemp2 - sutemp(i)
          if ( dabs(dtemp) < 1.0D-10 ) dtemp = dsign(1.0D-10,dtemp)
          tbar = (ptemp2+sutemp(i))*d_half
!
          ratioz = z10/zz0(i)
          logratio = dlog(ratioz)
          asq = 0.16D0/(logratio**d_two)
!
          if ( rib <= d_zero ) then
            aa = asq*9.4D0*dsqrt(ratioz)
            cm = 7.4D0*aa
            ch = 5.3D0*aa
            fm = d_one - (9.4D0*rib/(d_one+cm*dsqrt(dabs(rib))))
            fh = d_one - (9.4D0*rib/(d_one+ch*dsqrt(dabs(rib))))
          else
            fm = d_one/((d_one+4.7D0*rib)**d_two)
            fh = fm
          end if
!
          ustarsq = asq*ww**d_two*fm
          utstar = asq*ww*dtemp*fh/0.74D0
          ustar(i,j) = dsqrt(ustarsq)
          thstar = utstar/ustar(i,j)
!
          mol = tbar*ustarsq/(vonkar*egrav*thstar)
        end if
   
        kui = d_one/(vonkar*ustar(i,j))
   
!       **************************************************************
!       * compute the values of  ra                            *******
!       **************************************************************
   
        z = z10
        zl = z/mol
   
        if ( zl >= d_zero ) then
          ra(i,j) = kui*(0.74D0*dlog(z/zz0(i))+4.7D0*zl)
        else
          ra(i,j) = kui*0.74D0*(dlog(z/zz0(i))- &
                    d_two*dlog((d_one+dsqrt(d_one-9.0D0*zl))*d_half))
        end if
        ra(i,j) = dmax1(ra(i,j),0.99D0)
        ra(i,j) = dmin1(ra(i,j),999.9D0)
      end do
    end do
!    
!======================================================================
!   
!   *****************************************************
!   * the schmidt number is the ratio of the         ****
!   * kinematic viscosity of air to the particle     ****
!   * brownian diffusivity ===> sc=v/d               ****
!   *****************************************************
   
    do n = 1 , isize
      do l = 1 , ilev
        do i = il1 , il2
   
!         *****************************************************
!         * for now we will not consider the humidity      ****
!         * impact so we will set the variable frx1=1.0    ****
!         * i.e. only dry particles                        ****
!         *****************************************************
   
          frx1 = d_one
          rhsize(i,l,n) = avesize(n)*frx1
          anu(i) = amu(i,l)/roarow(i,l)
          amob = 6.0D0*mathpi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
          pdiff(i,l,n) = boltzk*throw(i,l)/amob
          schm(i) = anu(i)/pdiff(i,l,n)
   
!         ******************************************************
!         * for brownian diffusion, there is evidence that  ****
!         * its fromula depend on schmidt number as :       ****
!         * eb= schm x c^gama                               ****
!         * where gama is efficiency factor and its value   ****
!         * between 1/2 and 2/3 with larger values for      ****
!         * rougher surfaces                                ****
!         * ****************************************************
   
        end do
        if ( l == ilev ) then
          do k = 1 , luc ! luc  = 1 for the moment
            do i = il1 , il2
!
!======================================================================
!
!             find the right table index for the cell cover ( ocean and lake
!             are 0 in the ivegcov and 14-15 in the table )

              if ( ivegcov(i) == 0 ) then
                kcov = 14
              else
                kcov = ivegcov(i)
              end if

   
!             ******************************************************
!             * the parameter governing impaction processes is *****
!             * the stokes number,st, which has the form of    *****
!             * 1) st = vg x u* /g a for vegetated surefaces   *****
!             *    (slinn, 1982)                               *****
!             * 2) st = vg x u*2/anu for smothed surfaces or   *****
!             *    surfaces with bluff roughness elements      *****
!             *    (giorgi,1988)                               *****
!             ******************************************************
   
              st = taurel(i,l,n)*ustar(i,k)*ustar(i,k)/anu(i)
   
              eb = schm(i)**(-0.666667D0)
!             eim=(st/(st+aest(k)))**2
              eim = (st/(st+aest(kcov)))**d_two
   
   
              eim = dmin1(eim,0.6D0)
              ein = d_zero
!             if (arye(k) > 0.0001) then
!               ein = (1000.0*2.0*avesize(n)/arye(k))**1.5
!             end if
   
              if ( arye(kcov) > 0.0001D0 ) then
                ein = (d_1000*d_two*avesize(n)/arye(kcov))**1.5D0
              end if
   
              ein = dmin1(ein,d_half)
   
!             *****************************************************
!             * partickes larger than 5 micro may rebounded   *****
!             * after hitting a surface,this process may be   *****
!             * included by modifying the total collection    *****
!             * by the factor of r1, which represents the     *****
!             * fraction of particles sticking to the surface *****
!             * slinn (1982) suggested the following:         *****
!             * r = exp (- st^0.2)                            *****
!             *****************************************************
   
!             r1 = exp (-st**d_half
              r1 = dmax1(d_half,dexp(-st**d_half))
!             if (k >= 11 .and. r1 < d_half r1=d_half
              if ( kcov >= 11 .and. r1 < d_half ) r1 = d_half
              if ( r1 < 0.4D0 ) r1 = 0.4D0
   
!             ***************************************************
!             * calculation of rs: the surface resistance   *****
!             * which depends on the collection efficiency  *****
!             * of the surface and is determined by the     *****
!             * various deposition processes                *****
!             ***************************************************
   
!             rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
              rs(i,k,n) = d_one/d_three/ustar(i,k)/(eb+eim+ein)/r1
            end do
          end do
        end if
      end do
    end do
!
!======================================================================
!
!   note for the moment vegcover = 1 ( one type of cover ivegcov per
!   gride cell)
    do n = 1 , isize
      do jc = 1 , luc
        do i = il1 , il2
!         modify the surface layer
          vegcover(i,jc) = d_one
          pdepvsub(i,ilev,n) = pdepvsub(i,ilev,n) + vegcover(i,jc) *  &
                               d_one/(ra(i,jc)+rs(i,jc,n))
   
        end do
      end do
   
    end do  ! end nsize loop
! 
!   average deposition velocities on bin
! 
    do k = 1 , nbin
      tot = 0
      do lev = 1 , ilev
        do i = 1 , ilg
          pdepv(i,lev,k) = d_zero
        end do
      end do
      do n = 1 , isize
        if ( avesize(n)*1.0D6 >= trsize(k,1) .and. &
             avesize(n)*1.0D6 <  trsize(k,2) ) then
          do lev = 1 , ilev
            do i = 1 , ilg
              pdepv(i,lev,k) = pdepv(i,lev,k) + pdepvsub(i,lev,n)
            end do
          end do
          tot = tot + 1
        end if
      end do
      if ( tot > 0 ) then
        do lev = 1 , ilev
          do i = 1 , ilg
            pdepv(i,lev,k) = pdepv(i,lev,k)/dble(tot)
          end do
        end do
      end if
    end do
! 
  end subroutine chdrydep
!
end module mod_che_semdde
