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
 
      subroutine chdrydep(ilg,il1,il2,ilev,luc,nbin,ivegcov,throw,      &
                        & roarow,shj,pressg,temp2,sutemp,srad,rh10,     &
                        & wind10,zeff,trsize,pdepv)

!**************************************************************
!  dry depostion scheme for dust particles                *****
!  this scheme based on :                                 *****
! - zhang et al,(2001) : a size-segregated particle       *****
!   dry deposition scheme for an atmospheric aerosol      *****
!   module, atmos. env. 35, 549-560                       *****
!                                                         *****
! - giorgi, f. (1986): a particle dry deposition          *****
!      parameterization scheme for use in tracer          *****
!      transport models. jgr,91, 9794-9806                *****
!                                                         *****
!  input:                                                 *****
!  =====                                                  *****
! - throw : temperature in k                              *****
! - roarow : air density                                  *****
! - shj    : local mid-layer sigma value                  *****
! - pressg : grid row of surface pressure [pa]            *****
! - temp2  : temperature at 10m. (deg k)                  *****
! - sutemp : surface temperature (deg k)                  *****
! - srad   : solar irradiance at the ground(w/m**2)       *****
! - rh10   : relative humidity of air at 10m.             *****
! -vegcover: vegetation cover                             *****
! - wind10 : wind at 10 m                                 *****
! - xrow   : dust concentration      [kg/kg]              *****
! - ustar  : u*                                           *****
! - ilev   : number of model level                        *****
!                                                         *****
! output                                                  *****
! ======                                                  *****
! - rtdry  : dry deposition tendency                      *****
!                                                         *****
!**************************************************************

      use mod_constants , only : ep2 , gti , rgti , vonkar , mathpi ,   &
                               & tzero , boltzk , stdpmb
      use mod_aero_param

      implicit none
!
! Dummy arguments
!
      integer :: il1 , il2 , ilev , ilg , luc , nbin
      integer , dimension(ilg) :: ivegcov
      real(8) , dimension(ilg,ilev,nbin) :: pdepv
      real(8) , dimension(ilg) :: pressg , rh10 , srad , sutemp ,       &
                                & temp2 , wind10 , zeff
      real(8) , dimension(ilg,ilev) :: roarow , throw
      real(8) , dimension(ilev) :: shj
      real(8) , dimension(nbin,2) :: trsize
      intent (in) il1 , il2 , ilev , ilg , ivegcov , luc , nbin ,       &
                & pressg , rh10 , roarow , shj , srad , sutemp , temp2 ,&
                & throw , trsize , wind10 , zeff
      intent (inout) pdepv
!
! Local variables
!
      real(8) :: amfp , amob , asq , ch , cm , cun , dtemp , dthv , eb ,&
               & eim , ein , es , fh , fm , frx1 , kui , logratio ,     &
               & mol , pre , prii , priiv , psit , psiu , ptemp2 , qs , &
               & r1 , ratioz , aa , rib , st , tbar , thstar , tsv ,    &
               & tsw , ustarsq , utstar , vp , vptemp , wvpm , ww , x , &
               & y , z , z0water , zdl , zl
      real(8) , dimension(ilg,ilev) :: amu
      real(8) , dimension(ilg) :: anu , schm , zz0
      real(8) , dimension(ilg,ilev,isize) :: cfac , pdepvsub , pdiff ,  &
           & rhsize , taurel
      integer :: i , j , jc , k , kcov , l , lev , n , tot
      real(8) , dimension(ilg,luc) :: ra , ustar , vegcover
      real(8) , dimension(ilg,luc,isize) :: rs
 
      real(8) , parameter :: z10 = 10.0
 
      do n = 1 , isize
        avesize(n) = (aerosize(1,n)+aerosize(2,n))/2.0
      end do
 
!======================================================================
!     ********************************************************
!     *   aerosize - dry radius                           ****
!     *   rhop  - density for each aerosol type           ****
!     ********************************************************
      do n = 1 , isize
        do l = 1 , ilev
          do i = il1 , il2
 
!           ********************************************************
!           *  aerosol gravitational settling velocity          ****
!           *  and diffusion coefficient                        ****
!           *                                                   ****
!           * air's dynamic viscosity                           ****
!           ********************************************************
 
            amu(i,l) = a1*1.E-8*throw(i,l)**a2/(throw(i,l)+a3)
 
!           . . . . mid layer pressure in [pascal].
            pre = pressg(i)*shj(l)
!           ********************************************************
!           * mean molecular free path.                         ****
!           *     k.v. beard [1976], j atm. sci., 33            ****
!           ********************************************************
 
            amfp = c1*(amu(i,l)/c2)*(c3/pre)*(throw(i,l)/c4)**(1./2.)
            prii = 2./9.*gti/amu(i,l)
            priiv = prii*(rhop-roarow(i,l))
 
!           ********************************************************
!           * cunningham slip correction factor and             ****
!           * relaxation time = vg/grav.                        ****
!           ********************************************************
 
            cfac(i,l,n) = 1. + amfp/avesize(n)                          &
                        & *(aa1+aa2*exp(-aa3*avesize(n)/amfp))
            taurel(i,l,n) = dmax1(priiv*avesize(n)**2*cfac(i,l,n)*rgti, &
                          & 0.D0)
 
!           ********************************************************
!           * stokes friction                                  *****
!           pdepvsub(i,l,n) ' sellting dep. velocity = '
!           ********************************************************
 
            pdepvsub(i,l,n) = taurel(i,l,n)*gti
          end do
        end do
      end do
!======================================================================
 
!======================================================================
!     ****************************************************
!     * ra : is the aerodynamic resistance above the  ****
!     *      canopy and it is function in u* and      ****
!     *      z0: roughness length and the stability   ****
!     *      function                                 ****
!     * mol  - monin obukhov length (m) - calculated  ****
!     *           for each land use category          ****
!     * ptemp2 -potential temperature at z2  (deg. k) ****
!     * temp2 - temperature at 10m. (deg k)           ****
!     * z10   - 10 m.                                 ****
!     * sutemp -surface temperature (deg k)           ****
!     * srad   -solar irradiance at the ground(w/m**2)****
!     * rh10  - relative humidity of air at 10m.      ****
!     *           (0.0-1.0)                           ****
!     * stdpmb - sea level pressure (mb)               ****
!     ****************************************************
      do j = 1 , luc
        do i = il1 , il2
          ww = dmax1(wind10(i),1.0D0)
          zz0(i) = zeff(i)
 
! ***************************************************************
!         *  potential temperature at z2  (deg. k)                 
!         *****
 
! ***************************************************************
          ptemp2 = temp2(i) + z10*0.0098
 
! ***************************************************************
!         * for calculations over water compute values of critical 
!         ***** * profile variables: l and ustar                       
!         ***** *           ******begin for water***                   
!         *****
 
! ***************************************************************
          if ( ivegcov(i).eq.0 ) then
 
! **************************************************************
!           *  vp  - vapour pressure at z2                          
!           ***** *  wvpm- water vapour mixing ratio at  z2            
!           ***** *  vptemp- virtual potential temperature at z2 (deg.
!           k)  *****
! **************************************************************
            es = 6.108*exp(17.27*(temp2(i)-tzero)/(temp2(i)-35.86))
            vp = rh10(i)*es
            wvpm = ep2*vp/(stdpmb-vp)
            vptemp = ptemp2*(1.0+0.61*wvpm)
 
! **************************************************************
!           *  assume rh10 at water surface is 100%                 
!           ***** *   vp = es(tsw-tzero) !sat. vap press at surface   
!           ***** *   saturated vapour pressure at surface             
!           ***** *   saturated mixing ratio at surface                
!           ***** *   tsv - virtual potential temperature at surface   
!           ***** *           (deg. k)                                 
!           *****
 
! **************************************************************
            tsw = sutemp(i)
            vp = 6.108*exp(17.27*(tsw-tzero)/(tsw-35.86))
            qs = ep2*vp/(stdpmb-vp)
            tsv = tsw*(1.+0.61*qs)
            z0water = 1.0E-4
! **************************************************************
!           * scalet  :  not required if  z2 = 10m                  
!           *****
! **************************************************************
            dthv = (vptemp-tsv)
! **************************************************************
!           * calculate drag coefficient cun with neutral condition 
!           ***** * assumption  garratt (1977)                         
!           *****
 
! **************************************************************
            cun = 7.5E-4 + 6.7E-5*ww
            mol = 9999.0
 
            if ( abs(dthv).gt.1.0E-6 )                                  &
               & mol = vptemp*cun**1.5*ww**2/(5.096E-3*dthv)
            if ( mol.gt.0. .and. mol.lt.5.0 ) mol = 5.0
            if ( mol.gt.-5.0 .and. mol.lt.0 ) mol = -5.0
            zdl = z10/mol
!
            if ( zdl.lt.0.0 ) then
 
! **************************************************************
!             *                        wind speed                     
!             *****
 
! **************************************************************
              x = (1.0-15.0*zdl)**0.25
              psiu = 2.*dlog(0.5*(1.0+x)) + dlog(0.5*(1.0+x*x))         &
                   & - 2.0*atan(x) + 0.5*mathpi
 
! **************************************************************
!             *                       pot temp                        
!             *****
! **************************************************************
              y = sqrt(1.-9.*zdl)
              psit = 2.*0.74*dlog((1+y)/2.0)
            else
              psiu = -4.7*zdl
              psit = psiu
            end if
            z0water = 0.000002*ww**2.5
!
            ustar(i,j) = vonkar*ww/(dlog(z10/z0water)-psiu)
            thstar = vonkar*(ptemp2-sutemp(i))                          &
                   & /(0.74*dlog(z10/z0water)-psit)
            zz0(i) = z0water
!
          else
 
! **************************************************************
!           * compute ustar and l for land use categories other than 
!           **** * water use louis method. !pkk 7/16/85, find bulk     
!           **** * richardson number.                                  
!           ****
 
! **************************************************************
            rib = gti*z10*(ptemp2-sutemp(i))/(sutemp(i)*ww**2)
 
! ***************************************************************
!           * ensure that conditions over land are never stable when 
!           ***** * there is incoming solar radiatiom                  
!           *****
! ***************************************************************
            if ( srad(i).gt.0.0 .and. rib.gt.0.0 ) rib = 1.E-15
!
            dtemp = ptemp2 - sutemp(i)
            if ( dabs(dtemp).lt.1.E-10 ) dtemp = dsign(1.D-10,dtemp)
            tbar = 0.5*(ptemp2+sutemp(i))
!
            ratioz = z10/zz0(i)
            logratio = dlog(ratioz)
            asq = 0.16/(logratio**2)
!
            if ( rib.le.0.0 ) then
              aa = asq*9.4*sqrt(ratioz)
              cm = 7.4*aa
              ch = 5.3*aa
              fm = 1. - (9.4*rib/(1.+cm*sqrt(abs(rib))))
              fh = 1. - (9.4*rib/(1.+ch*sqrt(abs(rib))))
            else
              fm = 1./((1.+4.7*rib)**2)
              fh = fm
            end if
!
            ustarsq = asq*ww**2*fm
            utstar = asq*ww*dtemp*fh/0.74
            ustar(i,j) = sqrt(ustarsq)
            thstar = utstar/ustar(i,j)
!
            mol = tbar*ustarsq/(vonkar*gti*thstar)
          end if
 
          kui = 1.0/(vonkar*ustar(i,j))
 
!         **************************************************************
!         * compute the values of  ra                            *******
!         **************************************************************
 
          z = z10
          zl = z/mol
 
          if ( zl.ge.0. ) then
            ra(i,j) = kui*(0.74*dlog(z/zz0(i))+4.7*zl)
          else
            ra(i,j) = kui*0.74*(dlog(z/zz0(i))-2.0*dlog((1+sqrt(1-9.*zl)&
                    & )*0.5))
          end if
          ra(i,j) = dmax1(ra(i,j),0.99D0)
          ra(i,j) = dmin1(ra(i,j),999.9D0)
        end do
      end do
!======================================================================
 
!======================================================================
!     find the right table index for the cell cover ( ocean and lake
!     are 0 in the ivegcov and 14-15 in the table )
 
      if ( ivegcov(i).eq.0 ) then
        kcov = 14
      else
        kcov = ivegcov(i)
      end if
 
!     *****************************************************
!     * the schmidt number is the ratio of the         ****
!     * kinematic viscosity of air to the particle     ****
!     * brownian diffusivity ===> sc=v/d               ****
!     *****************************************************
 
      do n = 1 , isize
        do l = 1 , ilev
          do i = il1 , il2
 
!           *****************************************************
!           * for now we will not consider the humidity      ****
!           * impact so we will set the variable frx1=1.0    ****
!           * i.e. only dry particles                        ****
!           *****************************************************
 
            frx1 = 1.0
            rhsize(i,l,n) = avesize(n)*frx1
            anu(i) = amu(i,l)/roarow(i,l)
            amob = 6.*mathpi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
            pdiff(i,l,n) = boltzk*throw(i,l)/amob
            schm(i) = anu(i)/pdiff(i,l,n)
 
!           ******************************************************
!           * for brownian diffusion, there is evidence that  ****
!           * its fromula depend on schmidt number as :       ****
!           * eb= schm x c^gama                               ****
!           * where gama is efficiency factor and its value   ****
!           * between 1/2 and 2/3 with larger values for      ****
!           * rougher surfaces                                ****
!           * ****************************************************
 
          end do
          if ( l.eq.ilev ) then
            do k = 1 , luc ! luc  = 1 for the moment
              do i = il1 , il2
 
!               ******************************************************
!               * the parameter governing impaction processes is *****
!               * the stokes number,st, which has the form of    *****
!               * 1) st = vg x u* /g a for vegetated surefaces   *****
!               *    (slinn, 1982)                               *****
!               * 2) st = vg x u*2/anu for smothed surfaces or   *****
!               *    surfaces with bluff roughness elements      *****
!               *    (giorgi,1988)                               *****
!               ******************************************************
 
                st = taurel(i,l,n)*ustar(i,k)*ustar(i,k)/anu(i)
 
                eb = schm(i)**(-0.666667)
!               eim=(st/(st+aest(k)))**2
                eim = (st/(st+aest(kcov)))**2
 
 
                eim = dmin1(eim,0.6D0)
                ein = 0.0
!               if (arye(k) .gt. 0.0001) then
!               ein = (1000.0*2.0*avesize(n)/arye(k))**1.5
!               end if
 
                if ( arye(kcov).gt.0.0001 )                             &
                   & ein = (1000.0*2.0*avesize(n)/arye(kcov))**1.5
 
                ein = dmin1(ein,0.5D0)
 
!               *****************************************************
!               * partickes larger than 5 micro may rebounded   *****
!               * after hitting a surface,this process may be   *****
!               * included by modifying the total collection    *****
!               * by the factor of r1, which represents the     *****
!               * fraction of particles sticking to the surface *****
!               * slinn (1982) suggested the following:         *****
!               * r = exp (- st^0.2)                            *****
!               *****************************************************
 
!               r1 = exp (-st**0.5)
                r1 = dmax1(0.5D0,exp(-st**0.5))
!               if (k .ge. 11 .and. r1 .lt. 0.5 ) r1=0.5
                if ( kcov.ge.11 .and. r1.lt.0.5 ) r1 = 0.5
                if ( r1.lt.0.4 ) r1 = 0.4
 
!               ***************************************************
!               * calculation of rs: the surface resistance   *****
!               * which depends on the collection efficiency  *****
!               * of the surface and is determined by the     *****
!               * various deposition processes                *****
!               ***************************************************
 
!               rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
                rs(i,k,n) = 1.0/3.0/ustar(i,k)/(eb+eim+ein)/r1
              end do
            end do
          end if
        end do
      end do
!======================================================================
 
!     note for the moment vegcover = 1 ( one type of cover ivegcov per
 
!     gride cell)
      do n = 1 , isize
        do jc = 1 , luc
          do i = il1 , il2
 
!           modify the surface layer
 
            vegcover(i,jc) = 1.
            pdepvsub(i,ilev,n) = pdepvsub(i,ilev,n) + vegcover(i,jc)    &
                               & *1.0/(ra(i,jc)+rs(i,jc,n))
 
          end do
        end do
 
      end do  ! end nsize loop
 
!     average deposition velocities on bin
      do k = 1 , nbin
        tot = 0
        do lev = 1 , ilev
          do i = 1 , ilg
            pdepv(i,lev,k) = 0.
          end do
        end do
        do n = 1 , isize
          if ( avesize(n)*1.E6.ge.trsize(k,1) .and. avesize(n)          &
             & *1.E6.lt.trsize(k,2) ) then
            do lev = 1 , ilev
              do i = 1 , ilg
                pdepv(i,lev,k) = pdepv(i,lev,k) + pdepvsub(i,lev,n)
              end do
            end do
            tot = tot + 1
          end if
        end do
        if ( tot.gt.0 ) then
          do lev = 1 , ilev
            do i = 1 , ilg
              pdepv(i,lev,k) = pdepv(i,lev,k)/tot
            end do
          end do
        end if
      end do
 
      end subroutine chdrydep
