C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine chdrydep(ilg,il1,il2,ilev,luc,nbin, ivegcov,
     &         throw, roarow, shj, pressg,temp2, sutemp,
     &         srad, rh10,wind10,zeff,trsize,
     &         pdepv )   
c     ***************************************************************
c     *  dry depostion scheme for dust particles                *****
c     *  this scheme based on :                                 *****
c     * - zhang et al,(2001) : a size-segregated particle       *****
c     *   dry deposition scheme for an atmospheric aerosol      *****
c     *   module, atmos. env. 35, 549-560                       *****
c     *                                                         *****
c     * - giorgi, f. (1986): a particle dry deposition          *****
c     *      parameterization scheme for use in tracer          *****
c     *      transport models. jgr,91, 9794-9806                *****
c     *                                                         *****
c     *  input:                                                 *****
c     *  =====                                                  *****
c     * - throw : temperature in k                              *****
c     * - roarow : air density                                  *****
c     * - shj    : local mid-layer sigma value                  *****
c     * - pressg : grid row of surface pressure [pa]            *****
c     * - temp2  : temperature at 10m. (deg k)                  *****
c     * - sutemp : surface temperature (deg k)                  *****
c     * - srad   : solar irradiance at the ground(w/m**2)       *****
c     * - rh10   : relative humidity of air at 10m.             *****
c     * -vegcover: vegetation cover                             *****
c     * - wind10 : wind at 10 m                                 *****
c     * - xrow   : dust concentration      [kg/kg]              *****
c     * - ustar  : u*                                           *****
c     * - ilev   : number of model level                        *****
c     *                                                         *****
c     * output                                                  *****
c     * ======                                                  *****
c     * - rtdry  : dry deposition tendency                      *****
c     *                                                         *****
c     ***************************************************************

      implicit none
      integer n,l,i,k,j,jc
      integer isize,ilev,il1,il2,ilg
      integer luc,lev,nbin
      parameter(isize = 12)    
c subbbin
      real(kind=8)  shj(ilev),pressg(ilg),throw(ilg,ilev)
      real(kind=8)  roarow(ilg,ilev)
      real(kind=8)  pdepvsub(ilg,ilev,isize)
      real(kind=8)  taurel(ilg,ilev,isize)
      real(kind=8)  aerosize (2,isize)
      real(kind=8)  ustar(ilg,luc)
      real(kind=8)  cfac(ilg,ilev,isize)
      real(kind=8)  amu(ilg,ilev)
      real(kind=8)  schm(ilg)
      real(kind=8)  rs(ilg, luc,isize) 
      real(kind=8)  temp2(ilg), sutemp(ilg)
      real(kind=8)  wind10(ilg), rh10(ilg)
      real(kind=8)  srad(ilg)
      real(kind=8)  ra(ilg,luc),vegcover(ilg,luc)
!      real(kind=8)  arw
      real(kind=8)  avesize(isize)  
      integer ivegcov(ilg) ,tot
      real(kind=8)  zeff(ilg) 
      real(kind=8)  trsize(nbin,2)
      real(kind=8)  pdepv(ilg,ilev,nbin)

      data    aerosize/ 1.0e-08,   2.0e-08,
     &                  2.0e-08,   4.0e-08,
     &                  4.0e-08,   8.0e-08,
     &                  8.0e-08,   1.6e-07,
     &                  1.6e-07,   3.2e-07,
     &                  3.2e-07,   6.4e-07,
     &                  6.4e-07,   1.28e-06,
     &                  1.28e-06,  2.56e-06,
     &                  2.56e-06,  5.12e-06,
     &                  5.12e-06, 10.4e-06,
     &                  10.24e-06,20.48e-06,
     &                  20.48e-06,40.6e-06/ 

      real(kind=8)  pre,prii,priiv,amfp
      real(kind=8)  a1,a2,a3,c1,c2,c3,c4,aa1,aa2,aa3,g,boltzk,rhop,pi
      data    a1/145.8/, a2/1.5/, a3/110.4/
      data    c1/6.54e-8/,c2/1.818e-5/,c3/1.013e5/,c4/293.15/
      data    aa1/1.257/, aa2/0.4/, aa3/1.1/
      data    g/9.806160/,boltzk/1.3806044e-23/
      data    rhop/2650.0/,pi/3.14159265/

      real(kind=8)  arye(20),aest(20)
      data    aest/0.80,  0.80, 0.8, 0.8, 1.2,
     1             1.20,  2.00, 1.5, 1.5, 2.0,
     2            15.0 , 15.0 , 1.5, 1.5, 1.5,
     5            15.0 ,  1.20, 1.2, 1.2, 1.2 /
      data    arye/0.5   ,5.0   ,0.5 , 5.0 , 1.0,
     1             1.0   ,0.0001,5.0 ,10.0 ,10.0,
     2             0.0001,0.0001,0.56, 0.56, 0.56,
     3             0.56  ,0.56  ,0.56, 0.56, 0.56/
      integer kcov
      real(kind=8)  frx1,amob,eb,st, eim, ein, r1
      real(kind=8)  anu(ilg),rhsize(ilg,ilev,isize)
      real(kind=8)  pdiff(ilg,ilev,isize)

      real(kind=8)  ww, mol,ptemp2,tbar,ratioz,logratio,asq,dtemp,vp,es
      real(kind=8)  wvpm,ustarsq,utstar,thstar,kui,aa,cm,ch,fm,fh
      real(kind=8)  rib,vptemp,tsw,qs,tsv,z0water,dthv,cun,zdl
      real(kind=8)  psiu,x,z,zl,y,psit,zz0(ilg)

      real(kind=8)  krm,z10,pmb
      data    krm/0.4/,z10/10.0/,pmb /1000.0/

      do n=1,isize
        avesize(n) = (aerosize(1,n)+aerosize(2,n))/2.0
      end do

c======================================================================
c  ********************************************************
c  *   aerosize - dry radius                           ****
c  *   rhop  - density for each aerosol type           ****
c  ********************************************************
      do n=1,isize
         do l=1,ilev
           do i=il1,il2  

c  ********************************************************
c  *  aerosol gravitational settling velocity          ****
c  *  and diffusion coefficient                        ****
c  *                                                   ****
c  * air's dynamic viscosity                           ****
c  ********************************************************

             amu(i,l)=a1*1.e-8*throw(i,l)**a2/(throw(i,l)+a3)

c . . . . mid layer pressure in [pascal].
             pre = pressg(i)*shj(l)
c  ********************************************************
c  * mean molecular free path.                         ****
c  *     k.v. beard [1976], j atm. sci., 33            ****
c  ********************************************************

             amfp=c1*(amu(i,l)/c2)*(c3/pre)*(throw(i,l)/c4)**(1./2.)
             prii=2./9.*g/amu(i,l)
             priiv=prii*(rhop-roarow(i,l))

c  ********************************************************
c  * cunningham slip correction factor and             ****
c  * relaxation time = vg/grav.                        ****
c  ********************************************************

             cfac(i,l,n) = 1.+amfp/avesize(n)*(aa1+aa2*
     &            exp(-aa3*avesize(n)/amfp))
             taurel(i,l,n)=dmax1(priiv*avesize(n)**2*cfac(i,l,n)/g,0.d0)

c  ********************************************************
c  * stokes friction                                  *****  
c    pdepvsub(i,l,n) ' sellting dep. velocity = '
c  ********************************************************

             pdepvsub(i,l,n) =taurel(i,l,n)*g
           end do
         end do
      end do
c======================================================================

c======================================================================
c     ****************************************************
c     * ra : is the aerodynamic resistance above the  ****
c     *      canopy and it is function in u* and      ****
c     *      z0: roughness length and the stability   ****
c     *      function                                 ****
c     * mol  - monin obukhov length (m) - calculated  ****
c     *           for each land use category          ****
c     * ptemp2 -potential temperature at z2  (deg. k) ****
c     * temp2 - temperature at 10m. (deg k)           ****
c     * z10   - 10 m.                                 ****
c     * sutemp -surface temperature (deg k)           ****
c     * srad   -solar irradiance at the ground(w/m**2)****
c     * rh10  - relative humidity of air at 10m.      ****
c     *           (0.0-1.0)                           ****
c     * pmb   - sea level pressure (mb)               ****
c     ****************************************************
      do j=1,luc
        do i=il1,il2
          ww = dmax1(wind10(i),1.0d0)
          zz0(i) = zeff(i)              

c     ***************************************************************
c     *  potential temperature at z2  (deg. k)                  ***** 
c     ***************************************************************

          ptemp2 = temp2(i) + z10 * 0.0098

c     ***************************************************************
c     * for calculations over water compute values of critical  *****
c     * profile variables: l and ustar                          *****
c     *           ******begin for water***                      *****
c     ***************************************************************

          if ( ivegcov(i) .eq. 0) then 
   
c     **************************************************************
c     *  vp  - vapour pressure at z2                           *****
c     *  wvpm- water vapour mixing ratio at  z2                *****
c     *  vptemp- virtual potential temperature at z2 (deg. k)  *****
c     **************************************************************
            es = 6.108*exp(17.27*(temp2(i)-273.16)/(temp2(i)-35.86))
            vp = rh10(i) * es
            wvpm = 0.622 * vp/(pmb-vp)           
            vptemp = ptemp2 * (1.0 + 0.61 * wvpm)

c     **************************************************************
c     *  assume rh10 at water surface is 100%                  *****
c     *   vp = es(tsw-273.15) !sat. vap press at surface       *****
c     *   saturated vapour pressure at surface                 *****
c     *   saturated mixing ratio at surface                    *****
c     *   tsv - virtual potential temperature at surface       *****
c     *           (deg. k)                                     *****
c     **************************************************************
      
            tsw= sutemp(i)
            vp = 6.108*exp(17.27*(tsw - 273.16)/(tsw - 35.86))
            qs = 0.622 * vp / (pmb - vp)
            tsv = tsw * (1. + 0.61 * qs)
            z0water = 1.0e-4
c     **************************************************************
c     * scalet  :  not required if  z2 = 10m                   *****
c     **************************************************************
            dthv = (vptemp - tsv)
c     **************************************************************
c     * calculate drag coefficient cun with neutral condition  *****
c     * assumption  garratt (1977)                             ***** 
c     **************************************************************
                
            cun=7.5e-4+6.7e-5*ww
            mol=9999.0

      if (abs(dthv).gt.1.0e-6)mol=vptemp*cun**1.5*ww**2/(5.096e-3*dthv)
            if (mol.gt.0..and.mol.lt.5.0) mol = 5.0
            if (mol.gt.-5.0.and.mol.lt.0) mol =-5.0
            zdl = z10/mol
c
            if (zdl.lt.0.0) then

c     **************************************************************
c     *                        wind speed                      *****
c     ************************************************************** 

              x=(1.0 - 15.0*zdl)**0.25
              psiu=2.*dlog(0.5*(1.0+x))+dlog(0.5*(1.0+x*x))-
     &                     2.0*atan(x) + 0.5*pi

c     **************************************************************
c     *                       pot temp                         *****
c     **************************************************************
              y=sqrt(1.-9.*zdl)                   
              psit=2.*0.74*dlog((1+y)/2.0)
            else
              psiu = -4.7*zdl
              psit = psiu
            end if
            z0water = 0.000002 *ww**2.5
c
            ustar(i,j)= krm*ww/(dlog(z10/z0water) - psiu)
            thstar = krm*(ptemp2-sutemp(i))/(0.74*
     &                   dlog(z10/z0water)-psit)
            zz0(i)=z0water
c
          else

c     **************************************************************
c     * compute ustar and l for land use categories other than  ****
c     * water use louis method. !pkk 7/16/85, find bulk         **** 
c     * richardson number.                                      ****
c     **************************************************************

            rib = g*z10*(ptemp2 - sutemp(i))/(sutemp(i)*ww**2)

c     ***************************************************************
c     * ensure that conditions over land are never stable when  ***** 
c     * there is incoming solar radiatiom                       *****
c     *************************************************************** 
            if (srad(i).gt.0.0.and.rib.gt.0.0) rib = 1.e-15
c
            dtemp = ptemp2 - sutemp(i)
            if (dabs(dtemp).lt.1.e-10) dtemp=dsign(1.d-10,dtemp)
            tbar = 0.5*(ptemp2 + sutemp(i))
c
            ratioz = z10/zz0(i)
            logratio = dlog(ratioz)
            asq = 0.16/(logratio**2)
c
            if (rib.le.0.0)then
              aa = asq*9.4*sqrt(ratioz)
              cm = 7.4*aa
              ch = 5.3*aa
              fm = 1. - (9.4*rib/(1. + cm*sqrt(abs(rib))))
              fh = 1. - (9.4*rib/(1. + ch*sqrt(abs(rib))))
            else
              fm = 1./((1. + 4.7*rib)**2)
              fh = fm
            end if
c
            ustarsq = asq*ww**2*fm
            utstar = asq*ww*dtemp*fh/0.74
            ustar(i,j) = sqrt(ustarsq)
            thstar = utstar/ustar(i,j)
c
            mol = tbar*ustarsq/(0.4*9.8*thstar)
          end if

          kui=1.0/(krm*ustar(i,j))

c     **************************************************************  
c     * compute the values of  ra                            *******
c     **************************************************************
 
          z=z10
          zl=z/mol

          if (zl.ge.0.) then
            ra(i,j)=kui*(0.74*dlog(z/zz0(i))+4.7*zl)
          else
            ra(i,j)=kui*0.74*(dlog(z/zz0(i))-
     &                         2.0*dlog((1+sqrt(1-9.*zl))*0.5))
          end if
          ra(i,j)=dmax1(ra(i,j),0.99d0)
          ra(i,j)=dmin1(ra(i,j),999.9d0)
        end do
      end do 
c======================================================================

c======================================================================
c find the right table index for the cell cover ( ocean and lake are 0
c      in the ivegcov and 14-15 in the table )

      if (ivegcov(i) .eq. 0) then 
        kcov = 14
      else 
        kcov = ivegcov(i)
      end if
 
c     *****************************************************
c     * the schmidt number is the ratio of the         ****
c     * kinematic viscosity of air to the particle     ****
c     * brownian diffusivity ===> sc=v/d               ****
c     *****************************************************
         
      do 555 n=1,isize
        do 666 l=1,ilev
          do 777 i=il1,il2

c     *****************************************************
c     * for now we will not consider the humidity      ****
c     * impact so we will set the variable frx1=1.0    ****
c     * i.e. only dry particles                        ****
c     *****************************************************

            frx1 = 1.0
            rhsize(i,l,n)=avesize(n)*frx1
            anu(i)=amu(i,l)/roarow(i,l)
            amob=6.*pi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
            pdiff(i,l,n)=boltzk*throw(i,l)/amob
            schm(i)=anu(i)/pdiff(i,l,n)

c    ******************************************************
c    * for brownian diffusion, there is evidence that  ****
c    * its fromula depend on schmidt number as :       ****
c    * eb= schm x c^gama                               ****
c    * where gama is efficiency factor and its value   ****
c    * between 1/2 and 2/3 with larger values for      ****
c    * rougher surfaces                                ****
c    * ****************************************************
              
 777      continue
          if (l .eq. ilev) then
            do 400 k=1,luc ! luc  = 1 for the moment
              do 200 i = il1, il2
         
c    ******************************************************
c    * the parameter governing impaction processes is *****
c    * the stokes number,st, which has the form of    *****
c    * 1) st = vg x u* /g a for vegetated surefaces   *****
c    *    (slinn, 1982)                               *****
c    * 2) st = vg x u*2/anu for smothed surfaces or   *****
c    *    surfaces with bluff roughness elements      *****
c    *    (giorgi,1988)                               *****
c    ******************************************************

                st = taurel(i,l,n)* ustar(i,k)*ustar(i,k)/anu(i)

                eb = schm(i)**(-0.666667)
c               eim=(st/(st+aest(k)))**2 
                eim=(st/(st+aest(kcov)))**2 

              
                eim=dmin1(eim,0.6d0)
                ein = 0.0
c               if (arye(k) .gt. 0.0001) then
c                  ein = (1000.0*2.0*avesize(n)/arye(k))**1.5
c               end if

                if (arye(kcov) .gt. 0.0001) then
                  ein = (1000.0*2.0*avesize(n)/arye(kcov))**1.5
                end if

                ein = dmin1(ein,0.5d0)

c     *****************************************************
c     * partickes larger than 5 micro may rebounded   *****
c     * after hitting a surface,this process may be   *****
c     * included by modifying the total collection    *****
c     * by the factor of r1, which represents the     *****
c     * fraction of particles sticking to the surface *****
c     * slinn (1982) suggested the following:         *****
c     * r = exp (- st^0.2)                            *****
c     ***************************************************** 
              
c               r1 = exp (-st**0.5)
                r1= dmax1(0.5d0,  exp(-st**0.5))
c               if (k .ge. 11 .and. r1 .lt. 0.5 ) r1=0.5
                if ( kcov .ge. 11 .and. r1 .lt. 0.5 ) r1=0.5
                if (r1 .lt. 0.4) r1=0.4
                 
c      ***************************************************
c      * calculation of rs: the surface resistance   *****
c      * which depends on the collection efficiency  *****
c      * of the surface and is determined by the     *****
c      * various deposition processes                *****
c      ***************************************************

c               rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
                rs(i,k,n) = 1.0/3.0/ustar(i,k)/(eb+eim+ein)/r1
 200          continue                                    
 400        continue
          end if
 666    continue
 555  continue
c======================================================================

c note for the moment vegcover = 1 ( one type of cover ivegcov per gride cell) 

      do n = 1,isize 
        do 12 jc = 1,luc
          do 13  i= il1, il2

c modify the surface layer 

            vegcover(i,jc) =1.
            pdepvsub(i,ilev,n) =  pdepvsub(i,ilev,n) +
     &        vegcover(i,jc) *  1.0/(ra(i,jc)+rs(i,jc,n) )
 
13        continue
12      continue      

      end do  ! end nsize loop

c average deposition velocities on bin 
      do k = 1,nbin
        tot = 0
        do lev=1,ilev
        do i=1,ilg
          pdepv(i,lev,k) =0.
        end do
        end do
        do n=1,isize                      
          if (  avesize(n)*1.e6 .ge.trsize(k,1) .and. 
     &          avesize(n)*1.e6 .lt.trsize(k,2)) then
            do lev=1,ilev
            do i=1,ilg
              pdepv(i,lev,k)=pdepv(i,lev,k)+pdepvsub(i,lev,n) 
            end do
            end do
            tot = tot + 1                 
          end if 
        end do
        if (tot.gt.0) then
          do lev=1,ilev
          do i=1,ilg
            pdepv(i,lev,k) = pdepv(i,lev,k) / tot
          end do
          end do
        end if
      end do

      return
      end subroutine chdrydep
