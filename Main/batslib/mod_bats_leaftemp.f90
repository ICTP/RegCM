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
 
module mod_bats_leaftemp
!
!     Calculate leaf temperature, leaf fluxes, and net transpiration.
!     documented in NCAR Tech Note, Dickinson et al., 1986.
!     modifications by Klaus Blumel, 1988.
!
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_bats_common
!
  private
!
  public :: lfta , lftb , lftrs , lftra
  public :: allocate_mod_bats_leaftemp , lftemp , satur
!
  real(dp) , pointer , dimension(:,:) :: lfta , lftb
  real(dp) , pointer , dimension(:,:) :: lftra , lftrs
!
  real(dp) , pointer , dimension(:,:) :: cdrd , vpdc
  real(dp) , pointer , dimension(:,:) :: rppq , efe
  real(dp) , pointer , dimension(:,:) :: dcd , etrc
  real(dp) , pointer , dimension(:,:) :: qsatld
  real(dp) , pointer , dimension(:,:) :: dels
  real(dp) , pointer , dimension(:,:) :: efpot , tbef
  real(dp) , pointer , dimension(:,:) :: fsol0 , fsold
  real(dp) , pointer , dimension(:,:) :: radf , rmini
  real(dp) , pointer , dimension(:,:) :: trup , trupd
  real(dp) , pointer , dimension(:,:) :: cdrmin , dlstaf
  real(dp) , pointer , dimension(:,:) :: rib , rib1
!
  contains
!
  subroutine allocate_mod_bats_leaftemp
    implicit none

    call getmem2d(lfta,1,nnsg,1,iym1,'mod_leaftemp:lfta')
    call getmem2d(lftb,1,nnsg,1,iym1,'mod_leaftemp:lftb')
    call getmem2d(lftra,1,nnsg,1,iym1,'mod_leaftemp:lftra')
    call getmem2d(lftrs,1,nnsg,1,iym1,'mod_leaftemp:lftrs')
    call getmem2d(cdrd,1,nnsg,1,iym1,'mod_leaftemp:cdrd')
    call getmem2d(vpdc,1,nnsg,1,iym1,'mod_leaftemp:vpdc')
    call getmem2d(rppq,1,nnsg,1,iym1,'mod_leaftemp:rppq')
    call getmem2d(efe,1,nnsg,1,iym1,'mod_leaftemp:efe')
    call getmem2d(qsatld,1,nnsg,1,iym1,'mod_leaftemp:qsatld')
    call getmem2d(dcd,1,nnsg,1,iym1,'mod_leaftemp:dcd')
    call getmem2d(etrc,1,nnsg,1,iym1,'mod_leaftemp:etrc')
    call getmem2d(dels,1,nnsg,1,iym1,'mod_leaftemp:dels')
    call getmem2d(efpot,1,nnsg,1,iym1,'mod_leaftemp:efpot')
    call getmem2d(tbef,1,nnsg,1,iym1,'mod_leaftemp:tbef')
    call getmem2d(fsol0,1,nnsg,1,iym1,'mod_leaftemp:fsol0')
    call getmem2d(fsold,1,nnsg,1,iym1,'mod_leaftemp:fsold')
    call getmem2d(radf,1,nnsg,1,iym1,'mod_leaftemp:radf')
    call getmem2d(rmini,1,nnsg,1,iym1,'mod_leaftemp:rmini')
    call getmem2d(trup,1,nnsg,1,iym1,'mod_leaftemp:trup')
    call getmem2d(trupd,1,nnsg,1,iym1,'mod_leaftemp:trupd')
    call getmem2d(cdrmin,1,nnsg,1,iym1,'mod_leaftemp:cdrmin')
    call getmem2d(dlstaf,1,nnsg,1,iym1,'mod_leaftemp:dlstaf')
    call getmem2d(rib,1,nnsg,1,iym1,'mod_leaftemp:rib')
    call getmem2d(rib1,1,nnsg,1,iym1,'mod_leaftemp:rib1')

  end subroutine allocate_mod_bats_leaftemp
!
!=======================================================================
!l  based on: bats version 1e          copyright 18 august 1989
!=======================================================================
!
!     Calculate leaf temperature, leaf fluxes, and net transpiration.
!     documented in NCAR Tech Note, Dickinson et al., 1986.
!     modifications by Klaus Blumel, 1988.
!
!        f l o w   d i a g r a m   f o r   l e f t e m
!
!              lftemp ===> stomat
!                          frawat
!                            root
!                           satur
!                          lfdrag
!                          condch
!                          condcq
!                           deriv
!
!     cf = heat transfer coeff. from leaves; assumes same for moisture
!          (see e.g. d. gates' book) for laminar flow past leaf;
!          from ewing paper, cf has dimensions t**1/2 l**-1;
!          so to make it dimensionless, premultiply by 0.01 (si).
!     cgrnd = deriv. of soil energy flux with respect to soil temp.
!                             (used in tgrund)
!     delt = ts-taf
!     ef = transpiration rate
!     efpot = potential evaporation rate (kg/m**2/s)
!     fevpg = evaporative heat flux from ground
!     flnet = (temp gradient)*(d/dt(sig t**4)); hence, t**3 term
!     fseng = sensible heat flux from ground
!     iter = leaf temp iteration counter; runs from 1 to max of 100
!     ra = leaf aerodynamic resistance factor
!     taf = air temperature within foliage canopy
!     tbef = leaf temp at time step before current one
!     ts = air temperature of lowest model layer
!     uaf = mean wind within canopy
!
!     taf1d(n,i) = temperature of air in canopy
!     delt1d(n,i)= difference between temperature of overlying air
!                          and that in canopy
!     delq1d(n,i)= difference between humidity of overlying air
!                          and that in canopy
!
!     convergence of leaf temperature calculation is declared if
!     enough iterations (itmin) and change of temp small enough and
!     change of latent heat fluxes small enough, or if
!     maximum iteration reached (itmax).
!
!=======================================================================
!
  subroutine lftemp
!
  implicit none
!
  real(dp) :: dcn , delmax , efeb , eg1 , epss , fbare , qbare ,     &
           & qcan , qsatdg , rppdry , sf1 , sf2 , sgtg3 , vakb ,    &
           & xxkb
  integer :: iter , itfull , itmax , n , i
!
!=======================================================================
!l    1.   setup information
!=======================================================================
!
!l    1.1  get stress-free stomatal resistance
!     (1st guess at vapor pressure deficit)
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          vpdc(n,i) = d_10
          if ( lemiss ) then
            sgtg3 = emiss1d(n,i)*(sigm*tg1d(n,i)**d_three)
          else
            sgtg3 = sigm*tg1d(n,i)**d_three
          end if
          flneto(n,i) = d_four*sgtg3*(tlef1d(n,i)-tg1d(n,i))
        end if
      end if
    end do
  end do
  call stomat
!
!l    1.3  determine fraction of total and green canopy surface
!l    covered by water
  call frawat
!
!l    1.4  establish root function in terms of etrc = maximum
!l    sustainable transpiration rate
!     (routine also returns efpr, used in subr. water to
!     define upper soil layer transpiration)
  call root
 
!l    1.5  saturation specific humidity of leaf
  call satur(qsatl,tlef1d,p1d)
 
!=======================================================================
!l    2.   begin iteration for leaf temperature calculation
!=======================================================================
  iter = 0
  efeb = d_zero
  delmax = d_one
  itmax = 10
  itfull = itmax
!     itmax = 40
!     itfull = 20
 
  do iter = 0 , itmax
!
!l      2.1  recalc stability dependent canopy & leaf drag coeffs
    if ( iter == 0 ) call condch
    call lfdrag
    call condch
 
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            lftra(n,i) = d_one/(cf(n,i)*uaf(n,i))
            cn1(n,i) = wtlh(n,i)*rhs1d(n,i)
            df(n,i) = cn1(n,i)*cpd
 
!l              2.2  decrease foliage conductance for stomatal
!               resistance
            rppdry = lftra(n,i)*fdry(n,i)/(lftrs(n,i)+lftra(n,i))
            rpp(n,i) = rppdry + fwet(n,i)
 
!l              2.3  recalculate saturation vapor pressure
            eg1 = eg(n,i)
            eg(n,i) = c1es*dexp(lfta(n,i)*(tlef1d(n,i)-tzero)/      &
                     & (tlef1d(n,i)-lftb(n,i)))
            qsatl(n,i) = qsatl(n,i)*eg(n,i)/eg1
          end if
        end if
      end do
    end do
 
!l      2.4  canopy evapotranspiration
    if ( iter == 0 ) call condcq
 
    epss = 1.0D-10
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            efpot(n,i) = cn1(n,i)*(wtgaq(n,i)*qsatl(n,i) - &
                                      wtgq0(n,i)*qg1d(n,i) -  &
                                      wtaq0(n,i)*qs1d(n,i))
 
!as             if (efpot(n,i) >= 0.) then     !if 0 rpp could have
!               floating pt
            if ( efpot(n,i) > d_zero ) then
              etr(n,i) = efpot(n,i)*lftra(n,i)*fdry(n,i) / &
                         (lftrs(n,i)+lftra(n,i))
              rpp(n,i) = dmin1(rpp(n,i),(etr(n,i)+ldew1d(n,i)/      &
                        & dtbat)/efpot(n,i)-epss)
            else
              etr(n,i) = d_zero
              rpp(n,i) = d_one
            end if
 
            if ( ( efpot(n,i) >= d_zero ) .and. &
                 ( etr(n,i) >= etrc(n,i) ) )  then
!*                transpiration demand exceeds supply, stomat adjust
!                 demand
              rppdry = lftra(n,i)*fdry(n,i)/(lftrs(n,i)+lftra(n,i))
              rppdry = rppdry/(etr(n,i)/etrc(n,i))
              etr(n,i) = etrc(n,i)
!*                recalculate stomatl resistance and rpp
              lftrs(n,i) = lftra(n,i)*(fdry(n,i)/rppdry-d_one)
              rpp(n,i) = rppdry + fwet(n,i)
              rpp(n,i) = dmin1(rpp(n,i),(etr(n,i)+ldew1d(n,i)/      &
                        & dtbat)/efpot(n,i)-epss)
            end if
 
            rppq(n,i) = wlhv*rpp(n,i)
            efe(n,i) = rppq(n,i)*efpot(n,i)
            if ( efe(n,i)*efeb < d_zero ) &
              efe(n,i) = 0.1D0*efe(n,i)
          end if
        end if
      end do
    end do
!=======================================================================
!l      3.   solve for leaf temperature
!=======================================================================
!l      3.1  update conductances for changes in rpp and cdr
    call condcq
!
!l      3.2  derivatives of energy fluxes with respect to leaf
!l      temperature for newton-raphson calculation of
!l      leaf temperature.
!l      subr.  ii: rs,ra,cdrd,rppq,efe.
!l      subr. output: qsatld,dcd.
    if ( iter <= itfull ) call deriv
!
!l      3.3  compute dcn from dcd, output from subr. deriv
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            dcn = dcd(n,i)*tlef1d(n,i)
!
!l              1.2  radiative forcing for leaf temperature calculation
            if ( lemiss ) then
              sgtg3 = emiss1d(n,i)*(sigm*tg1d(n,i)**d_three)
            else
              sgtg3 = sigm*tg1d(n,i)**d_three
            end if
            sf1 = sigf(n,i)*(sabveg(i)-flw1d(i)-(d_one-sigf(n,i))* &
                  flneto(n,i)+d_four*sgtg3*tg1d(n,i))
            sf2 = d_four*sigf(n,i)*sgtg3 + df(n,i)*wtga(n,i) + &
                  dcd(n,i)
 
!l              3.4  iterative leaf temperature calculation
            tbef(n,i) = tlef1d(n,i)
            tlef1d(n,i) = (sf1+df(n,i)*(wta0(n,i)*ts1d(n,i)+        &
                   & wtg0(n,i)*tg1d(n,i))-efe(n,i)+dcn)/sf2
!
!l              3.5  chk magnitude of change; limit to max allowed value
            dels(n,i) = tlef1d(n,i) - tbef(n,i)
            if ( dabs(dels(n,i)) > delmax )                     &
              &   tlef1d(n,i) = tbef(n,i) + delmax*dels(n,i)/ &
              &   dabs(dels(n,i))
 
!l              3.6  update dependence of stomatal resistance
!l              on vapor pressure deficit
            qcan = wtlq0(n,i)*qsatl(n,i) + qg1d(n,i)*wtgq0(n,i)     &
                 & + qs1d(n,i)*wtaq0(n,i)
            vpdc(n,i) = (d_one-rpp(n,i))*(qsatl(n,i)-qcan)*d_1000/ep2
          end if
        end if
      end do
    end do
 
    call stomat
 
!l      3.8  end iteration
 
  end do
 
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
!=======================================================================
!l            4.   update dew accumulation (kg/m**2/s)
!=======================================================================
          ldew1d(n,i) = ldew1d(n,i) + (etr(n,i) - &
                        efe(n,i)/wlhv)*dtbat
 
!=======================================================================
!l            5.   collect parameters needed to evaluate
!l            sensible and latent fluxes
!=======================================================================
 
!l            5.1  canopy properties
          taf1d(n,i) = wtg0(n,i)*tg1d(n,i) + wta0(n,i)              &
                      & *ts1d(n,i) + wtl0(n,i)*tlef1d(n,i)
          delt1d(n,i) = wtgl(n,i)*ts1d(n,i) - (wtl0(n,i)*           &
                       & tlef1d(n,i)+wtg0(n,i)*tg1d(n,i))
          delq1d(n,i) = wtglq(n,i)*qs1d(n,i) - (wtlq0(n,i)*         &
                       & qsatl(n,i)+wtgq0(n,i)*qg1d(n,i))
          if ( lemiss ) then
            sgtg3 = emiss1d(n,i)*(sigm*tg1d(n,i)**d_three)
          else
            sgtg3 = sigm*tg1d(n,i)**d_three
          end if
          flnet(n,i) = sgtg3*(tlef1d(n,i)-tg1d(n,i))*d_four
          xxkb = dmin1(rough(lveg(n,i)),d_one)
          vakb = (d_one-sigf(n,i))*vspda(n,i) + sigf(n,i) &
               & *(xxkb*uaf(n,i)+(d_one-xxkb)*vspda(n,i))
          wtg2(n,i) = (d_one-sigf(n,i))*cdr(n,i)*vakb
          fbare = wtg2(n,i)*(tg1d(n,i)-ts1d(n,i))
          qbare = wtg2(n,i)*(qg1d(n,i)-qs1d(n,i))
 
!l            5.2  fluxes from soil
          fseng(n,i) = cpd*rhs1d(n,i)*(wtg(n,i)*((wta0(n,i)+        &
                      & wtl0(n,i))*tg1d(n,i)-wta0(n,i)*ts1d(n,i)-   &
                      & wtl0(n,i)*tlef1d(n,i))+fbare)
          fevpg(n,i) = rhs1d(n,i)*rgr(n,i)*(wtg(n,i)*((wtaq0(n,i)+  &
                      & wtlq0(n,i))*qg1d(n,i)-wtaq0(n,i)*qs1d(n,i)- &
                      & wtlq0(n,i)*qsatl(n,i))+qbare)
 
!l            5.3  deriv of soil energy flux with respect to soil temp
          qsatdg = qg1d(n,i)*rgr(n,i)*lfta(n,i)*(tzero-lftb(n,i))   &
                 & *(d_one/(tg1d(n,i)-lftb(n,i)))**d_two
          cgrnds(n,i) = rhs1d(n,i)*cpd*(wtg(n,i)*(wta0(n,i)+        &
                 & wtl0(n,i))+wtg2(n,i))
          cgrndl(n,i) = rhs1d(n,i)*qsatdg*((wta(n,i)+wtlq(n,i))*    &
                 & wtg(n,i)*wtsqi(n,i)+wtg2(n,i))
          cgrnd(n,i) = cgrnds(n,i) + cgrndl(n,i)*htvp(n,i)
 
!l            5.4  reinitialize cdrx
!!!           shuttleworth mods #3 removed here !!!!!!
          cdrx(n,i) = cdr(n,i)
!
!l            5.5  fluxes from canopy and soil to overlying air
          fbare = wtg2(n,i)*(tg1d(n,i)-ts1d(n,i))
          qbare = wtg2(n,i)*(qg1d(n,i)-qs1d(n,i))
          sent1d(n,i) = cpd*rhs1d(n,i)*(-wta(n,i)*delt1d(n,i)+fbare)
          evpr1d(n,i) = rhs1d(n,i)*(-wta(n,i)*delq1d(n,i)+          &
                 & rgr(n,i)*qbare)
        end if
      end if
    end do
  end do
 
  end subroutine lftemp
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     gives leaf stomatal resistance from environmental parameters
!             under conditions of no moisture stress
!
!     standard lai from xla=max & xlai0=min lai
!     seasb = fseas(tgb1d(n,i) (set in bndry) is a seasonal
!             factor for reduced winter lai and root water uptake
!        fc = light sensitivity for crops and grasses and has inverse
!             radiation units (m**2/watt)
!      rlai = sum of leaf and stem area indices
!     rmax0 = 5000. s/m (maximum resistance)
!     radu & radl = visible light intensity in upper & lower canopy
!     ft & fb = the fractional intercepted photo-active radiation
!             per unit (leaf & stem) area in the top (upper) and
!             bottom (lower) canopies, respectively
!     radfi = average of upper and lower canopy light factors
!        rs = stomatal resistance = min.res. * rad.factor * leaf factor
!      trup = transmission of the upper canopy, assumed to be the same
!             for the lower canopy,i.e., trup=dexp(-0.5*g*rlai/czen),
!             where g = attenuation factor
!
!     documented in NCAR Tech Note, Dickinson et al., 1986
!     improved stomatal shading, Dickinson, nov 88.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine stomat
!
  implicit none
!
  real(dp) :: difzen , g , radfi , seas , vpdf , rilmax
  integer :: il , ilmax , n , i
  real(dp) , dimension(10) :: rad , radd
!
!     ***** seasonal temperature factor
!     ***** g is average leaf crosssection per unit lai
!     ***** difzen is ave of inverse of cos of angle of diffuse vis
!     light ***** ilmax is number of canopy layers
!     ***** czen is cosine solar zenith angle for incident light
!     *****   (to spec from input data need a good treatment of diffuse
!     rad) ***** trup is transmission of direct beam light in one
!     canopy layer ***** trupd is transmission of diffuse light in one
!     canopy layer
  g = d_half
  difzen = d_two
  ilmax = 4
  rilmax = d_four
 
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
!             **********            zenith angle set in zenitm
          if ( (czen(i)/rilmax) > 0.001D0 ) then
            trup(n,i) = dexp(-g*rlai(n,i)/(rilmax*czen(i)))
            trupd(n,i) = dexp(-difzen*g*rlai(n,i)/(rilmax))
            fsold(n,i) = fracd(i)*solis(i)*fc(lveg(n,i))
            fsol0(n,i) = (d_one-fracd(i))*solis(i)*fc(lveg(n,i))
            rmini(n,i) = rsmin(lveg(n,i))/rmax0
          end if
        end if
      end if
    end do
  end do
 
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          if ( czen(i)/rilmax > 0.001D0 ) then
            rad(1) = (d_one-trup(n,i))*fsol0(n,i)*rilmax/rlai(n,i)
            radd(1) = (d_one-trupd(n,i))*fsold(n,i) * &
                      rilmax/rlai(n,i)
            do il = 2 , ilmax
              rad(il) = trup(n,i)*rad(il-1)
              radd(il) = trupd(n,i)*radd(il-1)
            end do
            radfi = d_zero
            do il = 1 , ilmax
              radfi = radfi + (rad(il)+radd(il)+rmini(n,i)) / &
                      (d_one+rad(il)+radd(il))
            end do
            radf(n,i) = rilmax/radfi
          end if
        end if
      end if
    end do
  end do
 
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          if ( (czen(i)/rilmax) > 0.001D0 ) then
            vpdf = d_one/dmax1(0.3D0,d_one-vpdc(n,i)*0.025D0)
            seas = d_one/(rmini(n,i)+fseas(tlef1d(n,i)))
            lftrs(n,i) = rsmin(lveg(n,i))*radf(n,i)*seas*vpdf
            lftrs(n,i) = dmin1(lftrs(n,i),rmax0)
          else
            lftrs(n,i) = rmax0
          end if
        end if
      end if
    end do
  end do
!
  contains

  function fseas(x)
    implicit none
    real(dp) :: fseas
    real(dp) , intent(in) :: x
    fseas = dmax1(d_zero,d_one-0.0016D0*dmax1(298.0D0-x,d_zero)**d_two)
  end function fseas
 
  end subroutine stomat
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  **** determines fraction of foliage covered by water fwet, and
!  **** the fraction of foliage that is dry transpiring leaf fdry.
!  note: their defns differ - fwet is the fraction of all veg surfaces
!  which are wet because mod_stems can evaporate, fdry is the fraction
!  of lai which is dry because only leaves can transpire
!
!  ldew1d(i) is in kg/m**2/s
!  fwet   = ratio of dew to max value to 2/3 power
!           ( 2/3 power comes from deardorff (1978) )
!              ** keep fwet le 1.0 **
!  dewmaxi = inverse of max allowed dew depth on leaf in mm
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine frawat

  implicit none
!
  integer :: n , i
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          fwet(n,i) = d_zero
          if ( ldew1d(n,i) > d_zero ) then
            fwet(n,i) = ((dewmaxi/vegt(n,i))*ldew1d(n,i))**twot
            fwet(n,i) = dmin1(fwet(n,i),d_one)
          end if
          fdry(n,i) = (d_one-fwet(n,i))*xlai(n,i)/xlsai(n,i)
        end if
      end if
    end do
  end do
 
  end subroutine frawat
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     this subroutine provides root function in terms of maximum
!     transpiration rate plants can sustain depending on soil moisture.
!
!     trsmx0 is a prescribed constant (kg/m**2/s).
!     trsmx is the maximum transpiration rate,
!        including a low temperature correction (=seasb)
!        and a correction for fractional vegetation (=sigf).
!
!     rotf is ratio of moisture extracton from top to total when
!           fully saturated
!     rootf is ratio of roots in upper soil layer
!                    to roots in root soil layer
!     bsw is the b param in clapp and hornberger
!
!     "wlt  " are ratios factors controlling the saturation
!                 cf wilting (see ewing paper)
!     wlttb (total) & wltub (upper) become 1 at the wilting point
!     (eqn 14 in ewing paper) n.b. etrc=etrmx in ewing paper
!
!     etrc= max poss transpiration given the soil moisture distributions
!     efpr = the relative contribution of upper soil layer to
!     evapotranspiration - need soil moist. budget (subrout water)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine root
  implicit none
!
  real(dp) :: bneg , rotf , trsmx , wlttb , wltub , wmli
  integer :: n , i
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
!             trsmx = trsmx0*sigf(n,i)*seasb(n,i)
          trsmx = trsmx0*sigf(n,i)
          rotf = rootf(lveg(n,i))
          bneg = -bsw(n,i)
          wmli = d_one/(wiltr(n,i)**bneg-d_one)
          wlttb = (watr(n,i)**bneg-d_one)*wmli
          wltub = (watu(n,i)**bneg-d_one)*wmli
          wlttb = dmin1(wlttb,d_one)
          wltub = dmin1(wltub,d_one)
          etrc(n,i) = trsmx*(d_one-(d_one-rotf)*wlttb-rotf*wltub)
          efpr(n,i) = trsmx*rotf*(d_one-wltub)
          if ( etrc(n,i) < 1.0D-12 ) then
            etrc(n,i) = 1.0D-12
            efpr(n,i) = d_one
          else
            efpr(n,i) = efpr(n,i)/etrc(n,i)
          end if
        end if
      end if
    end do
  end do
!
  end subroutine root
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****  calculates saturation vapor pressure (eg)
!           and qsat = saturated specific humidity (dimensionless)
!
!           uses tetens formula (1930) (ref. riegel,1974,jam,p606
!                                                 equation 1)
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine satur(qsat,t,p)
!
  implicit none
!
  real(dp) , dimension(nnsg,iym1) :: p , qsat , t
  intent (in) p , t
  intent (out) qsat
!
  integer :: n , i
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( t(n,i) <= tzero ) then
        lfta(n,i) = c3ies
        lftb(n,i) = c4ies
      else
        lfta(n,i) = c3les
        lftb(n,i) = c4les
      end if
      eg(n,i) = c1es*dexp(lfta(n,i)*(t(n,i)-tzero) / &
                                    (t(n,i)-lftb(n,i)))
      qsat(n,i) = ep2*eg(n,i)/(p(n,i)-0.378D0*eg(n,i))
    end do
  end do
 
  end subroutine satur
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     recalculate stability dependent drag coefficient for vegetation,
!     given the neutral drag coefficient.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine lfdrag
!
  implicit none
!
  real(dp) :: dthdz , ribi , sqrtf , tkb , u1 , u2 , zatild
  integer :: n , i
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          tkb = wta0(n,i)*ts1d(n,i) + wtl0(n,i)*tlef1d(n,i)         &
              & + wtg0(n,i)*tg1d(n,i)
          dlstaf(n,i) = ts1d(n,i) - sigf(n,i)                       &
                       & *tkb - (d_one-sigf(n,i))*tg1d(n,i)
          if ( dlstaf(n,i) <= d_zero ) then
            dthdz = (d_one-sigf(n,i))*tg1d(n,i) + sigf(n,i)         &
                  & *tkb - ts1d(n,i)
            u1 = wtur + d_two*dsqrt(dthdz)
            ribd(n,i) = us1d(i)**d_two + vs1d(i)**d_two + u1**d_two
          else
            u2 = wtur
            ribd(n,i) = us1d(i)**d_two + vs1d(i)**d_two + u2**d_two
          end if
          vspda(n,i) = dsqrt(ribd(n,i))
          if ( vspda(n,i) < d_one ) then
            vspda(n,i) = d_one
            ribd(n,i) = d_one
          end if
          zatild = (z1d(n,i)-displa(lveg(n,i)))*sigf(n,i)            &
                 & + z1d(n,i)*(d_one-sigf(n,i))
          rib1(n,i) = egrav*zatild/(ribd(n,i)*ts1d(n,i))
          rib(n,i) = rib1(n,i)*dlstaf(n,i)
          if ( rib(n,i) < d_zero ) then
            cdr(n,i) = cdrn(n,i)*(d_one+24.5D0*dsqrt(-cdrn(n,i)*    &
                  & rib(n,i)))
            sqrtf = dmin1(dsqrt(-cdrn(n,i)/rib(n,i)),11.5D0/12.25D0)
            cdrd(n,i) = cdrn(n,i)*12.25D0*wtl0(n,i)*rib1(n,i)      &
                       & *sigf(n,i)*sqrtf
          else
            ribi = d_one/(d_one+11.5D0*rib(n,i))
            cdr(n,i) = cdrn(n,i)*ribi
            cdrd(n,i) = cdr(n,i)*ribi*11.5D0*rib1(n,i)*wtl0(n,i)   &
                       & *sigf(n,i)
            cdrmin(n,i) = dmax1(cdrn(n,i)*d_rfour,6.0D-4)
          end if
          if ( (rib(n,i) >= d_zero) ) then
            if ( (cdr(n,i) < cdrmin(n,i)) ) then
              cdr(n,i) = cdrmin(n,i)
              cdrd(n,i) = d_zero
            end if
          end if
        end if
      end if
    end do
  end do
 
  end subroutine lfdrag
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional sensible heat conductances
!               for canopy and soil flux calculations
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine condch
!
  implicit none
!
  integer :: n , i
!
!     csoilc = constant drag coefficient for soil under canopy
!     symbols used for weights are:   wt : weight
!     a : air
!     l : leaf
!     i : inverse
!     s : sum
!     h : sensible heat
!     q : water vapor
!     0 : normalized (sums to one)
!     g : ground
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          uaf(n,i) = vspda(n,i)*dsqrt(cdr(n,i))
          cf(n,i) = 0.01D0*sqrtdi(lveg(n,i))/dsqrt(uaf(n,i))
          wta(n,i) = sigf(n,i)*cdr(n,i)*vspda(n,i)
          wtlh(n,i) = cf(n,i)*uaf(n,i)*vegt(n,i)
          wtg(n,i) = csoilc*uaf(n,i)*sigf(n,i)
          wtshi(n,i) = d_one/(wta(n,i)+wtlh(n,i)+wtg(n,i))
          wtl0(n,i) = wtlh(n,i)*wtshi(n,i)
          wtg0(n,i) = wtg(n,i)*wtshi(n,i)
          wtgl(n,i) = wtl0(n,i) + wtg0(n,i)
          wta0(n,i) = d_one - wtgl(n,i)
          wtga(n,i) = wta0(n,i) + wtg0(n,i)
        end if
      end if
    end do
  end do
!
  end subroutine condch
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     dimensional and non-dimensional latent heat conductances
!               for canopy and soil flux calculations
!
!     latent fluxes differ from sensible due to stomatal resistance
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine condcq
!
  implicit none
!
  integer :: n , i
!
!     symbols used for weights are:   wt : weight
!     a : air
!     l : leaf
!     i : inverse
!     s : sum
!     h : sensible heat
!     q : water vapor
!     0 : normalized (sums to one)
!     g : ground
!
  do i = 2 , iym1
    do n = 1 , nnsg
      if ( ldoc1d(n,i) /= 0 ) then
        if ( sigf(n,i) > 0.001D0 ) then
          rgr(n,i) = gwet1d(n,i)
          wtlq(n,i) = wtlh(n,i)*rpp(n,i)
          wtgq(n,i) = wtg(n,i)*rgr(n,i)
          wtsqi(n,i) = d_one/(wta(n,i)+wtlq(n,i)+wtgq(n,i))
          wtgq0(n,i) = wtgq(n,i)*wtsqi(n,i)
          wtlq0(n,i) = wtlq(n,i)*wtsqi(n,i)
          wtglq(n,i) = wtgq0(n,i) + wtlq0(n,i)
          wtaq0(n,i) = d_one - wtglq(n,i)
          wtgaq(n,i) = wtaq0(n,i) + wtgq0(n,i)
        end if
      end if
    end do
  end do
 
  end subroutine condcq
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     derivatives of energy fluxes with respect to leaf temperature for
!     newton-raphson calculation of leaf temperature.
!     input: rs,ra,cdrd,rppq,efe.    output: qsatld,dcd.
!
!     approximate by derivatives of cdr and ef.  many weaker
!     dependences on leaf temperature are omitted, as convergence
!     rate is not affected.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine deriv
!
  implicit none
!
    real(dp) :: dne , hfl , xkb
    integer :: n , i
!
    do i = 1 , iym1
      do n = 1 , nnsg
        if ( ldoc1d(n,i) /= 0 ) then
          if ( sigf(n,i) > 0.001D0 ) then
            dne = d_one/(tlef1d(n,i)-lftb(n,i))
            qsatld(n,i) = qsatl(n,i)*lfta(n,i) * &
                          (tzero-lftb(n,i))*dne**d_two
            xkb = cdrd(n,i)/cdr(n,i)
            hfl = df(n,i)*(wtga(n,i)*tlef1d(n,i) - &
                           wtg0(n,i)*tg1d(n,i)   - &
                           wta0(n,i)*ts1d(n,i))
            dcd(n,i) = cn1(n,i)*rppq(n,i)*wtgaq(n,i) *   &
                          qsatld(n,i) + (d_one-wtgaq(n,i)) *   &
                          efe(n,i) * xkb + (d_one-wtga(n,i)) * &
                          hfl * xkb
            dcd(n,i) = dmax1(dcd(n,i),d_zero)
            dcd(n,i) = dmin1(dcd(n,i),500.0D0)
          end if
        end if
      end do
    end do
 
  end subroutine deriv
!
end module mod_bats_leaftemp
