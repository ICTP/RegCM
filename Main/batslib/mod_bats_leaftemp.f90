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
  use mod_bats_internal
!
  private
!
  public :: lftemp , satur
!
  contains
!
!=======================================================================
!     based on: bats version 1e          copyright 18 august 1989
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
!     taf(n,j,i) = temperature of air in canopy
!     delt(n,j,i)= difference between temperature of overlying air
!                          and that in canopy
!     delq(n,j,i)= difference between humidity of overlying air
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
  integer :: iter , itfull , itmax , n , i , j
  !
  !=======================================================================
  ! 1.   setup information
  !=======================================================================
  !
  ! 1.1  get stress-free stomatal resistance
  !      (1st guess at vapor pressure deficit)
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          if ( sigf(n,j,i) > 0.001D0 ) then
            vpdc(n,j,i) = d_10
            if ( lemiss ) then
              sgtg3 = emiss(n,j,i)*(sigm*tgrd(n,j,i)**d_three)
            else
              sgtg3 = sigm*tgrd(n,j,i)**d_three
            end if
            flneto(n,j,i) = d_four*sgtg3*(tlef(n,j,i)-tgrd(n,j,i))
          end if
        end if
      end do
    end do
  end do

  call stomat
  !
  ! 1.3  determine fraction of total and green canopy surface
  !      covered by water
  call frawat
  !
  ! 1.4  establish root function in terms of etrc = maximum
  !      sustainable transpiration rate
  ! (routine also returns efpr, used in subr. water to define upper soil
  !  layer transpiration)
  !
  call root
  ! 1.5  saturation specific humidity of leaf
  call satur(qsatl,tlef,sfcp)
 
!=======================================================================
!l    2.   begin iteration for leaf temperature calculation
!=======================================================================
  iter = 0
  efeb = d_zero
  delmax = d_one
  itmax = 10
  itfull = itmax
! itmax = 40
! itfull = 20
 
  do iter = 0 , itmax
    !
    ! 2.1  recalc stability dependent canopy & leaf drag coeffs
    !
    if ( iter == 0 ) call condch

    call lfdrag
    call condch
 
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              lftra(n,j,i) = d_one/(cf(n,j,i)*uaf(n,j,i))
              cn1(n,j,i) = wtlh(n,j,i)*rhs(n,j,i)
              df(n,j,i) = cn1(n,j,i)*cpd
              ! 
              ! 2.2  decrease foliage conductance for stomatal
              !      resistance
              rppdry = lftra(n,j,i)*fdry(n,j,i)/(lftrs(n,j,i)+lftra(n,j,i))
              rpp(n,j,i) = rppdry + fwet(n,j,i)
              ! 
              ! 2.3  recalculate saturation vapor pressure
              !
              eg1 = eg(n,j,i)
              eg(n,j,i) = c1es*dexp(lfta(n,j,i)*(tlef(n,j,i)-tzero)/      &
                        (tlef(n,j,i)-lftb(n,j,i)))
              qsatl(n,j,i) = qsatl(n,j,i)*eg(n,j,i)/eg1
            end if
          end if
        end do
      end do
    end do
 
    ! 2.4  canopy evapotranspiration
    if ( iter == 0 ) call condcq
 
    epss = 1.0D-10
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              efpot(n,j,i) = cn1(n,j,i)*(wtgaq(n,j,i)*qsatl(n,j,i) - &
                           wtgq0(n,j,i)*qgrd(n,j,i) - wtaq0(n,j,i)*qs(n,j,i))
              if ( efpot(n,j,i) > d_zero ) then
                etr(n,j,i) = efpot(n,j,i)*lftra(n,j,i)*fdry(n,j,i) / &
                           (lftrs(n,j,i)+lftra(n,j,i))
                rpp(n,j,i) = dmin1(rpp(n,j,i),(etr(n,j,i)+ldew(n,j,i)/      &
                           dtbat)/efpot(n,j,i)-epss)
              else
                etr(n,j,i) = d_zero
                rpp(n,j,i) = d_one
              end if
              if ( ( efpot(n,j,i) >= d_zero ) .and. &
                   ( etr(n,j,i) >= etrc(n,j,i) ) )  then
                ! transpiration demand exceeds supply, stomat adjust
                ! demand
                rppdry = lftra(n,j,i)*fdry(n,j,i)/(lftrs(n,j,i)+lftra(n,j,i))
                rppdry = rppdry/(etr(n,j,i)/etrc(n,j,i))
                etr(n,j,i) = etrc(n,j,i)
                ! recalculate stomatl resistance and rpp
                lftrs(n,j,i) = lftra(n,j,i)*(fdry(n,j,i)/rppdry-d_one)
                rpp(n,j,i) = rppdry + fwet(n,j,i)
                rpp(n,j,i) = dmin1(rpp(n,j,i),(etr(n,j,i)+ldew(n,j,i)/      &
                           dtbat)/efpot(n,j,i)-epss)
              end if
              rppq(n,j,i) = wlhv*rpp(n,j,i)
              efe(n,j,i) = rppq(n,j,i)*efpot(n,j,i)
              if ( efe(n,j,i)*efeb < d_zero ) efe(n,j,i) = 0.1D0*efe(n,j,i)
            end if
          end if
        end do
      end do
    end do
    !=======================================================================
    !      3.   solve for leaf temperature
    !=======================================================================
    !  3.1  update conductances for changes in rpp and cdr
    call condcq
    !
    !  3.2  derivatives of energy fluxes with respect to leaf
    !       temperature for newton-raphson calculation of leaf temperature.
    !  subr.  ii: rs,ra,cdrd,rppq,efe.
    !  subr. output: qsatld,dcd.
    if ( iter <= itfull ) call deriv
    !
    !  3.3  compute dcn from dcd, output from subr. deriv
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              dcn = dcd(n,j,i)*tlef(n,j,i)
              ! 1.2  radiative forcing for leaf temperature calculation
              if ( lemiss ) then
                sgtg3 = emiss(n,j,i)*(sigm*tgrd(n,j,i)**d_three)
              else
                sgtg3 = sigm*tgrd(n,j,i)**d_three
              end if
              sf1 = sigf(n,j,i)*(sabveg(j,i)-flw(j,i)-(d_one-sigf(n,j,i))* &
                    flneto(n,j,i)+d_four*sgtg3*tgrd(n,j,i))
              sf2 = d_four*sigf(n,j,i)*sgtg3 + &
                         df(n,j,i)*wtga(n,j,i) + dcd(n,j,i)
              ! 3.4  iterative leaf temperature calculation
              tbef(n,j,i) = tlef(n,j,i)
              tlef(n,j,i) = (sf1+df(n,j,i)*(wta0(n,j,i)*sts(n,j,i)+        &
                            wtg0(n,j,i)*tgrd(n,j,i))-efe(n,j,i)+dcn)/sf2
              !
              ! 3.5  chk magnitude of change; limit to max allowed value
              dels(n,j,i) = tlef(n,j,i) - tbef(n,j,i)
              if ( dabs(dels(n,j,i)) > delmax ) then
                tlef(n,j,i) = tbef(n,j,i) + &
                                delmax*dels(n,j,i)/dabs(dels(n,j,i))
              end if
              ! 3.6  update dependence of stomatal resistance
              !      on vapor pressure deficit
              qcan = wtlq0(n,j,i)*qsatl(n,j,i) + qgrd(n,j,i)*wtgq0(n,j,i) + &
                     qs(n,j,i)*wtaq0(n,j,i)
              vpdc(n,j,i) = (d_one-rpp(n,j,i))*(qsatl(n,j,i)-qcan)*d_1000/ep2
            end if
          end if
        end do
      end do
    end do
    call stomat
    ! 3.8  end iteration
  end do
 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          if ( sigf(n,j,i) > 0.001D0 ) then
            !==================================================================
            ! 4.   update dew accumulation (kg/m**2/s)
            !==================================================================
            ldew(n,j,i) = ldew(n,j,i) + (etr(n,j,i) - efe(n,j,i)/wlhv)*dtbat
            !=================================================================
            ! 5.   collect parameters needed to evaluate
            !      sensible and latent fluxes
            !=================================================================
            !  5.1  canopy properties
            taf(n,j,i) = wtg0(n,j,i)*tgrd(n,j,i) + wta0(n,j,i)*sts(n,j,i) + &
                         wtl0(n,j,i)*tlef(n,j,i)
            delt(n,j,i) = wtgl(n,j,i)*sts(n,j,i) - &
                            (wtl0(n,j,i)*tlef(n,j,i) + &
                             wtg0(n,j,i)*tgrd(n,j,i))
            delq(n,j,i) = wtglq(n,j,i)*qs(n,j,i) - &
                         (wtlq0(n,j,i)*qsatl(n,j,i) + &
                          wtgq0(n,j,i)*qgrd(n,j,i))
            if ( lemiss ) then
              sgtg3 = emiss(n,j,i)*(sigm*tgrd(n,j,i)**d_three)
            else
              sgtg3 = sigm*tgrd(n,j,i)**d_three
            end if
            flnet(n,j,i) = sgtg3*(tlef(n,j,i)-tgrd(n,j,i))*d_four
            xxkb = dmin1(rough(lveg(n,j,i)),d_one)
            vakb = (d_one-sigf(n,j,i))*vspda(n,j,i) + sigf(n,j,i) * &
                   (xxkb*uaf(n,j,i)+(d_one-xxkb)*vspda(n,j,i))
            wtg2(n,j,i) = (d_one-sigf(n,j,i))*cdr(n,j,i)*vakb
            fbare = wtg2(n,j,i)*(tgrd(n,j,i)-sts(n,j,i))
            qbare = wtg2(n,j,i)*(qgrd(n,j,i)-qs(n,j,i))
            !  5.2  fluxes from soil
            fseng(n,j,i) = cpd*rhs(n,j,i)*(wtg(n,j,i)*((wta0(n,j,i)+     &
                         wtl0(n,j,i))*tgrd(n,j,i)-wta0(n,j,i)*sts(n,j,i)- &
                         wtl0(n,j,i)*tlef(n,j,i))+fbare)
            fevpg(n,j,i) = rhs(n,j,i)*rgr(n,j,i) * &
                        (wtg(n,j,i)*((wtaq0(n,j,i)+  &
                         wtlq0(n,j,i))*qgrd(n,j,i)-wtaq0(n,j,i)*qs(n,j,i)- &
                         wtlq0(n,j,i)*qsatl(n,j,i))+qbare)
            !  5.3  deriv of soil energy flux with respect to soil temp
            qsatdg = qgrd(n,j,i)*rgr(n,j,i)*lfta(n,j,i)*(tzero-lftb(n,j,i)) * &
                     (d_one/(tgrd(n,j,i)-lftb(n,j,i)))**d_two
            cgrnds(n,j,i) = rhs(n,j,i)*cpd*(wtg(n,j,i)*(wta0(n,j,i) + &
                          wtl0(n,j,i))+wtg2(n,j,i))
            cgrndl(n,j,i) = rhs(n,j,i)*qsatdg*((wta(n,j,i)+wtlq(n,j,i)) * &
                          wtg(n,j,i)*wtsqi(n,j,i)+wtg2(n,j,i))
            cgrnd(n,j,i) = cgrnds(n,j,i) + cgrndl(n,j,i)*htvp(n,j,i)
            !  5.4  reinitialize cdrx
            !     shuttleworth mods #3 removed here !!!!!!
            cdrx(n,j,i) = cdr(n,j,i)
            !
            !  5.5  fluxes from canopy and soil to overlying air
            fbare = wtg2(n,j,i)*(tgrd(n,j,i)-sts(n,j,i))
            qbare = wtg2(n,j,i)*(qgrd(n,j,i)-qs(n,j,i))
            sent(n,j,i) = cpd*rhs(n,j,i)*(-wta(n,j,i)*delt(n,j,i)+fbare)
            evpr(n,j,i) = rhs(n,j,i)*(-wta(n,j,i)*delq(n,j,i) + &
                            rgr(n,j,i)*qbare)
          end if
        end if
      end do
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
!     seasb = fseas(tgbrd(n,j,i) (set in bndry) is a seasonal
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
!             for the lower canopy,i.e., trup=dexp(-0.5*g*rlai/coszrs),
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
  integer :: il , ilmax , n , i , j
  real(dp) , dimension(10) :: rad , radd
!
!     ***** seasonal temperature factor
!     ***** g is average leaf crosssection per unit lai
!     ***** difzen is ave of inverse of cos of angle of diffuse vis
!     light ***** ilmax is number of canopy layers
!     ***** coszrs is cosine solar zenith angle for incident light
!     *****   (to spec from input data need a good treatment of diffuse
!     rad) ***** trup is transmission of direct beam light in one
!     canopy layer ***** trupd is transmission of diffuse light in one
!     canopy layer
  g = d_half
  difzen = d_two
  ilmax = 4
  rilmax = d_four
 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          if ( sigf(n,j,i) > 0.001D0 ) then
            ! zenith angle set in zenitm
            if ( (coszrs(j,i)/rilmax) > 0.001D0 ) then
              trup(n,j,i) = dexp(-g*rlai(n,j,i)/(rilmax*coszrs(j,i)))
              trupd(n,j,i) = dexp(-difzen*g*rlai(n,j,i)/(rilmax))
              if ( trup(n,j,i) < dlowval )  trup(n,j,i)  = dlowval
              if ( trupd(n,j,i) < dlowval ) trupd(n,j,i) = dlowval
              fsold(n,j,i) = fracd(j,i)*solis(j,i)*fc(lveg(n,j,i))
              fsol0(n,j,i) = (d_one-fracd(j,i))*solis(j,i)*fc(lveg(n,j,i))
              rmini(n,j,i) = rsmin(lveg(n,j,i))/rmax0
            end if
          end if
        end if
      end do
    end do
  end do
 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          if ( sigf(n,j,i) > 0.001D0 ) then
            if ( coszrs(j,i)/rilmax > 0.001D0 ) then
              rad(1) = (d_one-trup(n,j,i))*fsol0(n,j,i)*rilmax/rlai(n,j,i)
              radd(1) = (d_one-trupd(n,j,i))*fsold(n,j,i) * &
                        rilmax/rlai(n,j,i)
              do il = 2 , ilmax
                rad(il) = trup(n,j,i)*rad(il-1)
                radd(il) = trupd(n,j,i)*radd(il-1)
              end do
              radfi = d_zero
              do il = 1 , ilmax
                radfi = radfi + (rad(il)+radd(il)+rmini(n,j,i)) / &
                        (d_one+rad(il)+radd(il))
              end do
              radf(n,j,i) = rilmax/radfi
            end if
          end if
        end if
      end do
    end do
  end do
 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          if ( sigf(n,j,i) > 0.001D0 ) then
            if ( (coszrs(j,i)/rilmax) > 0.001D0 ) then
              vpdf = d_one/dmax1(0.3D0,d_one-vpdc(n,j,i)*0.025D0)
              seas = d_one/(rmini(n,j,i)+fseas(tlef(n,j,i)))
              lftrs(n,j,i) = rsmin(lveg(n,j,i))*radf(n,j,i)*seas*vpdf
              lftrs(n,j,i) = dmin1(lftrs(n,j,i),rmax0)
            else
              lftrs(n,j,i) = rmax0
            end if
          end if
        end if
      end do
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
!  ldew(i) is in kg/m**2/s
!  fwet   = ratio of dew to max value to 2/3 power
!           ( 2/3 power comes from deardorff (1978) )
!              ** keep fwet le 1.0 **
!  dewmaxi = inverse of max allowed dew depth on leaf in mm
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine frawat
    implicit none
    integer :: n , i , j
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              fwet(n,j,i) = d_zero
              if ( ldew(n,j,i) > d_zero ) then
                fwet(n,j,i) = ((dewmaxi/vegt(n,j,i))*ldew(n,j,i))**twot
                fwet(n,j,i) = dmin1(fwet(n,j,i),d_one)
              end if
              fdry(n,j,i) = (d_one-fwet(n,j,i))*xlai(n,j,i)/xlsai(n,j,i)
            end if
          end if
        end do
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
    real(dp) :: bneg , rotf , trsmx , wlttb , wltub , wmli
    integer :: n , i , j
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              ! trsmx = trsmx0*sigf(n,j,i)*seasb(n,j,i)
              trsmx = trsmx0*sigf(n,j,i)
              rotf = rootf(lveg(n,j,i))
              bneg = -bsw(n,j,i)
              wmli = d_one/(wiltr(n,j,i)**bneg-d_one)
              wlttb = (watr(n,j,i)**bneg-d_one)*wmli
              wltub = (watu(n,j,i)**bneg-d_one)*wmli
              wlttb = dmin1(wlttb,d_one)
              wltub = dmin1(wltub,d_one)
              etrc(n,j,i) = trsmx*(d_one-(d_one-rotf)*wlttb-rotf*wltub)
              efpr(n,j,i) = trsmx*rotf*(d_one-wltub)
              if ( etrc(n,j,i) < 1.0D-12 ) then
                etrc(n,j,i) = 1.0D-12
                efpr(n,j,i) = d_one
              else
                efpr(n,j,i) = efpr(n,j,i)/etrc(n,j,i)
              end if
            end if
          end if
        end do
      end do
    end do
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
    implicit none
    real(dp) , pointer , dimension(:,:,:) :: p , qsat , t
    intent (in) p , t
    intent (out) qsat
    integer :: n , i , j
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( t(n,j,i) <= tzero ) then
            lfta(n,j,i) = c3ies
            lftb(n,j,i) = c4ies
          else
            lfta(n,j,i) = c3les
            lftb(n,j,i) = c4les
          end if
          eg(n,j,i) = c1es*dexp(lfta(n,j,i)*(t(n,j,i)-tzero) / &
                                        (t(n,j,i)-lftb(n,j,i)))
          qsat(n,j,i) = ep2*eg(n,j,i)/(p(n,j,i)-0.378D0*eg(n,j,i))
        end do
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
    implicit none
!
    real(dp) :: dthdz , ribi , sqrtf , tkb , u1 , u2 , zatild
    integer :: n , i , j
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              tkb = wta0(n,j,i)*sts(n,j,i) + wtl0(n,j,i)*tlef(n,j,i) + &
                    wtg0(n,j,i)*tgrd(n,j,i)
              dlstaf(n,j,i) = sts(n,j,i) - sigf(n,j,i)*tkb - &
                            (d_one-sigf(n,j,i))*tgrd(n,j,i)
              if ( dlstaf(n,j,i) <= d_zero ) then
                dthdz = (d_one-sigf(n,j,i))*tgrd(n,j,i) + &
                         sigf(n,j,i)*tkb-sts(n,j,i)
                u1 = wtur + d_two*dsqrt(dthdz)
                ribd(n,j,i) = usw(j,i)**d_two + vsw(j,i)**d_two + u1**d_two
              else
                u2 = wtur
                ribd(n,j,i) = usw(j,i)**d_two + vsw(j,i)**d_two + u2**d_two
              end if
              vspda(n,j,i) = dsqrt(ribd(n,j,i))
              if ( vspda(n,j,i) < d_one ) then
                vspda(n,j,i) = d_one
                ribd(n,j,i) = d_one
              end if
              zatild = (zh(n,j,i)-displa(lveg(n,j,i)))*sigf(n,j,i) + &
                        zh(n,j,i)*(d_one-sigf(n,j,i))
              rib1(n,j,i) = egrav*zatild/(ribd(n,j,i)*sts(n,j,i))
              rib(n,j,i) = rib1(n,j,i)*dlstaf(n,j,i)
              if ( rib(n,j,i) < d_zero ) then
                cdr(n,j,i) = cdrn(n,j,i)*(d_one+24.5D0 * &
                                          dsqrt(-cdrn(n,j,i)*rib(n,j,i)))
                sqrtf = dmin1(dsqrt(-cdrn(n,j,i)/rib(n,j,i)),11.5D0/12.25D0)
                cdrd(n,j,i) = cdrn(n,j,i)*12.25D0* &
                            wtl0(n,j,i)*rib1(n,j,i)*sigf(n,j,i)*sqrtf
              else
                ribi = d_one/(d_one+11.5D0*rib(n,j,i))
                cdr(n,j,i) = cdrn(n,j,i)*ribi
                cdrd(n,j,i) = cdr(n,j,i)*ribi*11.5D0 * &
                              rib1(n,j,i)*wtl0(n,j,i)*sigf(n,j,i)
                cdrmin(n,j,i) = dmax1(cdrn(n,j,i)*d_rfour,6.0D-4)
              end if
              if ( (rib(n,j,i) >= d_zero) ) then
                if ( (cdr(n,j,i) < cdrmin(n,j,i)) ) then
                  cdr(n,j,i) = cdrmin(n,j,i)
                  cdrd(n,j,i) = d_zero
                end if
              end if
            end if
          end if
        end do
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
    implicit none
    integer :: n , i , j
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
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              uaf(n,j,i) = vspda(n,j,i)*dsqrt(cdr(n,j,i))
              cf(n,j,i) = 0.01D0*sqrtdi(lveg(n,j,i))/dsqrt(uaf(n,j,i))
              wta(n,j,i) = sigf(n,j,i)*cdr(n,j,i)*vspda(n,j,i)
              wtlh(n,j,i) = cf(n,j,i)*uaf(n,j,i)*vegt(n,j,i)
              wtg(n,j,i) = csoilc*uaf(n,j,i)*sigf(n,j,i)
              wtshi(n,j,i) = d_one/(wta(n,j,i)+wtlh(n,j,i)+wtg(n,j,i))
              wtl0(n,j,i) = wtlh(n,j,i)*wtshi(n,j,i)
              wtg0(n,j,i) = wtg(n,j,i)*wtshi(n,j,i)
              wtgl(n,j,i) = wtl0(n,j,i) + wtg0(n,j,i)
              wta0(n,j,i) = d_one - wtgl(n,j,i)
              wtga(n,j,i) = wta0(n,j,i) + wtg0(n,j,i)
            end if
          end if
        end do
      end do
    end do
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
    implicit none
    integer :: n , i , j
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
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              rgr(n,j,i) = gwet(n,j,i)
              wtlq(n,j,i) = wtlh(n,j,i)*rpp(n,j,i)
              wtgq(n,j,i) = wtg(n,j,i)*rgr(n,j,i)
              wtsqi(n,j,i) = d_one/(wta(n,j,i)+wtlq(n,j,i)+wtgq(n,j,i))
              wtgq0(n,j,i) = wtgq(n,j,i)*wtsqi(n,j,i)
              wtlq0(n,j,i) = wtlq(n,j,i)*wtsqi(n,j,i)
              wtglq(n,j,i) = wtgq0(n,j,i) + wtlq0(n,j,i)
              wtaq0(n,j,i) = d_one - wtglq(n,j,i)
              wtgaq(n,j,i) = wtaq0(n,j,i) + wtgq0(n,j,i)
            end if
          end if
        end do
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
    implicit none
    real(dp) :: dne , hfl , xkb
    integer :: n , i , j
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
            if ( sigf(n,j,i) > 0.001D0 ) then
              dne = d_one/(tlef(n,j,i)-lftb(n,j,i))
              qsatld(n,j,i) = qsatl(n,j,i)*lfta(n,j,i) * &
                            (tzero-lftb(n,j,i))*dne**d_two
              xkb = cdrd(n,j,i)/cdr(n,j,i)
              hfl = df(n,j,i)*(wtga(n,j,i)*tlef(n,j,i) - &
                             wtg0(n,j,i)*tgrd(n,j,i)   - &
                             wta0(n,j,i)*sts(n,j,i))
              dcd(n,j,i) = cn1(n,j,i)*rppq(n,j,i)*wtgaq(n,j,i) *   &
                         qsatld(n,j,i) + (d_one-wtgaq(n,j,i)) *   &
                         efe(n,j,i) * xkb + (d_one-wtga(n,j,i)) * hfl * xkb
              dcd(n,j,i) = dmax1(dcd(n,j,i),d_zero)
              dcd(n,j,i) = dmin1(dcd(n,j,i),500.0D0)
            end if
          end if
        end do
      end do
    end do
  end subroutine deriv
!
end module mod_bats_leaftemp
