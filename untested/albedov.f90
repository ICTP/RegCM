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
 
      subroutine albedov(j,iemiss)
 
      use mod_regcm_param
      use mod_bats , only : npts , albvs , albvl , aldirs , aldirl ,    &
                  & aldifs , aldifl , czen , sice1d , emiss1d , coszrs ,&
                  & ldoc1d , tgb1d , lveg , ts1d , scv1d , sag1d , wt , &
                  & scvk , veg1d , emiss2d , albvgs , albvgl , kolsol , &
                  & depuv , ssw1d
      use mod_constants , only : tmelt
      implicit none
!
! Dummy arguments
!
      integer :: iemiss , j
      intent (in) iemiss , j
!
! Local variables
!
      real(8) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               & albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               & conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               & dralbs , fsol1 , fsol2 , sfac , sical0 , sical1 , sl , &
               & sl2 , sli , snal0 , snal1 , tdiff , tdiffs , wet , x
      real(8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 & aldifs_s , aldirl_s , aldirs_s
      real(8) :: fseas
      integer :: kolour , n , np
      real(8) , dimension(8) :: solour
!
!     Albedo calculates fragmented albedos (direct and diffuse) in
!     wavelength regions split at 0.7um.
!
!     CM hands albedos to radiation package which computes
!     fsw1d(np) = net solar absorbed over full grid square
!     sabveg(np) = vegetation absorbed (full solar spectrum)
!     solis(np) = shortwave  solar incident
!
!     Here these are calculated at the end of albedo - they use only
!     direct albedos for now
!
!     in both versions :  lftemp uses sabveg
!     tgrund uses sabveg & fsw1d(np) to get
!     ground absorbed solar
!     photosynthesis uses solis - see subrouts
!     stomat and co2 (carbon)
!
!     For sea, sea-ice veg albedos are not set
!     these albedos are not treated as arrays here
!
!     (depuv/10.0)= the ratio of upper soil layer to total
!     root depth; used to compute "wet" for soil albedo
!     solour        =  soil albedo for different coloured soils
!
!     data    solour/ .12, .11, .10, .09, .08, .07, .06, .05 /
      data solour / 0.16D0 , 0.11D0 , 0.10D0 , 0.09D0 , 0.08D0 ,        &
                    0.07D0 , 0.06D0 , 0.05D0 /
 
!     =================================================================
!     1. set initial parameters
!     =================================================================
      fseas(x) = dmax1(0.D0,1.D0-0.0016D0*dmax1(298.D0-x,0.D0)**2)
!
!     1.1 constants
!
!     Solar flux partitioned at wavelength of 0.7micr
      fsol1 = 0.5D0
      fsol2 = 0.5D0
!     Short and long wave albedo for new snow
      snal0 = 0.95D0
      snal1 = 0.65D0
!     Short and long wave albedo for sea ice
      sical0 = 0.6D0
      sical1 = 0.4D0
 
!
!     In depth, wt is frac of grid square covered by snow;
!     depends on average snow depth, vegetation, etc.
!
      npts = nbmax
 
      call depth
 
!     1.2 set pointers
!     ***************************************************
!     *    set n"x"k params here  in ccm but not needed *
!     ***************************************************
!
!     1.3  set default vegetation and albedo
!     do loop 50 in ccm not used here )
 
      do np = np1 , npts
        czen(np) = dmax1(coszrs(np),0.D0)
        czeta = czen(np)
        do n = 1 , nnsg
          albgs = 0.0D0
          albgl = 0.0D0
          albgsd = 0.0D0
          albgld = 0.0D0
          albs = 0.0D0
          albl = 0.0D0
          albsd = 0.0D0
          albld = 0.0D0
 
          albvs_s(n) = 0.0D0
          albvl_s(n) = 0.0D0
 
!================================================================
!         2.   get albedo over land
!================================================================
!         can't use pointer "nalbk" here because not set - use nldock
!         instead tgb1d(np) used instead of tbelow
!
          if ( ldoc1d(n,np).gt.0.1D0 .and. sice1d(n,np).eq.0.D0 ) then
            sfac = 1.D0 - fseas(tgb1d(n,np))
 
!           **********  ccm tests here on land mask for veg and soils
 
 
!c          data *** reduces albedo at low temps !!!!!should respond to
!c          moisture too the following card inactivated (commented out)
!c          (pat, 27 oct 86)
 
!           veg1d(np)=vegc(lveg(np))-seasf(lveg(np))*sfac
            albs = albvgs(lveg(n,np))
            albl = albvgl(lveg(n,np))
 
!----------------------------------------------------------------------
            if ( (lveg(n,np).lt.12) .or. (lveg(n,np).gt.15) ) then
 
!             2.1  bare soil albedos
!             (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,np))
              wet = ssw1d(n,np)/depuv(lveg(n,np))
              alwet = dmax1((11.D0-40.D0*wet),0.D0)*0.01D0
              alwet = dmin1(alwet,solour(kolour))
              albg = solour(kolour) + alwet
!             if((lveg(n,np).eq.8)) albg=0.40      !Laura, cambiato il
!             DESERTO
              albgs = albg
              albgl = 2.D0*albg
!             **********            higher nir albedos
!             **********              set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl
 
!             Dec. 15   albzn=0.85+1./(1.+10.*czen(np))
!             Dec. 12, 2008
              albzn = 1.0D0
!             Dec. 15, 2008
 
!             **********            leafless hardwood canopy: no or
!             inverse zen dep
              if ( lveg(n,np).eq.5 .and. sfac.lt.0.1 ) albzn = 1.
!             **********            multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
 
!             **********            albedo over vegetation after zenith
!             angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
 
            else if ( lveg(n,np).eq.12 ) then
 
!             2.2   permanent ice sheet
              albgs = 0.8D0
              albgsd = 0.8D0
              albgl = 0.55D0
              albgld = 0.55D0
            else
 
!             2.3  inland water, swamps, rice paddies etc.
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgsd = albg
              albgl = albg
              albgld = albg
            end if
 
          else if ( sice1d(n,np).gt.0.D0 ) then
!====================================================================
!           3.  get albedo over sea ice
!====================================================================
!           **********          albedo depends on wave-length and ts.
!           the ts **********          dependence accounts for melt
!           water puddles.
            tdiffs = ts1d(n,np) - tmelt
            tdiff = dmax1(tdiffs,0.D0)
            tdiffs = dmin1(tdiff,20.D0)
            albgl = sical1 - 1.1E-2*tdiffs
            albgs = sical0 - 2.45E-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          else
          end if
! ===================================================================
!         4.  correct for snow cover
! ===================================================================
          if ( scv1d(n,np).gt.0.0D0 ) then
!           **********            snow albedo depends on  snow-age,
!           zenith angle, **********            and thickness of snow
 
!           **********            zenith angle set in zenitm
!           **********            snow albedoes for visible and ir
!           solar rad **********            visible albedo depends on
!           snow age **********            age gives reduction of
!           visible rad snow albedo **********              due to age
            cons = 0.2D0
            conn = 0.5D0
            age = (1.D0-1.D0/(1.D0+sag1d(n,np)))
!           **********            sl helps control albedo zenith
!           dependence
            sl = 2.0D0
            sli = 1.D0/sl
            sl2 = 2.D0*sl
!           **********            snal0= new snow albedo for vis rad,
!           sol zen le 6 **********            snal1= new snow albedo
!           for long-wave rad
            dfalbs = snal0*(1.D0-cons*age)
!           **********            czf corrects albedo of new snow for
!           solar zenith
            cf1 = ((1.D0+sli)/(1.D0+sl2*czen(np))-sli)
            cff = dmax1(cf1,0.D0)
            czf = 0.4D0*cff*(1.D0-dfalbs)
            dralbs = dfalbs + czf
            dfalbl = snal1*(1.D0-conn*age)
            czf = 0.4D0*cff*(1.D0-dfalbl)
            dralbl = dfalbl + czf
 
            if ( veg1d(n,np).gt.0.001D0 ) then
!             **********            effective albedo over vegetation
!             with snow
              albl = (1.D0-wt(n,np))*albl + dralbl*wt(n,np)
              albld = (1.D0-wt(n,np))*albld + dfalbl*wt(n,np)
              albs = (1.D0-wt(n,np))*albs + dralbs*wt(n,np)
              albsd = (1.D0-wt(n,np))*albsd + dfalbs*wt(n,np)
            end if
 
!----------------------------------------------------------------------
!           4.1  compute albedo for snow on bare ground
!----------------------------------------------------------------------
            albgs = (1.D0-scvk(n,np))*albgs + dralbs*scvk(n,np)
            albgl = (1.D0-scvk(n,np))*albgl + dralbl*scvk(n,np)
            albgsd = (1.D0-scvk(n,np))*albgsd + dfalbs*scvk(n,np)
            albgld = (1.D0-scvk(n,np))*albgld + dfalbl*scvk(n,np)
          end if
 
!=====================================================================
!         5.  albedo over open ocean
!=====================================================================
          if ( ldoc1d(n,np).eq.0.D0 ) then
!           *********   ocean albedo depends on zenith angle
            if ( czeta.ge.0.0D0 ) then
!             **********   albedo independent of wavelength
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            end if
          end if
 
!
!         ***************not part of albedo in the ccm ****************
!
          aldirs_s(n) = (1.D0-veg1d(n,np))*albgs + veg1d(n,np)*albs
          aldirl_s(n) = (1.D0-veg1d(n,np))*albgl + veg1d(n,np)*albl
          aldifs_s(n) = (1.D0-veg1d(n,np))*albgsd + veg1d(n,np)*albsd
          aldifl_s(n) = (1.D0-veg1d(n,np))*albgld + veg1d(n,np)*albld
        end do
        albvs(np) = albvs_s(1)
        albvl(np) = albvl_s(1)
        aldirs(np) = aldirs_s(1)
        aldirl(np) = aldirl_s(1)
        aldifs(np) = aldifs_s(1)
        aldifl(np) = aldifl_s(1)
        if ( iemiss.eq.1 ) emiss1d(np) = emiss2d(1,np,j)
        do n = 2 , nnsg
          albvs(np) = albvs(np) + albvs_s(n)
          albvl(np) = albvl(np) + albvl_s(n)
          aldirs(np) = aldirs(np) + aldirs_s(n)
          aldirl(np) = aldirl(np) + aldirl_s(n)
          aldifs(np) = aldifs(np) + aldifs_s(n)
          aldifl(np) = aldifl(np) + aldifl_s(n)
          if ( iemiss.eq.1 ) emiss1d(np) = emiss1d(np) + emiss2d(n,np,j)
        end do
        albvs(np) = albvs(np)/dble(nnsg)
        albvl(np) = albvl(np)/dble(nnsg)
        aldirs(np) = aldirs(np)/dble(nnsg)
        aldirl(np) = aldirl(np)/dble(nnsg)
        aldifs(np) = aldifs(np)/dble(nnsg)
        aldifl(np) = aldifl(np)/dble(nnsg)
        if ( iemiss.eq.1 ) emiss1d(np) = emiss1d(np)/dble(nnsg)
 
!       ******   fsw1d(np),sabveg(np),solis(np) computed in colrad
 
      end do
 
      end subroutine albedov
