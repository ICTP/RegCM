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

      subroutine albedoclm(j)
      use mod_regcm_param
      use mod_bats
      implicit none
!
! Dummy arguments
!
      integer :: j
!
! Local variables
!
      real(8) :: albg , albgl , albgld , albgs , albgsd , albl , albld ,&
               & albs , albsd , albzn , alwet , czeta , fsol1 , fsol2 , &
               & sfac , sical0 , sical1 , snal0 , snal1 , wet , x
      real(8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 & aldifs_s , aldirl_s , aldirs_s
      real(8) :: fseas
      integer :: kolour , n , i
!
! *******************************************************************
! albedo calculates fragmented albedos (direct and diffuse) in   
! wavelength regions split at 0.7um.                         
!
! ccm hands albedos to radiation package which computes         
! fsw1d(np) = net solar absorbed over full grid square   
! sabveg(np) = vegetation absorbed (full solar spectrum)  
! solis(np) = shortwave  solar incident                  
!
! here these are calculated at the end of albedo - they use
! only direct albedos for now                    
!
! in both versions :  lftemp uses sabveg tgrund uses sabveg fsw1d(np)
! to get ground absorbed solar photosynthesis uses solis - see
! subrouts stomat and co2 (carbon)  
!
! ******************************************************************** 
! for sea, sea-ice veg albedos are not set                       
! these albedos are not treated as arrays here   
!
! (depuv/10.0)= the ratio of upper soil layer to total         
! root depth; used to compute "wet" for soil albedo
!*********************************************************************
 
!     =================================================================
!l    1. set initial parameters
!     =================================================================

      fseas(x) = dmax1(0.D0,1.D0-0.0016D0*dmax1(298.D0-x,0.D0)**2)
!
!l    1.1 constants
!     **********    solar flux partitioned at wavelength of 0.7micr
      fsol1 = 0.5
      fsol2 = 0.5
!     **********    short and long wave albedo for new snow
      snal0 = 0.95
      snal1 = 0.65
!     **********    short and long wave albedo for sea ice
      sical0 = 0.6
      sical1 = 0.4
 
!     in depth, wt is frac of grid square covered by
!     snow; depends on average snow depth, vegetation, etc.

      call depth
 
!l    1.3  set default vegetation and albedo
!     ( do loop 50 in ccm not used here )
 
      do i = 2 , iym1
        czen(i) = dmax1(coszrs(i),0.D0)
        czeta = czen(i)
        do n = 1 , nnsg
          albgs = 0.0
          albgl = 0.0
          albgsd = 0.0
          albgld = 0.0
          albs = 0.0
          albl = 0.0
          albsd = 0.0
          albld = 0.0
          albvs_s(n) = 0.0
          albvl_s(n) = 0.0
 
!================================================================
!l        2.   get albedo over land
!================================================================
! can't use pointer "nalbk" here because not set - use nldock
! instead tgb1d(i) used instead of tbelow
!
!***      abt modified below
!***      Since CLM uses a value of 1 albedo when cosz <= 0 there is a
!***      potential for a floating point exception in radcsw.F routine
!***      when calculating sabveg and other variables such as fluxdn
!***      and fluxup. In order to by-pass this error we use below
!***      calculation for albedo during night conditions or whenever

!original albedo = 1 if(ldoc1d(n,i).gt.0.1.and.sice1d(n,i).eq.0.) then

          if ( ocld2d(n,i,j)==1. .and. (aldirs2d(i,j)==1 .or.           &
             & aldifs2d(i,j)==1. .or. aldirl2d(i,j)==1 .or.             &
             & aldifl2d(i,j)==1. ) ) then
 
            sfac = 1. - fseas(tgb1d(n,i))
 
! **********  ccm tests here on land mask for veg and soils
! data *** reduces albedo at low temps !!!!!should respond to
! moisture too the following card inactivated (commented out)
! (pat, 27 oct 86)
! veg1d(i)=vegc(lveg(i))-seasf(lveg(i))*sfac abt added
! lveg statement because lveg typically calculated in interf
! which CLM does not use

            lveg(n,i) = nint(veg2d1(n,i,j))
            albs = albvgs(lveg(n,i))
            albl = albvgl(lveg(n,i))
 
!----------------------------------------------------------------------
            if ( (lveg(n,i)<12) .or. (lveg(n,i)>15) ) then
 
!l            2.1  bare soil albedos (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,i))
              wet = ssw1d(n,i)/depuv(lveg(n,i))
              alwet = dmax1((11.D0-40.D0*wet),0.D0)*0.01
              alwet = dmin1(alwet,solour(kolour))
              albg = solour(kolour) + alwet
              albgs = albg
              albgl = 2.*albg
!             higher nir albedos set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl
 
              albzn = 0.85 + 1./(1.+10.*czen(i))
 
!             leafless hardwood canopy: no or inverse zen dep
              if ( lveg(n,i)==5 .and. sfac<0.1 ) albzn = 1.
!             multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
 
!             albedo over vegetation after zenith angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
 
            else if ( lveg(n,i)==12 ) then
 
!l            2.2   permanent ice sheet
              albgs = 0.8
              albgsd = 0.8
              albgl = 0.55
              albgld = 0.55
            else
 
!l            2.3  inland water, swamps, rice paddies etc.
              albg = 0.05/(czeta+0.15)
              albgs = albg
              albgsd = albg
              albgl = albg
              albgld = albg
            end if
          end if
 
!=====================================================================
!l        5.  albedo over open ocean
!=====================================================================
 
 
          if ( ocld2d(n,i,j)==0. ) then
!           *********   ocean albedo depends on zenith angle
            albg = 0.05/(czeta+0.15)
            albgs = albg
            albgl = albg
            albgsd = 0.08
            albgld = 0.08
          end if
          aldirs_s(n)=(1.-veg1d(n,i))*albgs +veg1d(n,i)*albs
          aldirl_s(n)=(1.-veg1d(n,i))*albgl +veg1d(n,i)*albl
          aldifs_s(n)=(1.-veg1d(n,i))*albgsd+veg1d(n,i)*albsd
          aldifl_s(n)=(1.-veg1d(n,i))*albgld+veg1d(n,i)*albld
        end do
        albvs(i)  = albvs_s(1)
        albvl(i)  = albvl_s(1)
        aldirs(i) = aldirs_s(1)
        aldirl(i) = aldirl_s(1)
        aldifs(i) = aldifs_s(1)
        aldifl(i) = aldifl_s(1)

        do n = 1 , nnsg
          albvs(i)  = albvs(i) +albvs_s(n)
          albvl(i)  = albvl(i) +albvl_s(n)
          aldirs(i) = aldirs(i)+aldirs_s(n)
          aldirl(i) = aldirl(i)+aldirl_s(n)
          aldifs(i) = aldifs(i)+aldifs_s(n)
          aldifl(i) = aldifl(i)+aldifl_s(n)
        end do
        albvs(i)  = albvs(i)/dble(nnsg)
        albvl(i)  = albvl(i)/dble(nnsg)
        aldirs(i) = aldirs(i)/dble(nnsg)
        aldirl(i) = aldirl(i)/dble(nnsg)
        aldifs(i) = aldifs(i)/dble(nnsg)
        aldifl(i) = aldifl(i)/dble(nnsg)
      end do
 
      end subroutine albedoclm
