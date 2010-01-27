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

      subroutine albedov(j, iemiss)

      implicit none

C
C  Albedo calculates fragmented albedos (direct and diffuse) in
C    wavelength regions split at 0.7um.
C
C  CM hands albedos to radiation package which computes
C        fsw1d(np) = net solar absorbed over full grid square
C       sabveg(np) = vegetation absorbed (full solar spectrum)
C        solis(np) = shortwave  solar incident
C
C  Here these are calculated at the end of albedo - they use only
C               direct albedos for now
C
C   in both versions :  lftemp uses sabveg
C                       tgrund uses sabveg & fsw1d(np) to get
C                                           ground absorbed solar
C                       photosynthesis uses solis - see subrouts
C                           stomat and co2 (carbon)
C
C  For sea, sea-ice veg albedos are not set
C                these albedos are not treated as arrays here
C
C  (depuv/10.0)= the ratio of upper soil layer to total
C                root depth; used to compute "wet" for soil albedo
C     solour        =  soil albedo for different coloured soils
C
#     include "parameter.inc"

#     include "include/param1.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/bats2.cb"
#     include "include/bats.cb"

      integer j, iemiss

      real(kind=8)  solour(8)
!     data    solour/ .12, .11, .10, .09, .08, .07, .06, .05 /
      data    solour/ .16, .11, .10, .09, .08, .07, .06, .05 /
      real(kind=8)  fsol1,fsol2,snal0,snal1,sical0,sical1
      real(kind=8)  czeta,albgs,albgl,albgsd,albgld,albs,albl,albsd
      real(kind=8)  albld,sfac,wet,alwet,albg,albzn,tdiffs,tdiff
      real(kind=8)  cons,conn,age,sl,sli,sl2,dfalbs,cf1,cff,czf
      real(kind=8)  dralbs,dfalbl,dralbl
      real(kind=8)  albvs_s(NNSG),aldirs_s(NNSG),aldifs_s(NNSG)
      real(kind=8)  albvl_s(NNSG),aldirl_s(NNSG),aldifl_s(NNSG)
      integer n, np, kolour

C  =================================================================
C  1. set initial parameters
C  =================================================================
      real(kind=8)  x, fseas
      fseas(x) = dmax1(0.d0,1.d0-0.0016d0*dmax1(298.d0-x,0.d0)**2)
C
C        1.1 constants
C
C Solar flux partitioned at wavelength of 0.7micr
      fsol1=0.5
      fsol2=0.5
C Short and long wave albedo for new snow
      snal0=0.95
      snal1=0.65
C Short and long wave albedo for sea ice
      sical0=0.6
      sical1=0.4

C
C In depth, wt is frac of grid square covered by snow;
C   depends on average snow depth, vegetation, etc.
C
      npts = nbmax

      call depth

C 1.2 set pointers
C ***************************************************
C *    set n"x"k params here  in ccm but not needed *
C ***************************************************
C
C 1.3  set default vegetation and albedo
C    do loop 50 in ccm not used here )

      do np=np1,npts
        czen(np) = dmax1( coszrs(np), 0.d0 )
        czeta    = czen(np)
        do n=1,NNSG
          albgs  = 0.0
          albgl  = 0.0
          albgsd = 0.0
          albgld = 0.0
          albs   = 0.0
          albl   = 0.0
          albsd  = 0.0
          albld  = 0.0

          albvs_s(n) = 0.0
          albvl_s(n) = 0.0

C================================================================
C        2.   get albedo over land
C================================================================
C  can't use pointer "nalbk" here because not set - use nldock instead
C                            tgb1d(np) used instead of tbelow
C
         if (ldoc1d(n,np).gt.0.1.and.sice1d(n,np).eq.0.) then
           sfac = 1.-fseas(tgb1d(n,np))

C **********  ccm tests here on land mask for veg and soils data


Cc *** reduces albedo at low temps !!!!!should respond to moisture too
Cc    the following card inactivated (commented out) (pat, 27 oct 86)
Cc       veg1d(np)=vegc(lveg(np))-seasf(lveg(np))*sfac

           albs=albvgs(lveg(n,np))
           albl=albvgl(lveg(n,np))

C----------------------------------------------------------------------
           if ((lveg(n,np).lt.12).or.(lveg(n,np).gt.15)) then

C        2.1  bare soil albedos
C               (soil albedo depends on moisture)
             kolour   = kolsol(lveg(n,np))
             wet=ssw1d(n,np)/depuv(lveg(n,np))
             alwet=dmax1((11.d0-40.d0*wet),0.d0) *0.01
             alwet=dmin1(alwet,solour(kolour))
             albg=solour(kolour)+alwet
C            if((lveg(n,np).eq.8)) albg=0.40      !Laura, cambiato il DESERTO
             albgs=albg
             albgl=2.*albg
C **********            higher nir albedos
C **********              set diffuse albedo
             albgld=albgl
             albgsd=albgs
             albsd=albs
             albld=albl

! Dec. 15   albzn=0.85+1./(1.+10.*czen(np))
! Dec. 12, 2008
             albzn=1.0d0
! Dec. 15, 2008

C **********            leafless hardwood canopy: no or inverse zen dep
            if (lveg(n,np).eq.5.and.sfac.lt.0.1) albzn=1.
C **********            multiply by zenith angle correction
              albs=albs*albzn
              albl=albl*albzn

C **********            albedo over vegetation after zenith angle corr
              albvs_s(n)=albs
              albvl_s(n)=albl

            else if (lveg(n,np).eq.12) then

C        2.2   permanent ice sheet
              albgs=0.8
              albgsd=0.8
              albgl=0.55
              albgld=0.55
            else

C        2.3  inland water, swamps, rice paddies etc.
              albg=0.05/(czeta+0.15)
              albgs=albg
              albgsd=albg
              albgl=albg
              albgld=albg
            end if

          else if (sice1d(n,np).gt.0.) then
C====================================================================
C        3.  get albedo over sea ice
C====================================================================
C **********          albedo depends on wave-length and ts. the ts
C **********          dependence accounts for melt water puddles.
            tdiffs=ts1d(n,np)-c(67)
            tdiff=dmax1(tdiffs,0.d0)
            tdiffs=dmin1(tdiff,20.d0)
            albgl=sical1-1.1e-2*tdiffs
            albgs=sical0-2.45e-2*tdiffs
            albg=fsol1*albgs+fsol2*albgl
            albgsd=albgs
            albgld=albgl
          end if
C ===================================================================
C        4.  correct for snow cover
C ===================================================================
          if (scv1d(n,np).gt.0.0) then
C **********            snow albedo depends on  snow-age, zenith angle,
C **********            and thickness of snow

C **********            zenith angle set in zenitm
C **********            snow albedoes for visible and ir solar rad
C **********            visible albedo depends on snow age
C **********            age gives reduction of visible rad snow albedo
C **********              due to age
            cons=0.2
            conn=0.5
            age=  (1.-1./(1.+sag1d(n,np)))
C **********            sl helps control albedo zenith dependence
            sl=2.0
            sli=1./sl
            sl2=2.*sl
C **********            snal0= new snow albedo for vis rad, sol zen le 6
C **********            snal1= new snow albedo for long-wave rad
            dfalbs=snal0*(1.-cons*age)
C **********            czf corrects albedo of new snow for solar zenith
            cf1 = ( (1.+sli) / (1.+sl2*czen(np)) - sli )
            cff=dmax1(cf1,0.d0)
            czf=0.4*cff*(1.-dfalbs)
            dralbs=dfalbs+czf
            dfalbl=snal1*(1.-conn*age)
            czf=0.4*cff*(1.-dfalbl)
            dralbl=dfalbl+czf

            if (veg1d(n,np).gt.0.001) then
C **********            effective albedo over vegetation with snow
              albl =(1.-wt(n,np))*albl +dralbl*wt(n,np)
              albld=(1.-wt(n,np))*albld+dfalbl*wt(n,np)
              albs =(1.-wt(n,np))*albs +dralbs*wt(n,np)
              albsd=(1.-wt(n,np))*albsd+dfalbs*wt(n,np)
            end if

C----------------------------------------------------------------------
C        4.1  compute albedo for snow on bare ground
C----------------------------------------------------------------------
            albgs  = (1.-scvk(n,np))*albgs +dralbs*scvk(n,np)
            albgl  = (1.-scvk(n,np))*albgl +dralbl*scvk(n,np)
            albgsd = (1.-scvk(n,np))*albgsd+dfalbs*scvk(n,np)
            albgld = (1.-scvk(n,np))*albgld+dfalbl*scvk(n,np)
          end if

C=====================================================================
C        5.  albedo over open ocean
C=====================================================================
          if (ldoc1d(n,np).eq.0.) then
C  *********   ocean albedo depends on zenith angle
            if (czeta.lt.0.0) go to 2
C  **********   albedo independent of wavelength
            albg=0.05/(czeta+0.15)
            albgs=albg
            albgl=albg
            albgsd=0.08
            albgld=0.08
          end if

2         continue

C
C  ***************not part of albedo in the ccm ****************
C
          aldirs_s(n)=(1.-veg1d(n,np))*albgs +veg1d(n,np)*albs
          aldirl_s(n)=(1.-veg1d(n,np))*albgl +veg1d(n,np)*albl
          aldifs_s(n)=(1.-veg1d(n,np))*albgsd+veg1d(n,np)*albsd
          aldifl_s(n)=(1.-veg1d(n,np))*albgld+veg1d(n,np)*albld
        end do
        albvs(np)  = albvs_s(1)
        albvl(np)  = albvl_s(1)
        aldirs(np) = aldirs_s(1)
        aldirl(np) = aldirl_s(1)
        aldifs(np) = aldifs_s(1)
        aldifl(np) = aldifl_s(1)
        if(iemiss.eq.1) emiss1d(np) = emiss2d(1,np,j)
        do n = 2, NNSG
          albvs(np)  = albvs(np) +albvs_s(n)
          albvl(np)  = albvl(np) +albvl_s(n)
          aldirs(np) = aldirs(np)+aldirs_s(n)
          aldirl(np) = aldirl(np)+aldirl_s(n)
          aldifs(np) = aldifs(np)+aldifs_s(n)
          aldifl(np) = aldifl(np)+aldifl_s(n)
          if (iemiss.eq.1) emiss1d(np) = emiss1d(np)+emiss2d(n,np,j)
        end do
        albvs(np)  = albvs(np)/dble(NNSG)
        albvl(np)  = albvl(np)/dble(NNSG)
        aldirs(np) = aldirs(np)/dble(NNSG)
        aldirl(np) = aldirl(np)/dble(NNSG)
        aldifs(np) = aldifs(np)/dble(NNSG)
        aldifl(np) = aldifl(np)/dble(NNSG)
        if (iemiss.eq.1) emiss1d(np)= emiss1d(np)/dble(NNSG)

C ******   fsw1d(np),sabveg(np),solis(np) computed in colrad

      end do

      return
      end subroutine albedov
