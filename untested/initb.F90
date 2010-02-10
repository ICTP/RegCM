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
 
      subroutine initb

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ***  provides constant initial fields to boundary subroutine
!     ***    soil textures are set in soilbc
!     ***    soil   colors are set in albedo
!
!     ***   units are si
!
      use mod_regcm_param
      use mod_param2
      use mod_bats , only : pptnc , pptc , veg2d , veg2d1 , ocld2d ,    &
                  & tg2d , tgb2d , taf2d , tlef2d , swt2d , srw2d ,     &
                  & ssw2d , dew2d , sag2d , scv2d , sice2d , gwet2d ,   &
                  & sena2d , evpa2d , rnos2d , rno2d , ocld2d , ircp2d ,&
                  & fsw2d , flw2d , sabv2d , sol2d , fswa2d , flwa2d ,  &
                  & prca2d , prnca2d , svga2d , sina2d , satbrt1 ,      &
                  & iexsol , xmopor , deptv , deprv , depuv
      use mod_main
      implicit none
! 
! Local variables
!
      integer :: ill , is , itex , jll , k , nlveg
      real(8) , dimension(20) :: slmo
!
!----------------------------------------------------------------------
!     slmo   : surface moisture availability in fraction of one.
!     data slmo/0.75,0.60,0.75,0.75,0.75,
!     &          0.80,0.60,0.10,0.90,0.90,
!     &          0.10,0.95,0.90,1.00,1.00,
!     &          0.60,0.75,0.85,0.85,0.85 /
      data slmo/0.65 , 0.45 , 0.60 , 0.60 , 0.65 , 0.65 , 0.55 , 0.10 , &
         & 0.90 , 0.80 , 0.20 , 0.90 , 0.90 , 1.00 , 1.00 , 0.50 ,      &
         & 0.50 , 0.65 , 0.60 , 0.60/       !BATS land types
 
!     ****** vgtran transforms mm4 surface types into bats surface types
!MM4  real(kind=8)  vgtran(13)
!MM4  data vgtran /5.,1.,2.,5.,3.,18.,0.,13.,11.,9.,12.,6.,7./
!MM4  integer ist
 
!     ****** typically resp=1.0 kg/m**2/s and changes by<10% in 10 days
!
#ifdef MPP1
      do jll = 1 , jendx
        do ill = 1 , ixm1
          pptnc(ill,jll) = 0.
          pptc(ill,jll) = 0.
!MM4      ist=nint(satbrt(ill,jll))
!MM4      if(ist.le.13)then
!MM4      ist=ist
!MM4      else
!MM4      ist=7
!MM4      end if
!MM4      veg2d(ill,jll)=vgtran(ist)
!eros     note:  when using bats dataset comment line above
          veg2d(ill,jll) = satbrt(ill,jll)
          if ( satbrt(ill,jll).gt.13.9 .and. satbrt(ill,jll).lt.15.1 )  &
             & veg2d(ill,jll) = 0.
        end do
      end do
      do jll = 1 , jendx
        do ill = 1 , ixm1
          do k = 1 , nnsg
            veg2d1(k,ill,jll) = satbrt1(k,ill,jll)
            if ( satbrt1(k,ill,jll).gt.13.9 .and. satbrt1(k,ill,jll)    &
               & .lt.15.1 ) veg2d1(k,ill,jll) = 0.
          end do
        end do
      end do
#else
      do jll = 1 , jx - 1
        do ill = 1 , ixm1
          pptnc(ill,jll) = 0.
          pptc(ill,jll) = 0.
!MM4      ist=nint(satbrt(ill,jll))
!MM4      if(ist.le.13)then
!MM4      ist=ist
!MM4      else
!MM4      ist=7
!MM4      end if
!MM4      veg2d(ill,jll)=vgtran(ist)
!eros     note:  when using bats dataset comment line above
          veg2d(ill,jll) = satbrt(ill,jll)
          if ( satbrt(ill,jll).gt.13.9 .and. satbrt(ill,jll).lt.15.1 )  &
             & veg2d(ill,jll) = 0.
        end do
      end do
      do jll = 1 , jx - 1
        do ill = 1 , ixm1
          do k = 1 , nnsg
            veg2d1(k,ill,jll) = satbrt1(k,ill,jll)
            if ( satbrt1(k,ill,jll).gt.13.9 .and. satbrt1(k,ill,jll)    &
               & .lt.15.1 ) veg2d1(k,ill,jll) = 0.
          end do
        end do
      end do
#endif
 
!     ******  initialize hostetler lake model
#ifdef MPP1
      if ( lakemod.eq.1 ) call initlk(veg2d1,(ixm1)*nsg,jxp*nsg)
      do jll = 1 , jendx
        do ill = 1 , ixm1
          do k = 1 , nnsg
            if ( veg2d1(k,ill,jll).lt.0.5 ) then
              ocld2d(k,ill,jll) = 0.
            else
              ocld2d(k,ill,jll) = 1.
            end if
            nlveg = nint(veg2d1(k,ill,jll))
            if ( nlveg.eq.0 ) then
              nlveg = 15
            else
              nlveg = nlveg
            end if
!sol        itex=int(text2d(k,ill,jll))
            itex = iexsol(nlveg)
            tg2d(k,ill,jll) = tgb(ill,jll)
            tgb2d(k,ill,jll) = tgb(ill,jll)
            taf2d(k,ill,jll) = tgb(ill,jll)
            tlef2d(k,ill,jll) = tgb(ill,jll)
 
!           ******  initialize soil moisture in the 3 layers
            is = nint(satbrt1(k,ill,jll))
            swt2d(k,ill,jll) = deptv(nlveg)*xmopor(itex)*slmo(is)
            srw2d(k,ill,jll) = deprv(nlveg)*xmopor(itex)*slmo(is)
            ssw2d(k,ill,jll) = depuv(nlveg)*xmopor(itex)*slmo(is)
 
            dew2d(k,ill,jll) = 0.
            sag2d(k,ill,jll) = 0.
            scv2d(k,ill,jll) = dmax1(snowc(k,ill,jll),0.D0)
            sice2d(k,ill,jll) = 0.
            gwet2d(k,ill,jll) = 0.5
            sena2d(k,ill,jll) = 0.
            evpa2d(k,ill,jll) = 0.
            rnos2d(k,ill,jll) = 0.
            rno2d(k,ill,jll) = 0.
            if ( sice2d(k,ill,jll).gt.0. ) ocld2d(k,ill,jll) = 2.
            ircp2d(k,ill,jll) = 0.
          end do
        end do
      end do
#else
      if ( lakemod.eq.1 ) call initlk(veg2d1,(ixm1)*nsg,(jxm1)*nsg)

      do jll = 1 , jx - 1
        do ill = 1 , ixm1
          do k = 1 , nnsg
            if ( veg2d1(k,ill,jll).lt.0.5 ) then
              ocld2d(k,ill,jll) = 0.
            else
              ocld2d(k,ill,jll) = 1.
            end if
            nlveg = nint(veg2d1(k,ill,jll))
            if ( nlveg.eq.0 ) then
              nlveg = 15
            else
              nlveg = nlveg
            end if
!sol        itex=int(text2d(k,ill,jll))
            itex = iexsol(nlveg)
            tg2d(k,ill,jll) = tgb(ill,jll)
            tgb2d(k,ill,jll) = tgb(ill,jll)
            taf2d(k,ill,jll) = tgb(ill,jll)
            tlef2d(k,ill,jll) = tgb(ill,jll)

!           ******  initialize soil moisture in the 3 layers
            is = nint(satbrt1(k,ill,jll))
            swt2d(k,ill,jll) = deptv(nlveg)*xmopor(itex)*slmo(is)
            srw2d(k,ill,jll) = deprv(nlveg)*xmopor(itex)*slmo(is)
            ssw2d(k,ill,jll) = depuv(nlveg)*xmopor(itex)*slmo(is)

            dew2d(k,ill,jll) = 0.
            sag2d(k,ill,jll) = 0.
            scv2d(k,ill,jll) = dmax1(snowc(k,ill,jll),0.D0)
            sice2d(k,ill,jll) = 0.
            gwet2d(k,ill,jll) = 0.5
            sena2d(k,ill,jll) = 0.
            evpa2d(k,ill,jll) = 0.
            rnos2d(k,ill,jll) = 0.
            rno2d(k,ill,jll) = 0.
            if ( sice2d(k,ill,jll).gt.0. ) ocld2d(k,ill,jll) = 2.
            ircp2d(k,ill,jll) = 0.
          end do
        end do
      end do
#endif
 
#ifdef MPP1
      do jll = 1 , jendx
        do ill = 1 , ixm1
          fsw2d(ill,jll) = 0.
          flw2d(ill,jll) = 0.
          sabv2d(ill,jll) = 0.
          sol2d(ill,jll) = 0.
          fswa2d(ill,jll) = 0.
          flwa2d(ill,jll) = 0.
          prca2d(ill,jll) = 0.
          prnca2d(ill,jll) = 0.
          svga2d(ill,jll) = 0.
          sina2d(ill,jll) = 0.
        end do
      end do
#else
      do jll = 1 , jx - 1
        do ill = 1 , ixm1
          fsw2d(ill,jll) = 0.
          flw2d(ill,jll) = 0.
          sabv2d(ill,jll) = 0.
          sol2d(ill,jll) = 0.
          fswa2d(ill,jll) = 0.
          flwa2d(ill,jll) = 0.
          prca2d(ill,jll) = 0.
          prnca2d(ill,jll) = 0.
          svga2d(ill,jll) = 0.
          sina2d(ill,jll) = 0.
        end do
      end do
#endif
 
      end subroutine initb
