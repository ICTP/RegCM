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
 
      subroutine outsav(iutl)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output for restart.            c
!                                                                     c
!     iutl : is the output unit number for large-domain variables.    c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : mdate0
      use mod_main
      use mod_bdycod
      use mod_bats
      use mod_date
      use mod_iunits
      use mod_radbuf
      use mod_pmoist
      use mod_mainchem
      use mod_trachem
      use mod_split
      use mod_rad
#ifdef MPP1
      use mod_mppio
#endif
#ifdef DIAG
      use mod_diagnosis
#endif
      implicit none
!
! Dummy arguments
!
      integer :: iutl
      intent (in) iutl
!
! Local variables
!
      integer :: depth , freeze , ilake , j , jlake , n
      real(8) :: eta , hi , hii , hs
      real(8) , dimension(400) :: tlake
!
!-----output large-domain variables:
!
      write (iutl) mdate0
      write (iutl) ktau , xtime , ldatez , lyear , lmonth , lday ,      &
                 & lhour , ntime
#ifdef MPP1
      if ( ehso4 ) then
        write (iutl) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io ,       &
                   & ts0_io , so0_io
      else
        write (iutl) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , ts0_io
      end if
      write (iutl) ua_io
      write (iutl) ub_io
      write (iutl) va_io
      write (iutl) vb_io
      write (iutl) ta_io
      write (iutl) tb_io
      write (iutl) qva_io
      write (iutl) qvb_io
      write (iutl) qca_io
      write (iutl) qcb_io
      write (iutl) psa_io , psb_io , satbrt_io , satbrt1_io , f_io
      write (iutl) ht_io , ht1_io , msfx_io , msfd_io , xlat_io ,       &
                 & xlong_io
      write (iutl) tga_io , tgb_io , rainc_io , rainnc_io
      if ( icup.eq.1 ) then
        write (iutl) rsheat_io , rswat_io
      else if ( icup.eq.3 ) then
        write (iutl) tbase_io , cldefi_io
      else if ( icup.eq.4 ) then
        write (iutl) cbmf2d_io
      else
      end if
      write (iutl) hfx_io , qfx_io , snowc_io , uvdrag_io
#ifdef DIAG
      write (iutl) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
      write (iutl) absnxt_io , abstot_io , emstot_io
      if ( ipptls.eq.1 ) write (iutl) fcc_io
      write (iutl) sol2d_io , solvd2d_io , solvs2d_io , flw2d_io ,      &
                 & flwd2d_io , fsw2d_io , sabv2d_io , sinc2d_io
      write (iutl) taf2d_io , tlef2d_io , tgbb_io , ssw2d_io ,          &
                 & srw2d_io , tg2d_io , tgb2d_io , swt2d_io , scv2d_io ,&
                 & gwet2d_io , veg2d_io , veg2d1_io , sag2d_io ,        &
                 & sice2d_io , dew2d_io , ircp2d_io , text2d_io ,       &
                 & col2d_io , ocld2d_io , heatrt_io , o3prof_io
      write (iutl) pptnc_io , pptc_io , prca2d_io , prnca2d_io
      if ( iocnflx.eq.2 ) write (iutl) zpbl_io
      if ( ichem.eq.1 ) then
        write (iutl) chia_io
        write (iutl) chib_io
        write (iutl) remlsc_io
        write (iutl) remcvc_io
        write (iutl) remdrd_io
        write (iutl) ssw2da_io
        write (iutl) sdeltk2d_io
        write (iutl) sdelqk2d_io
        write (iutl) sfracv2d_io
        write (iutl) sfracb2d_io
        write (iutl) sfracs2d_io
        write (iutl) svegfrac2d_io
#ifdef DIAG
        write (iutl) tchiad
        write (iutl) tchitb
        write (iutl) tchie
#endif
      end if
#else
      if ( ehso4 ) then
        write (iutl) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
      else
        write (iutl) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
      end if
      write (iutl) ua
      write (iutl) ub
      write (iutl) va
      write (iutl) vb
      write (iutl) ta
      write (iutl) tb
      write (iutl) qva
      write (iutl) qvb
      write (iutl) qca
      write (iutl) qcb
      write (iutl) psa , psb , satbrt , satbrt1 , f
      write (iutl) ht , ht1 , msfx , msfd , xlat , xlong
      write (iutl) tga , tgb , rainc , rainnc
      if ( icup.eq.1 ) then
        write (iutl) rsheat , rswat
      else if ( icup.eq.3 ) then
        write (iutl) tbase , cldefi
      else if ( icup.eq.4 ) then
        write (iutl) cbmf2d
      else
      end if
      write (iutl) hfx , qfx , snowc , uvdrag
#ifdef DIAG
      write (iutl) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
      write (iutl) absnxt , abstot , emstot
      if ( ipptls.eq.1 ) write (iutl) fcc
      write (iutl) sol2d , solvd2d , solvs2d , flw2d , flwd2d , fsw2d , &
                 & sabv2d , sinc2d
      write (iutl) taf2d , tlef2d , tgbb , ssw2d , srw2d , tg2d ,       &
                 & tgb2d , swt2d , scv2d , gwet2d , veg2d , veg2d1 ,    &
                 & sag2d , sice2d , dew2d , ircp2d , text2d , col2d ,   &
                 & ocld2d , heatrt , o3prof
      write (iutl) pptnc , pptc , prca2d , prnca2d
      if ( iocnflx.eq.2 ) write (iutl) zpbl

!chem2
      if ( ichem.eq.1 ) then
!       pronostic concentrartions
        write (iutl) chia
        write (iutl) chib
!       cumul removal terms (3d, 2d)
        write (iutl) remlsc
        write (iutl) remcvc
        write (iutl) remdrd
        write (iutl) ssw2da
        write (iutl) sdeltk2d
        write (iutl) sdelqk2d
        write (iutl) sfracv2d
        write (iutl) sfracb2d
        write (iutl) sfracs2d
        write (iutl) svegfrac2d
!       cumul ad, dif, emis terms ( scalar)
#ifdef  DIAG
        write (iutl) tchiad
        write (iutl) tchitb
        write (iutl) tchie
#endif
      end if
!chem2_
#endif
 
!---lake model
!
      if ( lakemod.eq.1 ) then
        write (iutl) numpts
        print * , 'writing lake model restart file. numpts = ' , numpts
        do n = 1 , numpts
          read (iin) ilake , jlake , depth , freeze , hi , hii , hs ,   &
                   & eta , (tlake(j),j=1,depth)
!         print *,' in outsav ilake jlake = ',ilake,jlake
 
          write (iutl) ilake , jlake , depth , freeze , hi , hii , hs , &
                     & eta , (tlake(j),j=1,depth)
        end do
        rewind (iin)
      end if
!
#ifdef MPP1
      write (iutl) dstor_io
      write (iutl) hstor_io
      write (iutl) uj1 , uj2 , ujlx , ujl
      write (iutl) ui1_io , ui2_io , uilx_io , uil_io
      write (iutl) vj1 , vj2 , vjlx , vjl
      write (iutl) vi1_io , vi2_io , vilx_io , vil_io
#else
      write (iutl) dstor
      write (iutl) hstor
      write (iutl) uj1 , uj2 , ujlx , ujl
      write (iutl) ui1 , ui2 , uilx , uil
      write (iutl) vj1 , vj2 , vjlx , vjl
      write (iutl) vi1 , vi2 , vilx , vil
#endif
      end subroutine outsav
