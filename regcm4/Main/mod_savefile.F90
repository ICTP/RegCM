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

      module mod_savefile

        use mod_dynparam
        use mod_message
        use mod_runparams
        use mod_bats
        use mod_pmoist
        use mod_main
        use mod_mainchem
        use mod_bdycod
        use mod_rad
        use mod_trachem
        use mod_date
        use mod_radiation
        use mod_lake
#ifndef BAND
        use mod_diagnosis
#endif
#ifdef MPP1
        use mod_mppio
#ifdef CLM
        use mod_clm
#endif
#endif
        private

        public :: read_savefile_part1 , read_savefile_part2
        public :: write_savefile

        integer :: isavlast
        integer :: iutrst
        logical :: lrp1

        data isavlast /-1/
        data iutrst   /-1/
        data lrp1 /.false./

        contains

        subroutine read_savefile_part1(idate)
          implicit none
          integer , intent(in) :: idate
          character(256) :: ffin
          character(16) :: fbname
          logical :: existing

          iutrst = 14
          write (fbname, '(a,i10)') 'SAV.', idate
          ffin = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
          inquire (file=ffin,exist=existing)
          if ( .not.existing ) then
            write (aline,*) 'The following SAV File does not exist: ' , &
                &            trim(ffin), ' please check location'
            call say
            call fatal(__FILE__,__LINE__, 'SAV FILE NOT FOUND')
          else
            open (iutrst,file=ffin,form='unformatted',status='old')
          end if

          read (iutrst) mdate0
          jyear0 = mdate0/1000000
          read (iutrst) ktau , xtime , ldatez , lyear , lmonth , lday ,&
                     & lhour , ntime
          jyear = lyear
          jyearr = jyear
          ktaur = ktau
#ifdef MPP1
          if ( ehso4 ) then
            read (iutrst) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io , so0_io
          else
            read (iutrst) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io
          end if
          read (iutrst) ua_io
          read (iutrst) ub_io
          read (iutrst) va_io
          read (iutrst) vb_io
          read (iutrst) ta_io
          read (iutrst) tb_io
          read (iutrst) qva_io
          read (iutrst) qvb_io
          read (iutrst) qca_io
          read (iutrst) qcb_io
          read (iutrst) psa_io , psb_io , satbrt_io , satbrt1_io , f_io
          read (iutrst) ht_io , ht1_io , msfx_io , msfd_io , xlat_io , &
                     & xlong_io
          read (iutrst) tga_io , tgb_io , rainc_io , rainnc_io
          if ( icup.eq.1 ) then
            read (iutrst) rsheat_io , rswat_io
          else if ( icup.eq.3 ) then
            read (iutrst) tbase_io , cldefi_io
          else if ( icup.eq.4 ) then
            read (iutrst) cbmf2d_io
          else
          end if
          read (iutrst) hfx_io , qfx_io , snowc_io , uvdrag_io
#ifndef BAND
          if (debug_level > 2) call restdiag(iutrst)
#endif
          read (iutrst) absnxt_io , abstot_io , emstot_io
          if ( ipptls.eq.1 ) read (iutrst) fcc_io
#ifdef CLM
          read (iutrst) sols2d_io
          read (iutrst) soll2d_io
          read (iutrst) solsd2d_io
          read (iutrst) solld2d_io
          read (iutrst) flwd2d_io
          read (iutrst) aldirs2d_io
          read (iutrst) aldirl2d_io
          read (iutrst) aldifs2d_io
          read (iutrst) aldifl2d_io
          read (iutrst) coszrs2d_io
          read (iutrst) ocld2d_io
          read (iutrst) heatrt_io
          read (iutrst) o3prof_io
          read (iutrst) tgbb_io
          read (iutrst) flw2d_io
          read (iutrst) swt2d_io
          read (iutrst) sinc2d_io
          read (iutrst) fsw2d_io
          read (iutrst) taf2d_io
#else
          read (iutrst) sol2d_io , solvd2d_io , solvs2d_io , flw2d_io ,&
                     & flwd2d_io , fsw2d_io , sabv2d_io , sinc2d_io
          read (iutrst) taf2d_io , tlef2d_io , tgbb_io , ssw2d_io ,    &
                     & srw2d_io , tg2d_io , tgb2d_io , swt2d_io ,     &
                     & scv2d_io , gwet2d_io , veg2d_io , veg2d1_io ,  &
                     & sag2d_io , sice2d_io , dew2d_io , ircp2d_io ,  &
                     & text2d_io , col2d_io , ocld2d_io , heatrt_io , &
                     & o3prof_io
#endif
          read (iutrst) pptnc_io , pptc_io , prca2d_io , prnca2d_io
          if ( iocnflx.eq.2 ) read (iutrst) zpbl_io
          if ( ichem.eq.1 ) then
            read (iutrst) chia_io
            read (iutrst) chib_io
!           cumul removal terms (3d, 2d)
            read (iutrst) remlsc_io
            read (iutrst) remcvc_io
            read (iutrst) remdrd_io
            read (iutrst) ssw2da_io
            read (iutrst) sdeltk2d_io
            read (iutrst) sdelqk2d_io
            read (iutrst) sfracv2d_io
            read (iutrst) sfracb2d_io
            read (iutrst) sfracs2d_io
            read (iutrst) svegfrac2d_io
!           cumul ad, dif, emis terms ( scalar)
#ifndef BAND
            if (debug_level > 2) call restchemdiag(iutrst)
#endif
          end if
#else
          if ( ehso4 ) then
            read (iutrst) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
          else
            read (iutrst) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
          end if
          read (iutrst) ua
          read (iutrst) ub
          read (iutrst) va
          read (iutrst) vb
          read (iutrst) ta
          read (iutrst) tb
          read (iutrst) qva
          read (iutrst) qvb
          read (iutrst) qca
          read (iutrst) qcb
          read (iutrst) psa , psb , satbrt , satbrt1 , f
          read (iutrst) ht , ht1 , msfx , msfd , xlat , xlong
          read (iutrst) tga , tgb , rainc , rainnc
          if ( icup.eq.1 ) then
            read (iutrst) rsheat , rswat
          else if ( icup.eq.3 ) then
            read (iutrst) tbase , cldefi
          else if ( icup.eq.4 ) then
            read (iutrst) cbmf2d
          else
          end if
          read (iutrst) hfx , qfx , snowc , uvdrag
#ifndef BAND
          if (debug_level > 2) call restdiag(iutrst)
#endif
          read (iutrst) absnxt , abstot , emstot
          if ( ipptls.eq.1 ) read (iutrst) fcc
          read (iutrst) sol2d , solvd2d , solvs2d , flw2d , flwd2d ,    &
                     & fsw2d , sabv2d , sinc2d
          read (iutrst) taf2d , tlef2d , tgbb , ssw2d , srw2d , tg2d ,  &
                     & tgb2d , swt2d , scv2d , gwet2d , veg2d , veg2d1 ,&
                     & sag2d , sice2d , dew2d , ircp2d , text2d ,       &
                     & col2d , ocld2d , heatrt , o3prof
          read (iutrst) pptnc , pptc , prca2d , prnca2d
          if ( iocnflx.eq.2 ) read (iutrst) zpbl
          if ( ichem.eq.1 ) then
            read (iutrst) chia
            read (iutrst) chib
!           cumul removal terms (3d, 2d)
            read (iutrst) remlsc
            read (iutrst) remcvc
            read (iutrst) remdrd
            read (iutrst) ssw2da
            read (iutrst) sdeltk2d
            read (iutrst) sdelqk2d
            read (iutrst) sfracv2d
            read (iutrst) sfracb2d
            read (iutrst) sfracs2d
            read (iutrst) svegfrac2d
!           cumul ad, dif, emis terms ( scalar)
#ifndef BAND
            if (debug_level > 2) call restchemdiag(iutrst)
#endif
          end if
#endif
          if ( lakemod.eq.1 ) then
            call lakesavread(iutrst)
          end if
          lrp1 = .true.
        end subroutine read_savefile_part1

        subroutine read_savefile_part2
          implicit none

          if (.not. lrp1) then
            write (6,*) 'Reading part2 before part1'
            call fatal(__FILE__,__LINE__, 'SAV FILE ERROR')
          end if

#ifdef MPP1
          read (iutrst) dstor_io
          read (iutrst) hstor_io
#ifndef BAND
          read (iutrst) uj1 , uj2 , ujlx , ujl
#endif
          read (iutrst) ui1_io , ui2_io , uilx_io , uil_io
#ifndef BAND
          read (iutrst) vj1 , vj2 , vjlx , vjl
#endif
          read (iutrst) vi1_io , vi2_io , vilx_io , vil_io
#else
          read (iutrst) dstor
          read (iutrst) hstor
#ifndef BAND
          read (iutrst) uj1 , uj2 , ujlx , ujl
#endif
          read (iutrst) ui1 , ui2 , uilx , uil
#ifndef BAND
          read (iutrst) vj1 , vj2 , vjlx , vjl
#endif
          read (iutrst) vi1 , vi2 , vilx , vil
#endif
          close(iutrst)
        end subroutine read_savefile_part2

        subroutine write_savefile(idate,ltmp)
          implicit none
          integer , intent(in) :: idate
          logical , intent(in) :: ltmp
          integer , parameter :: iutsav = 52
          character(256) :: ffout
          character(32) :: fbname
          logical :: existing
          if (ltmp) then
            write (fbname, '(a,i10)') 'TMPSAV.', idate
          else
            write (fbname, '(a,i10)') 'SAV.', idate
          end if
          ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
          open (iutsav,file=ffout,form='unformatted',status='replace')
          write (iutsav) mdate0

          inquire (file=ffout,exist=existing)
          if ( .not.existing ) then
            write (aline,*) 'The SAV File cannot be created: ' , &
                &            trim(ffout), ' please check directory'
            call say
            call fatal(__FILE__,__LINE__, 'SAV FILE WRITE ERROR')
          end if

          write (iutsav) ktau , xtime , ldatez , lyear , lmonth , lday ,&
                     & lhour , ntime
#ifdef MPP1
          if ( ehso4 ) then
            write (iutsav) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io , so0_io
          else
            write (iutsav) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io
          end if
          write (iutsav) ua_io
          write (iutsav) ub_io
          write (iutsav) va_io
          write (iutsav) vb_io
          write (iutsav) ta_io
          write (iutsav) tb_io
          write (iutsav) qva_io
          write (iutsav) qvb_io
          write (iutsav) qca_io
          write (iutsav) qcb_io
          write (iutsav) psa_io , psb_io , satbrt_io , satbrt1_io , f_io
          write (iutsav) ht_io , ht1_io , msfx_io , msfd_io , xlat_io , &
                     & xlong_io
          write (iutsav) tga_io , tgb_io , rainc_io , rainnc_io
          if ( icup.eq.1 ) then
            write (iutsav) rsheat_io , rswat_io
          else if ( icup.eq.3 ) then
            write (iutsav) tbase_io , cldefi_io
          else if ( icup.eq.4 ) then
            write (iutsav) cbmf2d_io
          else
          end if
          write (iutsav) hfx_io , qfx_io , snowc_io , uvdrag_io
#ifndef BAND
          if (debug_level > 2) call savediag(iutsav)
#endif
          write (iutsav) absnxt_io , abstot_io , emstot_io
          if ( ipptls.eq.1 ) write (iutsav) fcc_io
#ifdef CLM
          write (iutsav) sols2d_io
          write (iutsav) soll2d_io
          write (iutsav) solsd2d_io
          write (iutsav) solld2d_io
          write (iutsav) flwd2d_io
          write (iutsav) aldirs2d_io
          write (iutsav) aldirl2d_io
          write (iutsav) aldifs2d_io
          write (iutsav) aldifl2d_io
          write (iutsav) coszrs2d_io
          write (iutsav) ocld2d_io
          write (iutsav) heatrt_io
          write (iutsav) o3prof_io
          write (iutsav) tgbb_io
          write (iutsav) flw2d_io
          write (iutsav) swt2d_io
          write (iutsav) sinc2d_io
          write (iutsav) fsw2d_io
          write (iutsav) taf2d_io
#else
          write (iutsav) sol2d_io , solvd2d_io , solvs2d_io , flw2d_io ,&
                     & flwd2d_io , fsw2d_io , sabv2d_io , sinc2d_io
          write (iutsav) taf2d_io , tlef2d_io , tgbb_io , ssw2d_io ,    &
                     & srw2d_io , tg2d_io , tgb2d_io , swt2d_io ,     &
                     & scv2d_io , gwet2d_io , veg2d_io , veg2d1_io ,  &
                     & sag2d_io , sice2d_io , dew2d_io , ircp2d_io ,  &
                     & text2d_io , col2d_io , ocld2d_io , heatrt_io , &
                     & o3prof_io
#endif
          write (iutsav) pptnc_io , pptc_io , prca2d_io , prnca2d_io
          if ( iocnflx.eq.2 ) write (iutsav) zpbl_io
          if ( ichem.eq.1 ) then
            write (iutsav) chia_io
            write (iutsav) chib_io
!           cumul removal terms (3d, 2d)
            write (iutsav) remlsc_io
            write (iutsav) remcvc_io
            write (iutsav) remdrd_io
            write (iutsav) ssw2da_io
            write (iutsav) sdeltk2d_io
            write (iutsav) sdelqk2d_io
            write (iutsav) sfracv2d_io
            write (iutsav) sfracb2d_io
            write (iutsav) sfracs2d_io
            write (iutsav) svegfrac2d_io
!           cumul ad, dif, emis terms ( scalar)
#ifndef BAND
            if (debug_level > 2) call savechemdiag(iutsav)
#endif
          end if
#else
          if ( ehso4 ) then
            write (iutsav) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
          else
            write (iutsav) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
          end if
          write (iutsav) ua
          write (iutsav) ub
          write (iutsav) va
          write (iutsav) vb
          write (iutsav) ta
          write (iutsav) tb
          write (iutsav) qva
          write (iutsav) qvb
          write (iutsav) qca
          write (iutsav) qcb
          write (iutsav) psa , psb , satbrt , satbrt1 , f
          write (iutsav) ht , ht1 , msfx , msfd , xlat , xlong
          write (iutsav) tga , tgb , rainc , rainnc
          if ( icup.eq.1 ) then
            write (iutsav) rsheat , rswat
          else if ( icup.eq.3 ) then
            write (iutsav) tbase , cldefi
          else if ( icup.eq.4 ) then
            write (iutsav) cbmf2d
          else
          end if
          write (iutsav) hfx , qfx , snowc , uvdrag
#ifndef BAND
          if (debug_level > 2) call savediag(iutsav)
#endif
          write (iutsav) absnxt , abstot , emstot
          if ( ipptls.eq.1 ) write (iutsav) fcc
          write (iutsav) sol2d , solvd2d , solvs2d , flw2d , flwd2d ,   &
                     & fsw2d , sabv2d , sinc2d
          write (iutsav) taf2d , tlef2d , tgbb , ssw2d , srw2d , tg2d , &
                     & tgb2d , swt2d , scv2d , gwet2d , veg2d , veg2d1 ,&
                     & sag2d , sice2d , dew2d , ircp2d , text2d ,       &
                     & col2d , ocld2d , heatrt , o3prof
          write (iutsav) pptnc , pptc , prca2d , prnca2d
          if ( iocnflx.eq.2 ) write (iutsav) zpbl
          if ( ichem.eq.1 ) then
            write (iutsav) chia
            write (iutsav) chib
!           cumul removal terms (3d, 2d)
            write (iutsav) remlsc
            write (iutsav) remcvc
            write (iutsav) remdrd
            write (iutsav) ssw2da
            write (iutsav) sdeltk2d
            write (iutsav) sdelqk2d
            write (iutsav) sfracv2d
            write (iutsav) sfracb2d
            write (iutsav) sfracs2d
            write (iutsav) svegfrac2d
!           cumul ad, dif, emis terms ( scalar)
#ifndef BAND
            if (debug_level > 2) call savechemdiag(iutsav)
#endif
          end if
#endif
          if ( lakemod.eq.1 ) then
            call lakesavwrite(iutsav)
          end if
#ifdef MPP1
          write (iutsav) dstor_io
          write (iutsav) hstor_io
#ifndef BAND
          write (iutsav) uj1 , uj2 , ujlx , ujl
#endif
          write (iutsav) ui1_io , ui2_io , uilx_io , uil_io
#ifndef BAND
          write (iutsav) vj1 , vj2 , vjlx , vjl
#endif
          write (iutsav) vi1_io , vi2_io , vilx_io , vil_io
#else
          write (iutsav) dstor
          write (iutsav) hstor
#ifndef BAND
          write (iutsav) uj1 , uj2 , ujlx , ujl
#endif
          write (iutsav) ui1 , ui2 , uilx , uil
#ifndef BAND
          write (iutsav) vj1 , vj2 , vjlx , vjl
#endif
          write (iutsav) vi1 , vi2 , vilx , vil
#endif
          close(iutsav)

          if (ltmp) then
            if (isavlast > 0) then
              write (fbname, '(a,i10)') 'TMPSAV.', isavlast
              ffout = trim(dirglob)//pthsep//trim(domname)// &
                      '_'//trim(fbname)
              call unlink(ffout)
            end if
            isavlast = idate
          end if

        end subroutine write_savefile

      end module mod_savefile
