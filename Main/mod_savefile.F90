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
        use mod_lake, only : lakesav_i, lakesav_o
        use mod_pmoist
        use mod_main
        use mod_mainchem
        use mod_bdycod
        use mod_rad
        use mod_trachem
        use mod_date
        use mod_radiation
        use mod_cu_bm
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
#ifdef MPP1
          if ( ehso4 ) then
            read (iutrst) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io , so0_io
          else
            read (iutrst) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                       & ts0_io
          end if
          read (iutrst) atm1_io%u
          read (iutrst) atm1_io%v
          read (iutrst) atm1_io%t
          read (iutrst) atm1_io%qv
          read (iutrst) atm1_io%qc
          read (iutrst) atm2_io%u
          read (iutrst) atm2_io%v
          read (iutrst) atm2_io%t
          read (iutrst) atm2_io%qv
          read (iutrst) atm2_io%qc
          read (iutrst) psa_io , psb_io
          read (iutrst) tga_io , tgb_io , rainc_io , rainnc_io
          if ( icup.eq.1 ) then
            read (iutrst) rsheat_io , rswat_io
          end if
          if ( icup.eq.3 ) then
            read (iutrst) tbase_io , cldefi_io
          end if
          if ( icup.eq.4 .or. icup.eq.99 .or. icup.eq.98 ) then
            read (iutrst) cbmf2d_io
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
          end if
          read (iutrst) ssw2da_io
          read (iutrst) sdeltk2d_io
          read (iutrst) sdelqk2d_io
          read (iutrst) sfracv2d_io
          read (iutrst) sfracb2d_io
          read (iutrst) sfracs2d_io
          read (iutrst) svegfrac2d_io
!         cumul ad, dif, emis terms ( scalar)
#ifndef BAND
          if ( ichem.eq.1 ) then
            if (debug_level > 2) call restchemdiag(iutrst)
          end if
#endif
#else
          if ( ehso4 ) then
            read (iutrst) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
          else
            read (iutrst) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
          end if
          read (iutrst) atm1%u
          read (iutrst) atm1%v
          read (iutrst) atm1%t
          read (iutrst) atm1%qv
          read (iutrst) atm1%qc
          read (iutrst) atm2%u
          read (iutrst) atm2%v
          read (iutrst) atm2%t
          read (iutrst) atm2%qv
          read (iutrst) atm2%qc
          read (iutrst) sps1%ps , sps2%ps
          read (iutrst) sts1%tg , sts2%tg , sfsta%rainc , sfsta%rainnc
          if ( icup.eq.1 ) then
            read (iutrst) rsheat , rswat
          end if
          if ( icup.eq.3 ) then
            read (iutrst) tbase , cldefi
          end if
          if ( icup.eq.4 .or. icup.eq.99 .or. icup.eq.98 ) then
            read (iutrst) cbmf2d
          end if
          read (iutrst) sfsta%hfx , sfsta%qfx , snowc , sfsta%uvdrag
#ifndef BAND
          if (debug_level > 2) call restdiag(iutrst)
#endif
          read (iutrst) absnxt , abstot , emstot
          if ( ipptls.eq.1 ) read (iutrst) fcc
          read (iutrst) sol2d , solvd2d , solvs2d , flw2d , flwd2d ,    &
                     & fsw2d , sabv2d , sinc2d
          read (iutrst) taf2d , tlef2d , sfsta%tgbb , ssw2d , srw2d ,   &
                     & tg2d ,  tgb2d , swt2d , scv2d , gwet2d , veg2d , &
                     & veg2d1 , sag2d , sice2d , dew2d , ircp2d ,       &
                     & text2d , col2d , ocld2d , heatrt , o3prof
          read (iutrst) pptnc , pptc , prca2d , prnca2d
          if ( iocnflx.eq.2 ) read (iutrst) sfsta%zpbl
          if ( ichem.eq.1 ) then
            read (iutrst) chia
            read (iutrst) chib
!           cumul removal terms (3d, 2d)
            read (iutrst) remlsc
            read (iutrst) remcvc
            read (iutrst) remdrd
          end if
          read (iutrst) ssw2da
          read (iutrst) sdeltk2d
          read (iutrst) sdelqk2d
          read (iutrst) sfracv2d
          read (iutrst) sfracb2d
          read (iutrst) sfracs2d
          read (iutrst) svegfrac2d
!         cumul ad, dif, emis terms ( scalar)
#ifndef BAND
          if ( ichem.eq.1 ) then
            if (debug_level > 2) call restchemdiag(iutrst)
          end if
#endif
#endif
!------lake model
          if ( lakemod.eq.1 ) then
            call lakesav_i(iutrst)
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
          read (iutrst) spsav%dstor
          read (iutrst) spsav%hstor
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
          write (iutsav) atm1_io%u
          write (iutsav) atm1_io%v
          write (iutsav) atm1_io%t
          write (iutsav) atm1_io%qv
          write (iutsav) atm1_io%qc
          write (iutsav) atm2_io%u
          write (iutsav) atm2_io%v
          write (iutsav) atm2_io%t
          write (iutsav) atm2_io%qv
          write (iutsav) atm2_io%qc
          write (iutsav) psa_io , psb_io
          write (iutsav) tga_io , tgb_io , rainc_io , rainnc_io
          if ( icup.eq.1 ) then
            write (iutsav) rsheat_io , rswat_io
          end if
          if ( icup.eq.3 ) then
            write (iutsav) tbase_io , cldefi_io
          end if
          if ( icup.eq.4 .or. icup.eq.99 .or. icup.eq.98 ) then
            write (iutsav) cbmf2d_io
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
          end if
          write (iutsav) ssw2da_io
          write (iutsav) sdeltk2d_io
          write (iutsav) sdelqk2d_io
          write (iutsav) sfracv2d_io
          write (iutsav) sfracb2d_io
          write (iutsav) sfracs2d_io
          write (iutsav) svegfrac2d_io
!         cumul ad, dif, emis terms ( scalar)
#ifndef BAND
          if ( ichem.eq.1 ) then
            if (debug_level > 2) call savechemdiag(iutsav)
          end if
#endif
#else
          if ( ehso4 ) then
            write (iutsav) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
          else
            write (iutsav) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
          end if
          write (iutsav) atm1%u
          write (iutsav) atm1%v
          write (iutsav) atm1%t
          write (iutsav) atm1%qv
          write (iutsav) atm1%qc
          write (iutsav) atm2%u
          write (iutsav) atm2%v
          write (iutsav) atm2%t
          write (iutsav) atm2%qv
          write (iutsav) atm2%qc
          write (iutsav) sps1%ps , sps2%ps
          write (iutsav) sts1%tg , sts2%tg , sfsta%rainc , sfsta%rainnc
          if ( icup.eq.1 ) then
            write (iutsav) rsheat , rswat
          end if
          if ( icup.eq.3 ) then
            write (iutsav) tbase , cldefi
          end if
          if ( icup.eq.4 .or. icup.eq.99 .or. icup.eq.98 ) then
            write (iutsav) cbmf2d
          end if
          write (iutsav) sfsta%hfx , sfsta%qfx , snowc , sfsta%uvdrag
#ifndef BAND
          if (debug_level > 2) call savediag(iutsav)
#endif
          write (iutsav) absnxt , abstot , emstot
          if ( ipptls.eq.1 ) write (iutsav) fcc
          write (iutsav) sol2d , solvd2d , solvs2d , flw2d , flwd2d ,   &
                     & fsw2d , sabv2d , sinc2d
          write (iutsav) taf2d , tlef2d , sfsta%tgbb , ssw2d , srw2d ,  &
                     & tg2d , tgb2d , swt2d , scv2d , gwet2d , veg2d ,  &
                     & veg2d1 , sag2d , sice2d , dew2d , ircp2d ,       &
                     & text2d , col2d , ocld2d , heatrt , o3prof
          write (iutsav) pptnc , pptc , prca2d , prnca2d
          if ( iocnflx.eq.2 ) write (iutsav) sfsta%zpbl
          if ( ichem.eq.1 ) then
            write (iutsav) chia
            write (iutsav) chib
!           cumul removal terms (3d, 2d)
            write (iutsav) remlsc
            write (iutsav) remcvc
            write (iutsav) remdrd
          end if
          write (iutsav) ssw2da
          write (iutsav) sdeltk2d
          write (iutsav) sdelqk2d
          write (iutsav) sfracv2d
          write (iutsav) sfracb2d
          write (iutsav) sfracs2d
          write (iutsav) svegfrac2d
!         cumul ad, dif, emis terms ( scalar)
#ifndef BAND
          if ( ichem.eq.1 ) then
            if (debug_level > 2) call savechemdiag(iutsav)
          end if
#endif
#endif
          if ( lakemod.eq.1 ) then
            call lakesav_o(iutsav)
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
          write (iutsav) spsav%dstor
          write (iutsav) spsav%hstor
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

          write (6,*) 'SAV variables written at ', idate, xtime

          if (isavlast > 0) then
            write (fbname, '(a,i10)') 'TMPSAV.', isavlast
            ffout = trim(dirout)//pthsep//trim(domname)// &
                    '_'//trim(fbname)
            call unlink(ffout)
          end if
          if (ltmp) then
            isavlast = idate
          else
            isavlast = 0
          end if

        end subroutine write_savefile

      end module mod_savefile
