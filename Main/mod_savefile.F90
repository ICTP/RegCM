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

  use mod_runparams
  use mod_mpmessage
  use mod_atm_interface
  use mod_lm_interface
  use mod_che_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_bdycod
#ifndef BAND
  use mod_diagnosis
#endif
  use mod_mppio
#ifdef CLM
  use mod_clm
  use restFileMod, only : restFile_write, restFile_write_binary
  use restFileMod, only : restFile_filename
  use clm_varctl , only : filer_rest
  use clm_time_manager, only : get_step_size
#endif
  private

  public :: read_savefile_part1 , read_savefile_part2
  public :: write_savefile

  integer :: isavlast
  integer :: iutrst
  logical :: lrp1
#ifdef CLM
  character(len=256) :: thisclmrest
  character(len=256) :: lastclmrest
#endif

  data isavlast /-1/
  data iutrst   /-1/
  data lrp1 /.false./

  contains

  subroutine read_savefile_part1(idate)
    implicit none
    type (rcm_time_and_date) , intent(in) :: idate
    character(256) :: ffin
    character(16) :: fbname
    logical :: existing

    if ( myid == 0 ) then
      iutrst = 14
      write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      ffin = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
      inquire (file=ffin,exist=existing)
      if ( .not.existing ) then
        write (aline,*) 'The following SAV File does not exist: ', &
                         trim(ffin), ' please check location'
        call say
        call fatal(__FILE__,__LINE__, 'SAV FILE NOT FOUND')
      else
        open (iutrst,file=ffin,form='unformatted',status='old')
      end if

      read (iutrst) ktau, idatex, nbdytime
      if ( ehso4 ) then
        read (iutrst) ub0_io, vb0_io, qb0_io, tb0_io, ps0_io, ts0_io, so0_io
      else
        read (iutrst) ub0_io, vb0_io, qb0_io, tb0_io, ps0_io, ts0_io
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
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        read (iutrst) atm1_io%tke
        read (iutrst) atm2_io%tke
      end if
      read (iutrst) kpbl_io
      if ( icup == 1 ) then
        read (iutrst) rsheat_io , rswat_io
      end if
      if ( icup == 3 ) then
        read (iutrst) tbase_io , cldefi_io
      end if
      if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
        read (iutrst) cbmf2d_io
      end if
      read (iutrst) hfx_io , qfx_io , uvdrag_io
#ifndef BAND
      call restdiag(iutrst)
#endif
      read (iutrst) absnxt_io , abstot_io , emstot_io
      if ( ipptls == 1 ) read (iutrst) fcc_io
#ifdef CLM
      read (iutrst) sols2d_io
      read (iutrst) soll2d_io
      read (iutrst) solsd2d_io
      read (iutrst) solld2d_io
      read (iutrst) aldirs2d_io
      read (iutrst) aldirl2d_io
      read (iutrst) aldifs2d_io
      read (iutrst) aldifl2d_io
      read (iutrst) lndcat2d_io
#endif
      read (iutrst) sol2d_io
      read (iutrst) solvd2d_io
      read (iutrst) solvs2d_io
      read (iutrst) sabv2d_io
      read (iutrst) tlef2d_io
      read (iutrst) ssw2d_io
      read (iutrst) srw2d_io
      read (iutrst) tg2d_io
      read (iutrst) tgb2d_io
      read (iutrst) scv2d_io
      read (iutrst) gwet2d_io
      read (iutrst) sag2d_io
      read (iutrst) sice2d_io
      read (iutrst) dew2d_io
      read (iutrst) ircp2d_io
      read (iutrst) col2d_io
      read (iutrst) veg2d_io
      read (iutrst) ldmsk_io
      read (iutrst) veg2d1_io
      read (iutrst) heatrt_io
      read (iutrst) o3prof_io
      read (iutrst) tgbb_io
      read (iutrst) flw2d_io
      read (iutrst) flwd2d_io
      read (iutrst) fsw2d_io
      read (iutrst) swt2d_io
      read (iutrst) sinc2d_io
      read (iutrst) taf2d_io
      read (iutrst) ocld2d_io
      read (iutrst) emiss2d_io
      read (iutrst) pptnc_io, pptc_io, prca2d_io, prnca2d_io
      if ( iocnflx == 2 ) read (iutrst) zpbl_io
      if ( ichem == 1 ) then
        read (iutrst) chia_io
        read (iutrst) chib_io
!             cumul removal terms (3d, 2d)
        read (iutrst) remlsc_io
        read (iutrst) remcvc_io
        read (iutrst) remdrd_io
!             cumul ad, dif, emis terms ( scalar)
        read (iutrst) ssw2da_io
        read (iutrst) sdeltk2d_io
        read (iutrst) sdelqk2d_io
        read (iutrst) sfracv2d_io
        read (iutrst) sfracb2d_io
        read (iutrst) sfracs2d_io
        read (iutrst) svegfrac2d_io
#ifndef BAND
        call restchemdiag(iutrst)
#endif
      end if
!------lake model
      if ( lakemod == 1 ) then
        call lakesav_i(iutrst)
      end if
      lrp1 = .true.
    end if
  end subroutine read_savefile_part1

  subroutine read_savefile_part2
    implicit none

    if ( myid == 0 ) then
      if (.not. lrp1) then
        write (6,*) 'Reading part2 before part1'
        call fatal(__FILE__,__LINE__, 'SAV FILE ERROR')
      end if

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
      close(iutrst)
    end if
  end subroutine read_savefile_part2

  subroutine write_savefile(idate,ltmp)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: ltmp
    integer , parameter :: iutsav = 52
    character(256) :: ffout
    character(32) :: fbname
    logical :: existing
#ifdef CLM
    real(8) :: cdtime
#endif
    if ( myid == 0 ) then
      if (ltmp) then
        write (fbname, '(a,i10)') 'TMPSAV.', toint10(idate)
      else
        write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      end if
      ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
      open (iutsav,file=ffout,form='unformatted',status='replace')

      inquire (file=ffout,exist=existing)
      if ( .not.existing ) then
        write (aline,*) 'The SAV File cannot be created: ' , &
                         trim(ffout), ' please check directory'
        call say
        call fatal(__FILE__,__LINE__, 'SAV FILE WRITE ERROR')
      end if

      write (iutsav) ktau , idatex , nbdytime
      if ( ehso4 ) then
        write (iutsav) ub0_io , vb0_io , qb0_io , tb0_io , &
                       ps0_io , ts0_io , so0_io
      else
        write (iutsav) ub0_io , vb0_io , qb0_io , tb0_io , &
                       ps0_io , ts0_io
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
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        write (iutsav) atm1_io%tke
        write (iutsav) atm2_io%tke
      end if
      write (iutsav) kpbl_io
      if ( icup == 1 ) then
        write (iutsav) rsheat_io , rswat_io
      end if
      if ( icup == 3 ) then
        write (iutsav) tbase_io , cldefi_io
      end if
      if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
        write (iutsav) cbmf2d_io
      end if
      write (iutsav) hfx_io , qfx_io , uvdrag_io
#ifndef BAND
      call savediag(iutsav)
#endif
      write (iutsav) absnxt_io , abstot_io , emstot_io
      if ( ipptls == 1 ) write (iutsav) fcc_io
#ifdef CLM
      write (iutsav) sols2d_io
      write (iutsav) soll2d_io
      write (iutsav) solsd2d_io
      write (iutsav) solld2d_io
      write (iutsav) aldirs2d_io
      write (iutsav) aldirl2d_io
      write (iutsav) aldifs2d_io
      write (iutsav) aldifl2d_io
      write (iutsav) lndcat2d_io
#endif
      write (iutsav) sol2d_io
      write (iutsav) solvd2d_io
      write (iutsav) solvs2d_io
      write (iutsav) sabv2d_io
      write (iutsav) tlef2d_io
      write (iutsav) ssw2d_io
      write (iutsav) srw2d_io
      write (iutsav) tg2d_io
      write (iutsav) tgb2d_io
      write (iutsav) scv2d_io
      write (iutsav) gwet2d_io
      write (iutsav) sag2d_io
      write (iutsav) sice2d_io
      write (iutsav) dew2d_io
      write (iutsav) ircp2d_io
      write (iutsav) col2d_io
      write (iutsav) veg2d_io
      write (iutsav) ldmsk_io
      write (iutsav) veg2d1_io
      write (iutsav) heatrt_io
      write (iutsav) o3prof_io
      write (iutsav) tgbb_io
      write (iutsav) flw2d_io
      write (iutsav) flwd2d_io
      write (iutsav) fsw2d_io
      write (iutsav) swt2d_io
      write (iutsav) sinc2d_io
      write (iutsav) taf2d_io
      write (iutsav) ocld2d_io
      write (iutsav) emiss2d_io
      write (iutsav) pptnc_io , pptc_io , prca2d_io , prnca2d_io
      if ( iocnflx == 2 ) write (iutsav) zpbl_io
      if ( ichem == 1 ) then
        write (iutsav) chia_io
        write (iutsav) chib_io
!             cumul removal terms (3d, 2d)
        write (iutsav) remlsc_io
        write (iutsav) remcvc_io
        write (iutsav) remdrd_io
!             cumul ad, dif, emis terms ( scalar)
        write (iutsav) ssw2da_io
        write (iutsav) sdeltk2d_io
        write (iutsav) sdelqk2d_io
        write (iutsav) sfracv2d_io
        write (iutsav) sfracb2d_io
        write (iutsav) sfracs2d_io
        write (iutsav) svegfrac2d_io
#ifndef BAND
        call savechemdiag(iutsav)
#endif
      end if
      if ( lakemod == 1 ) then
        call lakesav_o(iutsav)
      end if
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
      close(iutsav)
    end if

#ifdef CLM
    cdtime = dble(get_step_size())
    filer_rest = restFile_filename(type='netcdf',offset=-idint(cdtime))
    call restFile_write(filer_rest)
    filer_rest = restFile_filename(type='binary',offset=-idint(cdtime))
    call restFile_write_binary(filer_rest)
    thisclmrest = filer_rest(1:256)
#endif
    if ( myid == 0 ) then
      print *, 'SAV variables written at ', tochar(idate)

      if (isavlast > 0) then
        write (fbname, '(a,i10)') 'TMPSAV.', isavlast
        ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
        call unlink(ffout)
#ifdef CLM
        call unlink(trim(lastclmrest))
        call unlink((trim(lastclmrest)//'.nc'))
#endif            
      end if
      if (ltmp) then
        isavlast = toint10(idate)
      else
        isavlast = 0
      end if
#ifdef CLM
      lastclmrest = thisclmrest
#endif
    end if
  end subroutine write_savefile

end module mod_savefile
