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
  use mod_mppparam
  use mod_atm_interface
  use mod_lm_interface
  use mod_che_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_bdycod
  use mod_mppio
#ifdef CLM
  use mod_clm
  use restFileMod, only : restFile_write, restFile_write_binary
  use restFileMod, only : restFile_filename
  use clm_varctl , only : filer_rest
  use clm_time_manager, only : get_step_size
#endif
  private

  public :: read_savefile
  public :: write_savefile

  integer(ik4) :: isavlast
  integer(ik4) :: iutrst
  integer(ik4) :: ncsave

  integer(ik4) , dimension(2) :: idimcross2d
  integer(ik4) , dimension(2) :: idimdot2d
  integer(ik4) , dimension(3) :: idimcross3d
  integer(ik4) , dimension(3) :: idimdot3d
  integer(ik4) , dimension(1) :: idimsigmafull
  integer(ik4) , dimension(1) :: idimsigmahalf
  integer(ik4) , dimension(3) :: idimsubgridcross2d

#ifdef CLM
  character(len=256) :: thisclmrest
  character(len=256) :: lastclmrest
#endif

  data isavlast /-1/
  data iutrst   /-1/

  contains

  subroutine read_savefile(idate)
    implicit none
    type (rcm_time_and_date) , intent(in) :: idate
    character(256) :: ffin
    character(16) :: fbname
    integer(ik8) :: idt1 , idt2
    real(rk8) :: odtsec
    logical :: existing

    if ( myid == iocpu ) then
      iutrst = 14
      write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      ffin = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
      inquire (file=ffin,exist=existing)
      if ( .not.existing ) then
        write (stderr,*) 'The following SAV File does not exist: ', &
                         trim(ffin), ' please check location'
        call fatal(__FILE__,__LINE__, 'SAV FILE NOT FOUND')
      else
        open (iutrst,file=ffin,form='unformatted',status='old')
      end if

      read (iutrst) ktau, odtsec, idatex
      idt1 = idnint(odtsec)
      idt2 = idnint(dtsec)
      if ( idt1 /= idt2 ) then
        write (stdout,*) 'Recalculating ktau for the new dt'
        write (stdout,*) 'Restart file ktau is       = ', ktau
        ktau = (ktau * idt1) / idt2
        write (stdout,*) 'Actual ktau with new dt is = ', ktau
        write (stdout,*) 'Done Recalculating ktau for the new dt'
      end if
      read (iutrst) atm1_io%u
      read (iutrst) atm1_io%v
      read (iutrst) atm1_io%t
      read (iutrst) atm1_io%qx
      read (iutrst) atm2_io%u
      read (iutrst) atm2_io%v
      read (iutrst) atm2_io%t
      read (iutrst) atm2_io%qx
      read (iutrst) sfs_io%psa , sfs_io%psb
      read (iutrst) sfs_io%tga , sfs_io%tgb
      read (iutrst) sfs_io%hfx , sfs_io%qfx
      read (iutrst) sfs_io%rainc , sfs_io%rainnc
      read (iutrst) sfs_io%tgbb , sfs_io%uvdrag
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        read (iutrst) atm1_io%tke
        read (iutrst) atm2_io%tke
        read (iutrst) kpbl_io
      end if
      if ( icup == 1 ) then
        read (iutrst) rsheat_io , rswat_io
      end if
      if ( icup == 3 ) then
        read (iutrst) tbase_io , cldefi_io
      end if
      if ( icup == 4 .or. icup == 99 .or. icup == 98 .or. icup == 97 ) then
        read (iutrst) cbmf2d_io
      end if
      if ( idcsst == 1 ) then
        read (iutrst) dtskin_io
        read (iutrst) deltas_io
        read (iutrst) tdeltas_io
      end if
      if ( irrtm == 0 ) then
        read (iutrst) gasabsnxt_io , gasabstot_io , gasemstot_io
      end if
      if ( ipptls > 0 ) then
        read (iutrst) fcc_io
        if ( ipptls == 2 ) then
          read (iutrst) sfs_io%snownc
        end if
      end if
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
      read (iutrst) solis_io
      read (iutrst) solvd_io
      read (iutrst) solvs_io
      read (iutrst) sabveg_io
      read (iutrst) tlef_io
      read (iutrst) ssw_io
      read (iutrst) rsw_io
      read (iutrst) tgrd_io
      read (iutrst) tgbrd_io
      read (iutrst) sncv_io
      read (iutrst) gwet_io
      read (iutrst) snag_io
      read (iutrst) sfice_io
      read (iutrst) ldew_io
      read (iutrst) ldmsk_io
      read (iutrst) heatrt_io
      read (iutrst) o3prof_io
      read (iutrst) flw_io
      read (iutrst) flwd_io
      read (iutrst) fsw_io
      read (iutrst) tsw_io
      read (iutrst) sinc_io
      read (iutrst) taf_io
      read (iutrst) ldmsk1_io
      read (iutrst) emiss_io
      if ( iocnflx == 2 ) then
        read (iutrst) zpbl_io
      end if
      if ( ichem == 1 ) then
        read (iutrst) chia_io
        read (iutrst) chib_io
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          read (iutrst) chemall_io
          read (iutrst) taucldsp_io
        end if
        ! cumul removal terms (3d, 2d)
        read (iutrst) remlsc_io
        read (iutrst) remcvc_io
        read (iutrst) remdrd_io
        ! cumul ad, dif, emis terms ( scalar)
        read (iutrst) ssw2da_io
        read (iutrst) sdeltk2d_io
        read (iutrst) sdelqk2d_io
        read (iutrst) sfracv2d_io
        read (iutrst) sfracb2d_io
        read (iutrst) sfracs2d_io
        read (iutrst) svegfrac2d_io
      end if
#ifndef CLM
      if ( lakemod == 1 ) then
        call lakesav_i(iutrst)
      end if
#endif
      read (iutrst) dstor_io
      read (iutrst) hstor_io
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        read (iutrst) qflux_restore_sst_io
        read (iutrst) stepcount
      end if
      close(iutrst)
    end if
  end subroutine read_savefile

  subroutine write_savefile(idate,ltmp)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    logical , intent(in) :: ltmp
    integer(ik4) , parameter :: iutsav = 52
    character(256) :: ffout
    character(32) :: fbname
    logical :: existing
#ifdef CLM
    integer(ik4) :: ioff
#endif
    if ( myid == iocpu ) then
      if (ltmp) then
        write (fbname, '(a,i10)') 'TMPSAV.', toint10(idate)
      else
        write (fbname, '(a,i10)') 'SAV.', toint10(idate)
      end if
      ffout = trim(dirout)//pthsep//trim(domname)//'_'//trim(fbname)
      open (iutsav,file=ffout,form='unformatted',status='replace')

      inquire (file=ffout,exist=existing)
      if ( .not.existing ) then
        write (stderr,*) 'The SAV File cannot be created: ' , &
                         trim(ffout), ' please check directory'
        call fatal(__FILE__,__LINE__, 'SAV FILE WRITE ERROR')
      end if

      write (iutsav) ktau , dtsec , idatex
      write (iutsav) atm1_io%u
      write (iutsav) atm1_io%v
      write (iutsav) atm1_io%t
      write (iutsav) atm1_io%qx
      write (iutsav) atm2_io%u
      write (iutsav) atm2_io%v
      write (iutsav) atm2_io%t
      write (iutsav) atm2_io%qx
      write (iutsav) sfs_io%psa , sfs_io%psb
      write (iutsav) sfs_io%tga , sfs_io%tgb
      write (iutsav) sfs_io%hfx , sfs_io%qfx
      write (iutsav) sfs_io%rainc , sfs_io%rainnc
      write (iutsav) sfs_io%tgbb , sfs_io%uvdrag
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        write (iutsav) atm1_io%tke
        write (iutsav) atm2_io%tke
        write (iutsav) kpbl_io
      end if
      if ( icup == 1 ) then
        write (iutsav) rsheat_io , rswat_io
      end if
      if ( icup == 3 ) then
        write (iutsav) tbase_io , cldefi_io
      end if
      if ( icup == 4 .or. icup == 99 .or. icup == 98 .or. icup == 97 ) then
        write (iutsav) cbmf2d_io
      end if
      if ( idcsst == 1 ) then
        write (iutsav) dtskin_io
        write (iutsav) deltas_io
        write (iutsav) tdeltas_io
      end if
      if ( irrtm == 0 ) then
        write (iutsav) gasabsnxt_io , gasabstot_io , gasemstot_io
      end if
      if ( ipptls > 0 ) then
        write (iutsav) fcc_io
        if ( ipptls == 2 ) then
          write (iutsav) sfs_io%snownc
        end if
      end if
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
      write (iutsav) solis_io
      write (iutsav) solvd_io
      write (iutsav) solvs_io
      write (iutsav) sabveg_io
      write (iutsav) tlef_io
      write (iutsav) ssw_io
      write (iutsav) rsw_io
      write (iutsav) tgrd_io
      write (iutsav) tgbrd_io
      write (iutsav) sncv_io
      write (iutsav) gwet_io
      write (iutsav) snag_io
      write (iutsav) sfice_io
      write (iutsav) ldew_io
      write (iutsav) ldmsk_io
      write (iutsav) heatrt_io
      write (iutsav) o3prof_io
      write (iutsav) flw_io
      write (iutsav) flwd_io
      write (iutsav) fsw_io
      write (iutsav) tsw_io
      write (iutsav) sinc_io
      write (iutsav) taf_io
      write (iutsav) ldmsk1_io
      write (iutsav) emiss_io
      if ( iocnflx == 2 ) then
        write (iutsav) zpbl_io
      end if
      if ( ichem == 1 ) then
        write (iutsav) chia_io
        write (iutsav) chib_io
        if ( igaschem == 1 .and. ichsolver > 0 ) then
          write (iutsav) chemall_io
          write (iutsav) taucldsp_io
        end if
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
      end if
#ifndef CLM
      if ( lakemod == 1 ) then
        call lakesav_o(iutsav)
      end if
#endif
      write (iutsav) dstor_io
      write (iutsav) hstor_io
      if ( islab_ocean == 1 .and. do_restore_sst ) then
        write (iutsav) qflux_restore_sst_io
        write (iutsav) stepcount
      end if
      close(iutsav)
    end if

#ifdef CLM
    ioff = dtsec*(ntsrf-1)
    filer_rest = restFile_filename(type='netcdf',offset=ioff)
    call restFile_write(filer_rest)
    filer_rest = restFile_filename(type='binary',offset=ioff)
    call restFile_write_binary(filer_rest)
    thisclmrest = filer_rest(1:256)
#endif
    if ( myid == iocpu ) then
      write(stdout,*) 'SAV variables written at ', tochar(idate)
      if ( isavlast > 0 ) then
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
