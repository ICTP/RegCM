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

      module mod_runparams

      use mod_dynparam
      use mod_message
      use mod_service 

      implicit none
 
      integer :: ifrabe , klake , nbatst
!
      real(8) :: dt , dt0 , dt2 , dtbat , dtlake , dtmin
      real(8) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
      real(8) :: abatm , abemh
      real(8) :: c200 , c203 
      real(8) :: fnudge , gnudge
      real(8) :: xkhmax , xkhz

      integer :: ibltyp , iboudy , ichem , icup , idirect ,      &
               & iemiss , igcc , iocnflx , ipgf , ipptls , kbats ,      &
               & kchem , lakemod , nradisp , ntrad , ntsave , nttape ,  &
               & idcsst , iseaice , idesseas , klak

      logical :: ifrest , rfstrt
 
      real(8) :: bdytim , prttim , radfrq , savtim , taptim , tbdybe
      integer :: ndbgfrq , nsavfrq , ntapfrq

      integer :: ispgd , ispgx , k700 , kchi , kclo , kcmd , kt , ncld
!
      real(8) :: r8pt
      real(8) :: akht1 , akht2
      real(8) :: high_nudge , medium_nudge , low_nudge

      real(8) , allocatable , dimension(:) :: dtau
      real(8) , allocatable , dimension(:) :: a , anudg , dsigma , qcon
      real(8) , allocatable , dimension(:) :: sigma
      real(8) , allocatable , dimension(:,:) :: twt
      real(8) , allocatable , dimension(:) :: wgtd
      real(8) , allocatable , dimension(:) :: wgtx

      character(len=3) :: scenario

      integer, private  :: ierr 
      real(8) , private :: total_allocation_size
      data total_allocation_size /0.0/

      contains

      subroutine allocate_mod_runparams
        implicit none
        character (len=10) :: myname='run-params'
        allocate(a(kz),stat=ierr)
        call check_alloc(ierr,myname,'a',size(a)*kind(a))
        allocate(anudg(kz),stat=ierr)
        call check_alloc(ierr,myname,'anudg',size(anudg)*kind(anudg))
        allocate(dsigma(kz),stat=ierr)
        call check_alloc(ierr,myname,'dsigma',size(dsigma)*kind(dsigma))
        allocate(qcon(kz),stat=ierr)
        call check_alloc(ierr,myname,'qcon',size(qcon)*kind(dsigma))
        allocate(sigma(kzp1),stat=ierr)
        call check_alloc(ierr,myname,'sigma',size(sigma)*kind(sigma))
        allocate(twt(kz,2),stat=ierr)
        call check_alloc(ierr,myname,'twt',size(twt)*kind(twt))
        allocate(wgtd(nspgd),stat=ierr)
        call check_alloc(ierr,myname,'wgtd',size(wgtd)*kind(wgtd))
        allocate(wgtx(nspgx),stat=ierr)
        call check_alloc(ierr,myname,'wgtx',size(wgtx)*kind(wgtx))
        allocate(dtau(nsplit),stat=ierr)
        call check_alloc(ierr,myname,'dtau',size(dtau)*kind(dtau))

        a = 0.0D0
        anudg = 0.0D0
        dsigma = 0.0D0
        qcon = 0.0D0
        sigma = 0.0D0
        twt = 0.0D0
        wgtd = 0.0D0
        wgtx = 0.0D0
        dtau = 0.0D0
        
      call report_alloc('allocate_mod_run_params') 

      end subroutine allocate_mod_runparams
!
      subroutine check_alloc(ierr,where,what,isize)
        implicit none
        integer , intent(in) :: ierr , isize
        character(len=*) :: what , where
#ifdef DEBUG
        character(len=50) :: buffer
#endif

        real(8) :: storage 
        if (ierr /= 0) then
          call fatal(__FILE__,__LINE__,what//' CANNOT BE allocated')
        end if
        storage=isize*kind(wgtd(1))
#ifdef DEBUG 
        if (debug_level.gt.3) then 
           write (buffer,*)  what, &
             '  allocated succesfully: global size is'
           CALL write_info(where,buffer,isize)
        end if 
#endif           
        total_allocation_size = total_allocation_size + isize
      end subroutine check_alloc
!
      subroutine report_alloc(where)
        implicit none
        character(len=*) ::  where
        write(aline,*) where, &
          ': total allocation (in Kbyte)=', total_allocation_size*8/1024
        call say
      end subroutine report_alloc
!
      end module mod_runparams
