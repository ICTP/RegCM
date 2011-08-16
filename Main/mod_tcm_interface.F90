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

module mod_tcm_interface

  use mod_memutil
  use mod_constants
  use mod_runparams
  use mod_slice
  use mod_rad
  use mod_atm_interface
  use mod_gridfuncs
  use mod_service
  use mod_pbldim
  use m_realkinds

  private

  type , public :: tcm_state
    !
    ! TKE*ps (m^2/s^2 * cb)
    !
    real(dp) , pointer , dimension(:,:,:) :: tkeps
    !
    ! Coupled TKE Advective Tendency (m^2/s^3 * cb) 
    !
    real(dp) , pointer , dimension(:,:,:) :: advtke
    !
    ! Vertical momentum diffusivity (m^2/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: kzm
    !
    ! Vertical scalar diffusivity (m^2/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: kth
    !
    ! Boundary layer height (m)
    !
    real(dp) , pointer , dimension(:,:) :: zpbl
    !
    ! Surface layer TKE (m^2/s^2)
    !
    real(dp) , pointer , dimension(:,:) :: srftke
    !
  end type tcm_state

  !
  ! Custom types that hold the state of the host model
  !
  type , public :: host_atm_state
    !
    ! Zonal Wind (m/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: u
    !
    ! Meridional Wind (m/s)
    !
    real(dp) , pointer , dimension(:,:,:) :: v
    !
    ! Temperature (K)
    !
    real(dp) , pointer , dimension(:,:,:) :: t
    !
    ! Specific Humidity (kg/kg)
    !
    real(dp) , pointer , dimension(:,:,:) :: qv
    !
    ! Specific cloud water ratio (kg/kg)
    !
    real(dp) , pointer , dimension(:,:,:) :: qc
    !
    ! Turbulent Kinetic Energy (m^2/s^2)
    !
    real(dp) , pointer , dimension(:,:,:) :: tke 
    !
  end type host_atm_state

  type , public :: host_srf_state
    !
    ! Pressure (cb)
    !
    real(dp) , pointer , dimension(:,:) :: ps
    !
    ! Ground temperature (K)
    !
    real(dp) , pointer , dimension(:,:) :: tg
    !
    ! Surface Latent Enthalpy Flux (W/m^2)
    !
    real(dp) , pointer , dimension(:,:) :: qfx
    !
    ! Surface Sensible Enthalpy Flux (W/m^2)
    !
    real(dp) , pointer , dimension(:,:) :: hfx
    !
    ! Surface Drag (~momentum flux) (kg/m^2/s [?])
    !
    real(dp) , pointer , dimension(:,:) :: uvdrag
    !
  end type host_srf_state

  type , public :: host_rad_state
    !
    ! Radiative heating (cooling) rate
    real(dp) , pointer , dimension(:,:,:) :: heatrt
  end type host_rad_state

  type,public :: host_domain
    !
    ! Sigma at full levels
    !
    real(dp) , pointer , dimension(:) :: sigma
    !
    ! Sigma at half levels
    !
    real(dp) , pointer , dimension(:) :: a
    !
    ! Pressure at top of model (cb)
    !
    real(dp) :: ptop
    !
  end type host_domain

  integer , public :: itcmstart , itcmend
  integer , public :: jtcmstart , jtcmend

  !
  ! Pointers to point to the TCM's state variable
  !
  type(tcm_state) , public , pointer :: tcmstatea , tcmstateb

  real(dp) , public :: dttke ! TKE time step
  real(dp) , public :: tkemin

  !
  ! Specific instances of the model's state variables (at the b time step)
  !
  type(host_atm_state) , public :: atmstateb
  type(host_srf_state) , public :: srfstateb
  type(host_rad_state) , public :: radstateb
  !
  ! Parameters for the host model's domain
  !
  type(host_domain) , public :: hdomain

  public :: init_tcm_interface , end_tcm_interface
  public :: allocate_tcm_state , deallocate_tcm_state
  public :: get_data_from_tcm
  public :: hadvTKE , vadvTKE , set_tke_bc , check_conserve_qt

  contains

  subroutine init_tcm_interface
    implicit none

    itcmstart = 2
    itcmend = iy-1 
    if ( myid == 0 ) then
      jtcmstart = 2
    else
      jtcmstart = 1
    end if

    if ( myid == nproc-1 ) then
      jtcmend = jxp-1 
    else
      jtcmend = jxp 
    end if

    atmstateb%u => atms%ubx3d
    atmstateb%v => atms%vbx3d
    atmstateb%t => atms%tb3d
    atmstateb%qv => atms%qvb3d
    atmstateb%qc => atms%qcb3d
    atmstateb%tke => atm2%tke

    srfstateb%ps => sps2%ps
    srfstateb%tg => sts2%tg
    srfstateb%qfx => sfsta%qfx
    srfstateb%hfx => sfsta%hfx
    srfstateb%uvdrag => sfsta%uvdrag

    radstateb%heatrt => heatrt

    hdomain%sigma => sigma
    hdomain%a => a
    hdomain%ptop = r8pt 

    dttke = dt
    tkemin = 1.0D-8

  end subroutine init_tcm_interface

  subroutine allocate_tcm_state(tcmstate,lpar)
    implicit none
    type(tcm_state) , intent(out) :: tcmstate
    logical , intent(in) :: lpar
    !
    ! Allocate the tcm state variables
    !
    if ( lpar ) then
      call getmem3d(tcmstate%tkeps,1,iy,1,kzp1,1,jxp,'mod_tcm_interface:tkeps')
      call getmem3d(tcmstate%advtke,1,iy,1,kzp1,1,jxp, &
                    'mod_tcm_interface:advtke')
      call getmem3d(tcmstate%kzm,1,iy,1,kzp1,1,jxp,'mod_tcm_interface:kzm')
      call getmem3d(tcmstate%kth,1,iy,1,kzp1,1,jxp,'mod_tcm_interface:kth')
      call getmem2d(tcmstate%zpbl,1,iy,1,jxp,'mod_tcm_interface:zpbl')
      call getmem2d(tcmstate%srftke,1,iy,1,jxp,'mod_tcm_interface:srftke')
    else
      call getmem3d(tcmstate%tkeps,1,iy,1,kzp1,1,jx,'mod_tcm_interface:tkeps')
      call getmem3d(tcmstate%advtke,1,iy,1,kzp1,1,jx,'mod_tcm_interface:advtke')
      call getmem3d(tcmstate%kzm,1,iy,1,kzp1,1,jx,'mod_tcm_interface:kzm')
      call getmem3d(tcmstate%kth,1,iy,1,kzp1,1,jx,'mod_tcm_interface:kth')
      call getmem2d(tcmstate%zpbl,1,iy,1,jx,'mod_tcm_interface:zpbl')
      call getmem2d(tcmstate%srftke,1,iy,1,jx,'mod_tcm_interface:srftke')
    end if
  end subroutine allocate_tcm_state

  subroutine end_tcm_interface()
    implicit none
    nullify(atmstateb%u)
    nullify(atmstateb%v)
    nullify(atmstateb%t)
    nullify(atmstateb%qv)
    nullify(atmstateb%qc)
    nullify(atmstateb%tke)

    nullify(srfstateb%ps)
    nullify(srfstateb%tg)
    nullify(srfstateb%qfx)
    nullify(srfstateb%hfx)
    nullify(srfstateb%uvdrag)

    nullify(radstateb%heatrt)

    nullify(hdomain%sigma)
    nullify(hdomain%a)
  end subroutine end_tcm_interface

  subroutine get_data_from_tcm(tcmstate,tcmtend,bRegridWinds)
    implicit none 
    type(atmstate) , intent(inout) :: tcmtend
    type(tcm_state) , intent(inout) :: tcmstate
    logical , intent(in) :: bRegridWinds
    integer :: k

    ! Don't update the model variables if we are the diagnostic mode
    ! (Holtslag running, and UW updating TKE)
    if ( ibltyp /= 99 ) then
      !
      ! Put the t and qv tendencies in to difft and diffq for
      ! application of the sponge boundary conditions (see mod_tendency)
      !
      diffq(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                  diffq(itcmstart:itcmend,:,jtcmstart:jtcmend) +  &
                  tcmtend%qv(itcmstart:itcmend,:,jtcmstart:jtcmend)
      difft(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                  difft(itcmstart:itcmend,:,jtcmstart:jtcmend) +  &
                  tcmtend%t(itcmstart:itcmend,:,jtcmstart:jtcmend)

      ! Put the cloud water tendency in aten
      aten%qc(itcmstart:itcmend,:,jtcmstart:jtcmend) =  &
                 aten%qc(itcmstart:itcmend,:,jtcmstart:jtcmend) +   &
                 tcmtend%qc(itcmstart:itcmend,:,jtcmstart:jtcmend)

      if ( bRegridWinds ) then
        !
        ! If the TCM calculations were all done on the cross grid, then
        ! the u and v tendencies need to be regridded to the dot grid
        !
        call uvcross2dot(tcmtend,aten)
      else
        ! Otherwise, simply set the appropriate variables in aten
        aten%u(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                 aten%u(itcmstart:itcmend,:,jtcmstart:jtcmend) +    &
                 tcmtend%u(itcmstart:itcmend,:,jtcmstart:jtcmend)
        aten%v(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                 aten%v(itcmstart:itcmend,:,jtcmstart:jtcmend) +    &
                 tcmtend%v(itcmstart:itcmend,:,jtcmstart:jtcmend)
      end if

      sfsta%zpbl(itcmstart:itcmend,jtcmstart:jtcmend) = &
                tcmstate%zpbl(itcmstart:itcmend,jtcmstart:jtcmend)
    end if

!   !
!   ! Interpolate kzm and kth from the interfaces to the midpoints
!   ! if the diffusivities are to be used in the holtslag model
!   !
!   if ( ibltyp == 99 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,k,:) = sqrt(tcmstate%kzm(:,k,:)*tcmstate%kzm(:,k+1,:))
!       tcmstate%kth(:,k,:) = sqrt(tcmstate%kth(:,k,:)*tcmstate%kth(:,k+1,:))
!     end do
!   end if
!   !
!   ! Shift kth and kzm for output if the UW model is running
!   !
!   if ( ibltyp == 2 ) then
!     do k = 1 , kz
!       tcmstate%kzm(:,k,:) = tcmstate%kzm(:,k+1,:)
!       tcmstate%kth(:,k,:) = tcmstate%kth(:,k+1,:)
!     end do
!   end if

    aten%tke(itcmstart:itcmend,:,jtcmstart:jtcmend) = &
                tcmtend%tke(itcmstart:itcmend,:,jtcmstart:jtcmend)
    !
    ! Set the surface TKE (diagnosed)
    !
    atm1%tke(itcmstart:itcmend,kzp1,jtcmstart:jtcmend) =  &
               tcmstate%srftke(itcmstart:itcmend,jtcmstart:jtcmend)
    atm2%tke(itcmstart:itcmend,kzp1,jtcmstart:jtcmend) =  &
               tcmstate%srftke(itcmstart:itcmend,jtcmstart:jtcmend)

!   tcmtend%qc = 0.0d0
!   tcmtend%qv = 0.0d0
!   tcmtend%t = 0.0d0
!   tcmtend%u = 0.0d0
!   tcmtend%v = 0.0d0
!   tcmtend%tke = 0.0d0


  end subroutine get_data_from_tcm

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the horizontal flux-divergence terms   c
!     for TKE.  Second-order difference is used.                      c
!                                                                     c
!     dxx    : is the horizontal distance.                            c
!     j      : is the j'th slice of f anf ften.                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine hadvTKE(tcmstate,dxx,j)
    implicit none
    integer , intent(in) :: j
    real(dp) , intent(in) :: dxx
    type(tcm_state) , pointer :: tcmstate
    character (len=50) :: subroutine_name='hadvTKE'
    integer :: idindx = 0
    integer :: i , k , jm1 , jp1
!
    call time_begin(subroutine_name,idindx)

    jm1 = j - 1
    jp1 = j + 1

    do k = 2 , kz
      do i = 2 , iym1

        !
        ! Interpoalte the winds to the full sigma levels
        ! while the advection term is calculated
        !
        tcmstate%advtke(i,k,j)= tcmstate%advtke(i,k,j)  &
           -(((atm1%u(i+1,k-1,jp1)+atm1%u(i,k-1,jp1))*twt(k,2)  &
             +(atm1%u(i+1,k,jp1)  +atm1%u(i,k,jp1))*twt(k,1)) &
             *( atm1%tke(i,k,j)+atm1%tke(i,k,jp1))  &
            -((atm1%u(i+1,k-1,j)+atm1%u(i,k-1,j))*twt(k,2)  &
             +(atm1%u(i+1,k,j)  +atm1%u(i,k,j))*twt(k,1)) &
             *( atm1%tke(i,k,j)+atm1%tke(i,k,jm1))  &
            +((atm1%v(i+1,k-1,j)+atm1%v(i+1,k-1,jp1))*twt(k,2)  &
             +(atm1%v(i+1,k,j)  +atm1%v(i+1,k,jp1))*twt(k,1)) &
             *( atm1%tke(i,k,j)+atm1%tke(i+1,k,j))  &
            -((atm1%v(i,k-1,j)+atm1%v(i,k-1,jp1))*twt(k,2)  &
             +(atm1%v(i,k,j)  +atm1%v(i,k,jp1))*twt(k,1)) &
             *( atm1%tke(i,k,j)+atm1%tke(i-1,k,j))) &
             /(dxx*mddom%msfx(i,j)*mddom%msfx(i,j))

!TAO Debug:
!       tcmstate%advtke(i,k,j)= tcmstate%advtke(i,k,j)  + 1d-8*atm1%tke(i,k,jp1)

      end do
    end do

    call time_end(subroutine_name,idindx)

  end subroutine hadvTKE

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     j      : jth slice of variable fa.                              c
!                                                                     c
!     ind = 1 : Bretherton's vertical advection method                c
!           2 : Alternate vertical advection method (unknown origin)  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine vadvTKE(tcmstate,j,ind)
    implicit none
    integer , intent(in) :: j , ind
    type(tcm_state) , pointer :: tcmstate
    real(dp) , dimension(iy,kz) :: dotqdot , ftmp
    character (len=50) :: subroutine_name='vadvTKE'
    integer :: idindx = 0
!
    integer :: i , k

    call time_begin(subroutine_name,idindx)
        
    !
    ! Use Bretherton's method for TKE advection
    !
    if ( ind == 1 ) then
      do k = 1 , kz
        do i = 2 , iym1
          dotqdot(i,k) = (qdot(i,k,j) + qdot(i,k+1,j))
          ftmp(i,k) = 0.5D0*(tcmstate%tkeps(i,k,j) + tcmstate%tkeps(i,k+1,j))
        end do
      end do

      do k = 2 , kz
        do i = 2 , iym1
          tcmstate%advtke(i,k,j) = tcmstate%advtke(i,k,j) - &
                                    (dotqdot(i,k)*ftmp(i,k) -  &
                                    dotqdot(i,k-1)*ftmp(i,k-1))/ &
                                    (dsigma(k)+dsigma(k-1))
        end do
      end do
    !
    ! Use an alternative method (this came from where?)
    !
    else
      do i = 2 , iym1
        tcmstate%advtke(i,1,j)=tcmstate%advtke(i,1,j)-  &
                              qdot(i,2,j)*tcmstate%tkeps(i,2,j)/dsigma(1)
      end do
      do k = 2 , kzm1
        do i = 2 , iym1
          tcmstate%advtke(i,k,j)=tcmstate%advtke(i,k,j) &
                              -(qdot(i,k+1,j)*tcmstate%tkeps(i,k+1,j)   &
                                -qdot(i,k,j)*tcmstate%tkeps(i,k,j))/dsigma(k)
        end do
      end do
      do i = 2 , iym1
        tcmstate%advtke(i,kz,j)=tcmstate%advtke(i,kz,j)+  &
                              qdot(i,kz,j)*tcmstate%tkeps(i,kz,j)/dsigma(kz)
      end do
    end if

    call time_end(subroutine_name,idindx)

  end subroutine vadvTKE

  subroutine set_tke_bc()
    implicit none
    !
    ! Set TKE boundary conditions
    !
    if ( myid == 0 ) then
      atm1%tke(:,:,0:1) = tkemin ! East boundary
      atm2%tke(:,:,-1:1) = tkemin ! East boundary
    end if
    if ( myid == nproc-1 ) then
      atm1%tke(:,:,jxp:jxp+1) = tkemin ! West boundary
      atm2%tke(:,:,jxp:jxp+2) = tkemin ! West boundary
    end if

    atm1%tke(1,:,:) = tkemin  !North boundary
    atm2%tke(1,:,:) = tkemin  !North boundary
    atm1%tke(iy,:,:) = tkemin  !South boundary
    atm2%tke(iy,:,:) = tkemin  !South boundary
    !
    ! End set TKE boundary conditions
    !
  end subroutine set_tke_bc

  subroutine check_conserve_qt(rcmqvten,rcmqcten,tcmtend,dom,tcmstate,kmax)
    implicit none
    type(atmstate) , intent(in) :: tcmtend
    type(tcm_state) :: tcmstate
    type(host_domain) , intent(in) :: dom
    real(dp) , dimension(:,:,:) , intent(in) :: rcmqvten , rcmqcten
    integer , intent(in) :: kmax
    real(dp) , dimension(kmax) :: rho1d , rhobydpdz1d
    real(dp) :: qwtcm , qwrcm , qwanom , dtops , xps , ps2 , dza
    integer :: i , j , k

    do i = 2 , iym1
      do j = jbegin , jendx

        do k = 1 , kzm1
          xps = (dom%a(k)*sps2%ps(i,j)+dom%ptop)
          ps2 = (dom%a(k+1)*sps2%ps(i,j)+dom%ptop)
          dza = za(i,k,j) - za(i,k+1,j)
          rhobydpdz1d(k) = d_1000*(ps2-xps)/(egrav*dza)

        end do
        rhobydpdz1d(kz) = rhox2d(i,j)
        dtops = dt/sps2%ps(i,j)

        rho1d = d_1000*(dom%a*sps2%ps(i,j) + dom%ptop) / &
                    ( rgas * atms%tb3d(i,:,j) *  &
                    (d_one + ep1* atms%qvb3d(i,:,j) - atms%qcb3d(i,:,j))  )


        qwtcm = sum((tcmtend%qv(i,:,j) + tcmtend%qc(i,:,j))  &
                      *rho1d*dzq(i,:,j))*dtops
        qwrcm = sum((rcmqvten(i,:,j) + rcmqcten(i,:,j))   &
                      *rho1d*dzq(i,:,j))*dtops
!                     *rhobydpdz1d*dzq(i,:,j))*dtops

!       qwanom = qwtcm - qwrcm
!       qwanom = qwrcm-sfsta%qfx(i,j)*dt
        qwanom = qwtcm - dt*sfsta%qfx(i,j)
        tcmstate%kzm(i,1,j) = qwanom

      end do
    end do
  end subroutine check_conserve_qt

end module mod_tcm_interface
