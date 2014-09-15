program cumclouds
  use netcdf
  implicit none
  integer , parameter :: nk = 23     ! Num vert layers
  integer , parameter :: nds = 8     ! Num testing resolutions

  integer , parameter :: maxtopcloud = 6 ! Max cloud top found in model

  real(8) , parameter :: dsinc = 25.0D0 ! Increment of ds
  real(8) , parameter :: clfrcvmax = 1.0D0 ! Maximum convective cloud fraction

  real(8) , parameter :: maxcloud_dp = 60.0D0 ! In cb, cloud depth to have
                                              ! Maximum cloud for this grid
                                              ! resolution

  real(8) , parameter :: deltap = 5.0 ! Testing deltap

  integer :: k , ktop , kbot , ids , kclth , ikh , irec
  integer :: istat , ncid
  integer , dimension(2) :: idims , istart , icount
  integer , dimension(5) :: ivars

  real(8) :: dp , ds , dxtemc , clfrcv , scalep , totcf , akclth
  real(8) :: scalef , totcf_old
  real(8) , dimension(1) :: helpme
  
  real(8) :: afracl = 0.250D0
  real(8) :: afracs = clfrcvmax

  ! A coarse grid is 200 km size , a high resolution is 10 km
  real(8) :: dlargc = 200.0D0
  real(8) :: dsmalc = 10.0D0

  real(8) , dimension(10) :: fixed_cld_profile
  real(8) , dimension(10) :: cld_profile
  real(8) , dimension(10) :: rnum
  real(8) , dimension(nk) :: rcldfra
  real(8) , dimension(nk) :: old_rcldfra
  real(8) , dimension(nk) :: plevs
  real(8) , dimension(nk) :: sigma

  integer , parameter :: gun = 100
  integer :: nseed
  character(len=256) :: fname , gname
  logical :: addnoise = .false.
  real :: cputime
  integer , dimension(:) , allocatable :: iseed

  if ( addnoise ) then
    ! Initialize a random number generator
    call random_seed(size=nseed)
    call cpu_time(cputime)
    allocate(iseed(nseed))
    iseed = int(cputime) + 37*(/(k-1,k=1,nseed)/)
    call random_seed(put=iseed)
  end if

  fixed_cld_profile(1)  = 0.130D0
  fixed_cld_profile(2)  = 0.125D0
  fixed_cld_profile(3)  = 0.120D0
  fixed_cld_profile(4)  = 0.080D0
  fixed_cld_profile(5)  = 0.080D0
  fixed_cld_profile(6)  = 0.080D0
  fixed_cld_profile(7)  = 0.085D0
  fixed_cld_profile(8)  = 0.085D0
  fixed_cld_profile(9)  = 0.105D0
  fixed_cld_profile(10) = 0.110D0

  plevs(nk) = 963.0D0
  call getsigma(sigma)
  do k = 1 , nk
    plevs(k) = sigma(k)*plevs(nk) + 50.0D0
  end do

  do ids = 0 , nds

    write(fname,'(a,i0.4,a)') 'profile_', int(ds), '.nc'
    istat = nf90_create(fname,nf90_clobber,ncid)
    istat = nf90_def_dim(ncid,'klev',nk,idims(1))
    istat = nf90_def_dim(ncid,'record',nf90_unlimited,idims(2))
    istat = nf90_def_var(ncid,'klev',nf90_float,idims(1:1),ivars(1))
    istat = nf90_def_var(ncid,'cld',nf90_float,idims,ivars(2))
    istat = nf90_def_var(ncid,'totcf',nf90_float,idims(2:2),ivars(3))
    istat = nf90_def_var(ncid,'oldcld',nf90_float,idims,ivars(4))
    istat = nf90_def_var(ncid,'oldtotcf',nf90_float,idims(2:2),ivars(5))
    istat = nf90_put_att(ncid,nf90_global,'ds',ds)

    ! This part is in mod_param :
    !  For a given grid size , fix a Maximum Convective Cloud Cover
    if ( ids == 0 ) then
      ds = 10.0D0
    else
      ds = dsinc*dble(ids)
    end if
    write(gname,'(a,i0.4,a)') 'gnuplot_', int(ds), '.dat'
    open(file=gname,form='formatted',unit=gun)
    dxtemc = dmin1(dmax1(ds,dsmalc),dlargc)
    clfrcv = afracl + (afracs-afracl)*((dlargc-dxtemc)/(dlargc-dsmalc))**2
    clfrcv = dmin1(clfrcv,1.0D0)
    istat = nf90_put_att(ncid,nf90_global,'clfrcv',clfrcv)
    istat = nf90_enddef(ncid)

    istat = nf90_put_var(ncid,idims(1),plevs)
    ! Test the cloud generation.
    ! Input are kbot and ktop (maximum cloud top and bottom)
    ! Have them assume all possible ranges.

    irec = 1
    do kbot = nk-2 , 4 , -1
      do ktop = kbot-1 , 3 , -1

        rcldfra(:) = 0.0D0
        old_rcldfra(:) = 0.0D0

        ! The fixed shape for cloud profile can be perturbed
        if ( addnoise ) then
          call random_number(rnum)
          cld_profile = (0.75D0+(rnum/2.0D0))*fixed_cld_profile
        else
          cld_profile = fixed_cld_profile
        end if

        kclth = kbot-ktop
        dp = plevs(kbot) - plevs(ktop)
        scalep = min(dp/maxcloud_dp,1.0D0)
        akclth = 1.0D0/dble(kclth)

        ! This is function of the grid.
        ! In regcm3 this was disabled.
        scalef = (1.0D0-clfrcv)
        ! and this was used:
        scalef = 0.5D0

        do k = ktop , kbot
          ikh = max(1,min(10,int((dble(k-ktop+1)/dble(kclth))*10.0D0)))
          ! The icumcloud = 2
          rcldfra(k) = cld_profile(ikh)*clfrcv*scalep
          ! The old icumcloud = 1
          old_rcldfra(k) = 1.0D0 - scalef**akclth
        end do

        istart(1) = 1
        istart(2) = irec
        icount(1) = nk
        icount(2) = 1
        istat = nf90_put_var(ncid,ivars(2),rcldfra,istart,icount)
        istat = nf90_put_var(ncid,ivars(4),old_rcldfra,istart,icount)
        ! Compute the cloud fraction as seen by the radiative model
        totcf = 1.0D0
        totcf_old = 1.0D0
        do k = 1 , nk
          totcf = totcf * ( 1.0D0 - rcldfra(k) )
          totcf_old = totcf_old * ( 1.0D0 - old_rcldfra(k) )
        end do
        totcf = 1.0D0 - totcf
        write(gun,'(i4,3f11.6)') kclth , dp , totcf , totcf_old
        istart(1) = irec
        icount(1) = 1
        helpme(1) = totcf
        istat = nf90_put_var(ncid,ivars(3),helpme,istart,icount)
        helpme(1) = totcf_old
        istat = nf90_put_var(ncid,ivars(5),helpme,istart,icount)
        irec = irec + 1
      end do
    end do
    istat = nf90_close(ncid)
    close(gun)
  end do

  contains

    subroutine getsigma(sigma)
      implicit none
      real(8) , intent(out) , dimension(:) :: sigma
      integer :: kz
      kz = size(sigma)
      if ( kz == 14 ) then                      ! RegCM2
        sigma(1) = 0.0D0
        sigma(2) = 0.04D0
        sigma(3) = 0.10D0
        sigma(4) = 0.17D0
        sigma(5) = 0.25D0
        sigma(6) = 0.35D0
        sigma(7) = 0.46D0
        sigma(8) = 0.56D0
        sigma(9) = 0.67D0
        sigma(10) = 0.77D0
        sigma(11) = 0.86D0
        sigma(12) = 0.93D0
        sigma(13) = 0.97D0
        sigma(14) = 0.99D0
        sigma(15) = 1.0D0
      else if ( kz == 18 ) then                 ! RegCM3, default
        sigma(1) = 0.0D0
        sigma(2) = 0.05D0
        sigma(3) = 0.10D0
        sigma(4) = 0.16D0
        sigma(5) = 0.23D0
        sigma(6) = 0.31D0
        sigma(7) = 0.39D0
        sigma(8) = 0.47D0
        sigma(9) = 0.55D0
        sigma(10) = 0.63D0
        sigma(11) = 0.71D0
        sigma(12) = 0.78D0
        sigma(13) = 0.84D0
        sigma(14) = 0.89D0
        sigma(15) = 0.93D0
        sigma(16) = 0.96D0
        sigma(17) = 0.98D0
        sigma(18) = 0.99D0
        sigma(19) = 1.0D0
      else if ( kz == 23 ) then                 ! MM5V3
        sigma(1) = 0.0D0
        sigma(2) = 0.05D0
        sigma(3) = 0.1D0
        sigma(4) = 0.15D0
        sigma(5) = 0.2D0
        sigma(6) = 0.25D0
        sigma(7) = 0.3D0
        sigma(8) = 0.35D0
        sigma(9) = 0.4D0
        sigma(10) = 0.45D0
        sigma(11) = 0.5D0
        sigma(12) = 0.55D0
        sigma(13) = 0.6D0
        sigma(14) = 0.65D0
        sigma(15) = 0.7D0
        sigma(16) = 0.75D0
        sigma(17) = 0.8D0
        sigma(18) = 0.85D0
        sigma(19) = 0.89D0
        sigma(20) = 0.93D0
        sigma(21) = 0.96D0
        sigma(22) = 0.98D0
        sigma(23) = 0.99D0
        sigma(24) = 1.0D0
      end if
    end subroutine getsigma

end program cumclouds
