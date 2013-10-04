program cumclouds
  implicit none
  integer , parameter :: nk = 23
  integer , parameter :: nds = 8
  integer , parameter :: maxtopcloud = 6
  real(8) , parameter :: dsinc = 25.0D0
  real(8) , parameter :: clfrcvmax = 1.0D0
  real(8) , parameter :: maxcloud_dp = 45.0D0 ! In cb
  real(8) , parameter :: deltap = 5.0

  integer :: k , ktop , kbot , ids , kclth , ikh
  real(8) :: dp , ds , dxtemc , clfrcv , scalep , totcf , akclth
  real(8) :: scalef , totcf_old
  real(8) :: afracl = 0.250D0
  real(8) :: afracs = clfrcvmax
  real(8) :: dlargc = 200.0D0
  real(8) :: dsmalc = 10.0D0
  real(8) , dimension(10) :: fixed_cld_profile
  real(8) , dimension(10) :: cld_profile
  real(8) , dimension(10) :: rnum
  real(8) , dimension(nk) :: rcldfra
  real(8) , dimension(nk) :: old_rcldfra
  integer , parameter :: iun = 99 , gun = 100
  integer , dimension(12) :: tval
  integer :: issize
  character(len=256) :: fname , gname
  logical :: addnoise = .false.

  call date_and_time(values=tval(1:8))
  tval(8:12) = tval(4:8)
  issize = 10
  call random_seed(put=tval)

  fixed_cld_profile(1)  = 0.195D0
  fixed_cld_profile(2)  = 0.175D0
  fixed_cld_profile(3)  = 0.155D0
  fixed_cld_profile(4)  = 0.105D0
  fixed_cld_profile(5)  = 0.085D0
  fixed_cld_profile(6)  = 0.075D0
  fixed_cld_profile(7)  = 0.065D0
  fixed_cld_profile(8)  = 0.055D0
  fixed_cld_profile(9)  = 0.045D0
  fixed_cld_profile(10) = 0.045D0

  do ids = 0 , nds
    if ( ids == 0 ) then
      ds = 10.0D0
    else
      ds = dsinc*dble(ids)
    end if
    write(fname,'(a,i0.4,a)') 'profile_', int(ds), '.txt'
    write(gname,'(a,i0.4,a)') 'gnuplot_', int(ds), '.dat'
    open(file=fname,form='formatted',unit=iun)
    open(file=gname,form='formatted',unit=gun)
    dxtemc = dmin1(dmax1(ds,dsmalc),dlargc)
    clfrcv = afracl + (afracs-afracl)*((dlargc-dxtemc)/(dlargc-dsmalc))**2
    clfrcv = dmin1(clfrcv,1.0D0)
    write(iun,'(a,f11.6)') 'Horizontal Resolution ds = ',ds
    write(iun,*) 'Convective Cloud Cover parameters after resolution scaling'
    write(iun,'(a,f11.6)') '  Maximum Convective Cloud Cover : ',clfrcv
    do ktop = maxtopcloud , nk-3
      do kbot = ktop+1 , nk-2
        rcldfra(:) = 0.0D0
        old_rcldfra(:) = 0.0D0
        if ( addnoise ) then
          call random_number(rnum)
          cld_profile = (0.75D0+(rnum/2.0D0))*fixed_cld_profile
        else
          cld_profile = fixed_cld_profile
        end if
        kclth = kbot-ktop
        dp = dble(kclth)*deltap
        scalep = min(dp/maxcloud_dp,1.0D0)
        akclth = 1.0D0/dble(kclth)
        scalef = (1.0D0-clfrcv)
        scalef = 0.50
        do k = ktop , kbot
          ikh = max(1,min(10,int((dble(k-ktop+1)/dble(kclth))*10.0D0)))
          rcldfra(k) = cld_profile(ikh)*clfrcv*scalep
          old_rcldfra(k) = 1.0D0 - scalef**akclth
        end do
        totcf = 1.0D0
        totcf_old = 1.0D0
        do k = 1 , nk
          write(iun,'(i3,2f11.6)') k , rcldfra(k) , old_rcldfra(k)
          totcf = totcf * ( 1.0D0 - rcldfra(k) )
          totcf_old = totcf_old * ( 1.0D0 - old_rcldfra(k) )
        end do
        totcf = 1.0D0 - totcf
        write(iun,'(a,3f11.6)') 'Total cloud fraction ', dp , totcf , totcf_old
        write(gun,'(i4,3f11.6)') kclth , dp , totcf , totcf_old
      end do
    end do
    close(iun)
    close(gun)
  end do

end program cumclouds
