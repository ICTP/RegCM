program cumclouds
  implicit none
  integer , parameter :: nk = 23
  integer , parameter :: nds = 8
  real(8) , parameter :: dsinc = 25.0D0
  real(8) , parameter :: clfrcvmax = 12.0D0
  real(8) , parameter :: maxcloud_dp = 105.0D0 ! In cb

  integer :: k , ktop , kbot , ids , kclth , ikh
  real(8) :: dp , ds , dxtemc , clfrcv , scalep , totcf
  real(8) :: afracl = 0.25D0
  real(8) :: afracs = clfrcvmax
  real(8) :: dlargc = 200.0D0
  real(8) :: dsmalc = 10.0D0
  real(8) , dimension(10) :: cld_profile
  real(8) , dimension(nk) :: rcldfra
  integer , parameter :: iun = 99
  character(len=256) :: fname

  cld_profile(1) = 3.0D0/dble(nk)
  cld_profile(2) = 2.0D0/dble(nk)
  do k = 3 , 9
    cld_profile(k) = 1.00D0/dble(nk)
  end do
  cld_profile(10) = 2.0D0/dble(nk)

  do ids = 0 , nds
    rcldfra(:) = 0.0D0
    if ( ids == 0 ) then
      ds = 10.0D0
    else
      ds = dsinc*dble(ids)
    end if
    write(fname,'(a,i0.4,a)') 'cumcloud_', int(ds), '.txt'
    open(file=fname,form='formatted',unit=iun)
    dxtemc = dmin1(dmax1(ds,dsmalc),dlargc)
    clfrcv = afracl + (afracs-afracl)*((dlargc-dxtemc)/(dlargc-dsmalc))**2
    clfrcv = dmin1(clfrcv,1.0D0)
    write(iun,'(a,f11.6)') 'Horizontal Resolution ds = ',ds
    write(iun,*) 'Convective Cloud Cover parameters after resolution scaling'
    write(iun,'(a,f11.6)') '  Maximum Convective Cloud Cover : ',clfrcv
    do kclth = 2 , nk
      kbot = nk
      ktop = nk - kclth + 1
      dp = dble(kclth)*(maxcloud_dp/dble(nk))
      scalep = min(dp/maxcloud_dp,1.0D0)
      do k = ktop , kbot
        ikh = max(1,min(10,int((dble(k-ktop+1)/dble(kclth))*10.0D0)))
        rcldfra(k) = cld_profile(ikh)*clfrcv*scalep
      end do
      totcf = 1.0D0
      do k = 1 , nk
        write(iun,'(i3,f11.6)') k , rcldfra(k)
        totcf = totcf * ( 1.0D0 - rcldfra(k) )
      end do
      totcf = 1.0D0 - totcf
      write(iun,'(a,2f11.6)') 'Total cloud fraction ', dp , totcf
    end do
    close(iun)
  end do

end program cumclouds
