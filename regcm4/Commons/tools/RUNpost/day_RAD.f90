      program day_rad
      implicit none
!
! Local variables
!
      real(4) , dimension(398,330) :: b
      real(4) , dimension(398,330,81) :: c
      integer :: i , j , mday , month , mrec , n , n5 , nday , nrec ,   &
               & nrecord
      character(14) , dimension(72) :: in_fn
      character(12) , dimension(12) :: outfn
!
      data in_fn/'RAD.1990010100' , 'RAD.1990010600' ,                  &
         & 'RAD.1990011100' , 'RAD.1990011600' , 'RAD.1990012100' ,     &
          &'RAD.1990012600' , 'RAD.1990020100' , 'RAD.1990020600' ,     &
          &'RAD.1990021100' , 'RAD.1990021600' , 'RAD.1990022100' ,     &
          &'RAD.1990022600' , 'RAD.1990030100' , 'RAD.1990030600' ,     &
          &'RAD.1990031100' , 'RAD.1990031600' , 'RAD.1990032100' ,     &
          &'RAD.1990032600' , 'RAD.1990040100' , 'RAD.1990040600' ,     &
          &'RAD.1990041100' , 'RAD.1990041600' , 'RAD.1990042100' ,     &
          &'RAD.1990042600' , 'RAD.1990050100' , 'RAD.1990050600' ,     &
          &'RAD.1990051100' , 'RAD.1990051600' , 'RAD.1990052100' ,     &
          &'RAD.1990052600' , 'RAD.1990060100' , 'RAD.1990060600' ,     &
          &'RAD.1990061100' , 'RAD.1990061600' , 'RAD.1990062100' ,     &
          &'RAD.1990062600' , 'RAD.1990070100' , 'RAD.1990070600' ,     &
          &'RAD.1990071100' , 'RAD.1990071600' , 'RAD.1990072100' ,     &
          &'RAD.1990072600' , 'RAD.1990080100' , 'RAD.1990080600' ,     &
          &'RAD.1990081100' , 'RAD.1990081600' , 'RAD.1990082100' ,     &
          &'RAD.1990082600' , 'RAD.1990090100' , 'RAD.1990090600' ,     &
          &'RAD.1990091100' , 'RAD.1990091600' , 'RAD.1990092100' ,     &
          &'RAD.1990092600' , 'RAD.1990100100' , 'RAD.1990100600' ,     &
          &'RAD.1990101100' , 'RAD.1990101600' , 'RAD.1990102100' ,     &
          &'RAD.1990102600' , 'RAD.1990110100' , 'RAD.1990110600' ,     &
          &'RAD.1990111100' , 'RAD.1990111600' , 'RAD.1990112100' ,     &
          &'RAD.1990112600' , 'RAD.1990120100' , 'RAD.1990120600' ,     &
          &'RAD.1990121100' , 'RAD.1990121600' , 'RAD.1990122100' ,     &
          &'RAD.1990122600'/
      data outfn/'RAD.19900101' , 'RAD.19900201' , 'RAD.19900301' ,     &
          &'RAD.19900401' , 'RAD.19900501' , 'RAD.19900601' ,           &
          &'RAD.19900701' , 'RAD.19900801' , 'RAD.19900901' ,           &
          &'RAD.19901001' , 'RAD.19901101' , 'RAD.19901201'/
      do month = 11 , 12
        open (20,file=outfn(month),form='unformatted',recl=398*330*4,   &
            & access='direct')
        mrec = 0
        do n5 = 1 , 6
          open (10,file=in_fn((month-1)*6+n5),form='unformatted',       &
              & recl=398*330*4,access='direct')
          nrec = 0
          if ( month==1 .and. n5==1 ) nrec = 81
          nrecord = 5
          if ( n5==6 ) then
            if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.&
               & month==8 .or. month==10 .or. month==12 ) then
              nrecord = 6
            else if ( month==4 .or. month==6 .or. month==9 .or.         &
                    & month==11 ) then
              nrecord = 5
            else
!normal       Year     nrecord=3
              nrecord = 3
 
!             Leap  Year     nrecord=4
!             nrecord=4
            end if
          end if
          do nday = 1 , nrecord
            do n = 1 , 81
              do j = 1 , 330
                do i = 1 , 398
                  c(i,j,n) = 0.0
                end do
              end do
            end do
            do mday = 1 , 4
              do n = 1 , 81
                nrec = nrec + 1
                read (10,rec=nrec) b
                do j = 1 , 330
                  do i = 1 , 398
                    c(i,j,n) = c(i,j,n) + b(i,j)
                  end do
                end do
              end do
            end do
            do n = 1 , 81
              do j = 1 , 330
                do i = 1 , 398
                  c(i,j,n) = c(i,j,n)/4.
                end do
              end do
              mrec = mrec + 1
              write (20,rec=mrec) ((c(i,j,n),i=1,398),j=1,330)
            end do
          end do
          close (10)
        end do
        close (20)
      end do
      end program day_rad
