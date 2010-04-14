      program day_srf2
      implicit none
!
! Local variables
!
      real(4) , dimension(398,330,27) :: b , c
      integer :: i , j , mday , month , mrec , n , n5 , nday , nrec ,   &
               & nrecord
      character(14) , dimension(72) :: in_fn
      character(12) , dimension(12) :: outfn
!
      data in_fn/'SRF.1990010100' , 'SRF.1990010600' ,                  &
         & 'SRF.1990011100' , 'SRF.1990011600' , 'SRF.1990012100' ,     &
          &'SRF.1990012600' , 'SRF.1990020100' , 'SRF.1990020600' ,     &
          &'SRF.1990021100' , 'SRF.1990021600' , 'SRF.1990022100' ,     &
          &'SRF.1990022600' , 'SRF.1990030100' , 'SRF.1990030600' ,     &
          &'SRF.1990031100' , 'SRF.1990031600' , 'SRF.1990032100' ,     &
          &'SRF.1990032600' , 'SRF.1990040100' , 'SRF.1990040600' ,     &
          &'SRF.1990041100' , 'SRF.1990041600' , 'SRF.1990042100' ,     &
          &'SRF.1990042600' , 'SRF.1990050100' , 'SRF.1990050600' ,     &
          &'SRF.1990051100' , 'SRF.1990051600' , 'SRF.1990052100' ,     &
          &'SRF.1990052600' , 'SRF.1990060100' , 'SRF.1990060600' ,     &
          &'SRF.1990061100' , 'SRF.1990061600' , 'SRF.1990062100' ,     &
          &'SRF.1990062600' , 'SRF.1990070100' , 'SRF.1990070600' ,     &
          &'SRF.1990071100' , 'SRF.1990071600' , 'SRF.1990072100' ,     &
          &'SRF.1990072600' , 'SRF.1990080100' , 'SRF.1990080600' ,     &
          &'SRF.1990081100' , 'SRF.1990081600' , 'SRF.1990082100' ,     &
          &'SRF.1990082600' , 'SRF.1990090100' , 'SRF.1990090600' ,     &
          &'SRF.1990091100' , 'SRF.1990091600' , 'SRF.1990092100' ,     &
          &'SRF.1990092600' , 'SRF.1990100100' , 'SRF.1990100600' ,     &
          &'SRF.1990101100' , 'SRF.1990101600' , 'SRF.1990102100' ,     &
          &'SRF.1990102600' , 'SRF.1990110100' , 'SRF.1990110600' ,     &
          &'SRF.1990111100' , 'SRF.1990111600' , 'SRF.1990112100' ,     &
          &'SRF.1990112600' , 'SRF.1990120100' , 'SRF.1990120600' ,     &
          &'SRF.1990121100' , 'SRF.1990121600' , 'SRF.1990122100' ,     &
          &'SRF.1990122600'/
      data outfn/'SRF.19900101' , 'SRF.19900201' , 'SRF.19900301' ,     &
          &'SRF.19900401' , 'SRF.19900501' , 'SRF.19900601' ,           &
          &'SRF.19900701' , 'SRF.19900801' , 'SRF.19900901' ,           &
          &'SRF.19901001' , 'SRF.19901101' , 'SRF.19901201'/
      do month = 11 , 12
        open (20,file=outfn(month),form='unformatted',recl=398*330*4,   &
            & access='direct')
        mrec = 0
        do n5 = 1 , 6
          open (10,file=in_fn((month-1)*6+n5),form='unformatted',       &
              & recl=398*330*4,access='direct')
          nrec = 0
          if ( month==1 .and. n5==1 ) nrec = 27
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
            do n = 1 , 21
              do j = 1 , 330
                do i = 1 , 398
                  c(i,j,n) = 0.0
                end do
              end do
            end do
            do j = 1 , 330
              do i = 1 , 398
                c(i,j,22) = -1.E20
                c(i,j,23) = 1.E20
                c(i,j,24) = -1.E20
                c(i,j,25) = 1.E20
                c(i,j,26) = -1.E20
                c(i,j,27) = 1.E20
              end do
            end do
            do mday = 1 , 24
              do n = 1 , 27
                nrec = nrec + 1
                read (10,rec=nrec) ((b(i,j,n),i=1,398),j=1,330)
              end do
              do n = 1 , 21
                do j = 1 , 330
                  do i = 1 , 398
                    c(i,j,n) = c(i,j,n) + b(i,j,n)
                  end do
                end do
              end do
              do j = 1 , 330
                do i = 1 , 398
                  if ( c(i,j,22)<b(i,j,22) ) c(i,j,22) = b(i,j,22)
                  if ( c(i,j,23)>b(i,j,23) ) c(i,j,23) = b(i,j,23)
                  if ( c(i,j,24)<b(i,j,24) ) c(i,j,24) = b(i,j,24)
                  if ( c(i,j,25)>b(i,j,25) ) c(i,j,25) = b(i,j,25)
                  if ( c(i,j,26)<b(i,j,26) ) c(i,j,26) = b(i,j,26)
                  if ( c(i,j,27)>b(i,j,27) ) c(i,j,27) = b(i,j,27)
                end do
              end do
            end do
            do n = 1 , 21
              do j = 1 , 330
                do i = 1 , 398
                  c(i,j,n) = c(i,j,n)/24.
                end do
              end do
            end do
            do n = 1 , 27
              mrec = mrec + 1
              write (20,rec=mrec) ((c(i,j,n),i=1,398),j=1,330)
            end do
          end do
          close (10)
        end do
        close (20)
      end do
      end program day_srf2
