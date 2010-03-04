      module mod_datewk
      implicit none

      integer , dimension(427+1045) :: wkday

      contains

      subroutine headwk
      implicit none
!
! Local variables
!
      integer :: i , mday , month , myear
!
      wkday(1) = 19811029
      do i = 2 , 427
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      wkday(428) = 19891231
      do i = 429 , 427 + 1045
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      end subroutine headwk

      end module mod_datewk
