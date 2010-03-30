!!>
!!c  ROUTINE: HEADER
!!c  ACTION: print out initial information
!!c          about machine/time,compilation stamp 
!!c          and compilation flags 
!!c  
!!c  
!!************************************************************
!!<
!#ifdef INTEL
!  include 'ifport.f90'
!#endif

SUBROUTINE header(myid)

  IMPLICIT NONE 
  !! local variables:
  INTEGER len,len_strim
  INTEGER ihost,idir
!#ifdef INTEL
!  INTEGER hostnm
!#else
!  INTEGER gethostname 
!#endif
  INTEGER getcwd
  CHARACTER (len=24) :: data='?'
  CHARACTER (len=32) :: hostname='?' 
  CHARACTER (len=30) :: user='?' 
  CHARACTER (len=100) :: directory='?'
  INTEGER :: nrite=6, myid
  
  if (myid.eq.1)  then 

     !!open here the output file:
!     OPEN ( 20,file='OUTPUT')
     WRITE (nrite,"(/,2x,'This is RegCM version 4 ')")
     WRITE (nrite,100)  SVN_REV, __DATE__ , __TIME__   
100  FORMAT(2x,' SVN Revision: ',A,' compiled at: data : ',A,'  time: ',A,/)

#ifdef IBM
     hostname='ibm platform '
     user= 'Unknown'
     call fdate_(data)
#else
     Ihost = hostnm(hostname)
     call getlog(user)
     CALL fdate(data)
#endif 

     Idir=GETCWD(directory)


     WRITE(nrite,*) ": this run start at    : ",data
     len=len_strim(user)
     WRITE(nrite,*) ": it is submitted by   : ",user(1:len)
     len=len_strim(hostname)
     WRITE(nrite,*) ": it is running on     : ",hostname(1:len-1)
     len=len_strim(directory) 

     WRITE(nrite,*) ": in directory         : ",directory(1:len)
     WRITE(nrite,*) "                      " 
     end if 
  RETURN 
END SUBROUTINE header

!!>
  !!   ROUTINE : LEN_STRIM
  !!   PACKAGE VERSION : DLPROTEIN-2.1
  !!   ACTION : trimming of string
  !!< 
  FUNCTION len_strim (string) RESULT (len_trim_RESULT)
    IMPLICIT NONE 
    CHARACTER  (len=*), INTENT(IN) :: string
    INTEGER :: len_trim_RESULT,k

    len_trim_RESULT =0 
    DO k= LEN(string),1,-1
       IF (string(k:k) /=' ') THEN
          len_trim_RESULT =k
          EXIT
       END IF
    END DO

  END FUNCTION len_strim
