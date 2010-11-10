      integer, parameter :: GPTLwall           = 1
      integer, parameter :: GPTLcpu            = 2
      integer, parameter :: GPTLabort_on_error = 3
      integer, parameter :: GPTLoverhead       = 4
      integer, parameter :: GPTLdepthlimit     = 5

      integer, parameter :: GPTLgettimeofday   = 6
      integer, parameter :: GPTLnanotime       = 7
      integer, parameter :: GPTLrtc            = 8
      integer, parameter :: GPTLmpiwtime       = 9
      integer, parameter :: GPTLclockgettime   = 10
      integer, parameter :: GPTLpapitime       = 11

! Externals

      integer gptlsetoption
      integer gptlinitialize
      integer gptlstart
      integer gptlstop
      integer gptlstamp 
      integer gptlpr
      integer gptlreset 
      integer gptlfinalize
      integer gptlprint_memusage
      integer gptlget_memusage
      integer gptlsetutr
      integer gptlenable
      integer gptldisable
      integer gptlquery
      integer gptl_papiprinttable
      integer gptl_papiname2id

      external gptlsetoption
      external gptlinitialize
      external gptlstart
      external gptlstop
      external gptlstamp 
      external gptlpr
      external gptlreset 
      external gptlfinalize
      external gptlprint_memusage
      external gptlget_memusage
      external gptlsetutr
      external gptlenable
      external gptldisable
      external gptlquery
      external gptl_papiprinttable
      external gptl_papiname2id

