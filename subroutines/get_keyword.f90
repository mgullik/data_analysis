! !---------------------------------------------------------------------!
!     function get_keyword_string(unit,keyword,status)
!       implicit none
!       integer :: unit,status
!       character (len=30) keyword,comment,error_description
!       character (len=100) get_keyword_string

!       call ftgkyd(unit, keyword, get_keyword_string, comment, status)

!       if(status .ne. 0) then 
!          call ftgerr(status,error_description)
!          if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
!          write(*,'(A,A,A)') '   The keyword: ', trim(keyword),' does not match any of the keyword in the extension, try another one'
!          read(*,*) keyword
!          status = 0 
!          call ftgkyd(unit,keyword,get_keyword_string,comment,status)
!          if(status .ne. 0) then 
!             call ftgerr(status,error_description)
!             if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
!             write(*,*) '   Neither this one, EXIT...'
!          endif
!       endif
!     end function get_keyword_string
! !---------------------------------------------------------------------!

!---------------------------------------------------------------------!
    function get_keyword_double(unit,keyword,status)
      implicit none
      integer :: unit,status
      character (len=30) keyword,comment,error_description
      double precision get_keyword_double

      call ftgkyd(unit, keyword, get_keyword_double, comment, status)

      if(status .ne. 0) then 
         call ftgerr(status,error_description)
         if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
         write(*,'(A,A,A)') '   The keyword: ', trim(keyword),' does not match any of the keyword in the extension, try another one'
         read(*,*) keyword
         status = 0 
         call ftgkyd(unit,keyword,get_keyword_double,comment,status)
         if(status .ne. 0) then 
            call ftgerr(status,error_description)
            if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
            write(*,*) '   Neither this one, EXIT...'
         endif
      endif
    end function get_keyword_double
!---------------------------------------------------------------------!
