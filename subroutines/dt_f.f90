!---------------------------------------------------------------------!
    function dt_f(unit,keyword,status)
      implicit none
      integer :: unit,status
      character (len=30) keyword,comment,error_description
      double precision dt_f

      call ftgkyd(unit, keyword, dt_f, comment, status)

      if(status .ne. 0) then 
         call ftgerr(status,error_description)
         if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
         write(*,*) '   The keyword does not match, try another one'
         read(*,*) keyword
         status = 0 
         call ftgkyd(unit,keyword,dt_f,comment,status)
         if(status .ne. 0) then 
            call ftgerr(status,error_description)
            if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
            write(*,*) '   Neither this one, EXIT...'
         endif
      endif
    end function dt_f
!---------------------------------------------------------------------!
