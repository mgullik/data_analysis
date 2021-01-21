!---------------------------------------------------------------------!
    function colnum_f(unit,col_name_temp,status)
! PURPOSE: Get column number with a paricular name col_name_temp
!   This code check if there more than one column with a similar name and warns the user
      implicit none
      integer :: unit,status,colnum_f
      logical :: casesen,logic
      character (len=30) col_name_temp,col_name,error_description
      
      casesen = .false. !not case sensitive (if yes, it is)
      call ftgcnn(unit,casesen,col_name_temp,col_name,colnum_f,status)      
      call ftgerr(status,error_description)
      do while(status .ne. 0)
         if (status .eq. 219) then
            call ftgerr(status,error_description)
            write(*,*) '!! ATTENTION !! ',error_description
            write(*,*) "   Change the name of the column (type 'no' to avoid)"
            read(*,*)  col_name_temp
            if (col_name_temp .eq. 'no') then
               colnum_f = -1
               goto 11
            else
               logic = .false.
               status = 0
            endif
         else        
            write(*,*) '!! ATTENTION !! ', error_description
            write(*,'(A,A8,A)') '   Do you want to use column: ', col_name,'? (T for yes - F for the next one)'
            read(*,*) logic
         endif
      
11       continue
         
         if (logic) then
            status = 0
         else
            call ftgcnn(unit,casesen,col_name_temp,col_name,colnum_f,status)
         end if
         
      enddo
    end function colnum_f
!---------------------------------------------------------------------!
      
