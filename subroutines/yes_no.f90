!----------------------------------------------------------------------- 
     function yes_no(string)
!This function return a logical T or F asking a yes or no question
       implicit none 
       character (*), intent(IN) :: string
       character (len = 10)      :: answer
       character                 :: y_n, y, n 
       logical                   :: yes_no

       write(*,*) trim(string), ' [y/n]'
       read(*,*)  answer

       do 
          y_n = answer(1 : 1)
          y = 'y'
          n = 'n'
       
          if (y_n .eq. y) then 
             yes_no = .true.
             exit
          else if (y_n .eq. n) then 
             yes_no = .false.
             exit
          else
             write(*,*)
             write(*,*) '    PLAESE, answer yes or no'
             read(*,*)  answer
          endif
       enddo

     end function yes_no
!----------------------------------------------------------------------- 
