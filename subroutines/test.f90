subroutine run_test()
  implicit none
  write(*,*) 
  write(*,*) '****************************************************'
  write(*,*) 
  write(*,*) '   ----- STARTING TEST -----'
  write(*,*) 

  call test_lc_parameters()
  call test_print_lcs()

  write(*,*) 
  write(*,*) '   ----- ENDING TEST -----'
  write(*,*) 
  write(*,*) '****************************************************'
  write(*,*) 

end subroutine run_test


subroutine test_lc_parameters()
  use dyn_lc
  implicit none
  character (len=200) :: filename 

  filename = 'tests/test_parameters.txt'
  open(70, file = filename)

  write(70, *) 'light curve dimension', dim_lc
  write(70, *) 'number of GTI intervals ', dim_GTI
  write(70, *) 'number of intervals',  int_number 
  
end subroutine test_lc_parameters

subroutine test_print_lcs()
  use dyn_lc
  implicit none

  integer :: i, j
  real    :: mean_lc1, mean_lc2
  character (len=200) :: filename1, filename2 
  

  filename1 = 'tests/test_print_lc1.qdp'
  open(71, file = filename1)
  filename2 = 'tests/test_print_lc2.qdp'
  open(72, file = filename2)

   write(71,*) 'skip on'
   write(72,*) 'skip on'
  mean_lc1 = 0.0
  mean_lc2 = 0.0
  do i = 1, int_number
     do j = 1, int_len_dim
        write(71, *) time_int(i, j), lc_freq1(i, j)
        write(72, *) time_int(i, j), lc_freq2(i, j)  
        mean_lc1 = mean_lc1 + lc_freq1(i, j)
        mean_lc2 = mean_lc2 + lc_freq2(i, j)
     enddo
     write(71, *) 'no no'
     write(72, *) 'no no'
  enddo
  close(71)
  close(72)

  mean_lc1 = mean_lc1 / real(int_number * int_len_dim)
  mean_lc2 = mean_lc2 / real(int_number * int_len_dim)

  write(*,*) 'Light curve written on files: ', trim(filename1), '; ' , trim(filename2)
  write(*,*) 

  
  do i = 1, int_number
     do j = 1, int_len_dim
        if (abs(lc_freq1(i, j) -  mean_lc1) .gt. 10 * mean_lc1 ) then
           write(*,*) '   Very high value in the light curve 1: ', i, j, lc_freq1(i, j)
           write(*,*)
        endif
        if (abs(lc_freq2(i, j) -  mean_lc2) .gt. 10 * mean_lc2 ) then
           write(*,*) '   Very high value in the light curve 2: ', i, j, lc_freq2(i, j)
           write(*,*)
        endif

        if (lc_freq1(i, j) .gt. 1e5 ) then
           write(*,*) '   Very high value in the light curve 1: ', i, j, lc_freq1(i, j)
           write(*,*)
        endif
        if (lc_freq2(i, j) .gt. 1e5 ) then
           write(*,*) '   Very high value in the light curve 2: ', i, j, lc_freq2(i, j)
           write(*,*)
        endif

     enddo
  enddo
  
end subroutine test_print_lcs
