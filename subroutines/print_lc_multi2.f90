subroutine print_lc_multi2(time_obs)
  use dyn_lc
  implicit none

  double precision, intent(IN) ::  time_obs(int_number_max, int_len_dim_max, obs_num)

  
  integer :: o, k, i, j
  character (len=500) :: filename 

  
  
  do o = 1, obs_num
     do k = 1, en_num
        write(filename, '(A,I0,A,I0,A)') 'products/printed_lcs/lc_obs', o,'_en',k,'.dat'
        ! write(*,*) 'writing the lc ', trim(filename)
        open(3, file = trim(filename))
        write(3,*) 'skip on'
        do i = 1, int_number_obs(o)
           do j = 1, int_len_dim_obs(o)
              write(3,*) j, time_obs(i, j, o), lc_en_obs (i, j, k, o)
           enddo
           write(3,*) 'no no'
        enddo
        close(3)
        
     enddo
  enddo
  

end subroutine print_lc_multi2
  
