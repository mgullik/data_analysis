!-----------------------------------------------------------------!
     function get_p_noise_pw(pw, dim, freq, freq_limit)
       implicit none 
       integer          , intent(IN) :: dim
       double precision , intent(IN) :: pw(dim), freq(dim), freq_limit
       double precision              :: get_p_noise_pw

       integer                       :: i, count
       double precision              :: mean 

       count = 0
       mean = 0.d0
       do i = 1, dim
          if (freq(i) .ge. freq_limit) then 
             mean = mean + pw(i)
             count = count + 1
          endif 
       enddo
       ! write(*,*) 'p noise: sum, count', mean, count
       
       get_p_noise_pw = mean / real(count)
     end function get_p_noise_pw
!---------------------------------------------------------------!

