 subroutine make_freq_array()
   use dyn_lc
   implicit none
   integer :: j
  
! !Frequency array
   df = 0.5d0 / (dt * int_len_dim)
   if (.not. allocated(freq)) allocate(freq(int_len_dim / 2))
   ! print *, dt, int_len_dim
   
   do j = 1, int_len_dim / 2 
      freq(j) = real(j) / (dt * int_len_dim)
      ! write(*,*) freq(j)
   enddo

  
 end subroutine make_freq_array

!-------------------------------------------------------------------
 subroutine make_freq_array_obs()
   use dyn_lc
   implicit none
   integer :: j, o
   
   ! !Frequency array

   if(.not. allocated(df_obs)) allocate(df_obs(obs_num))
   if (.not. allocated(freq_obs)) allocate(freq_obs(int_len_dim_max / 2, obs_num))
   do o = 1, obs_num
      df_obs(o) = 0.5d0 / (dt * int_len_dim_obs(o))
      
      do j = 1, int_len_dim_obs(o) / 2 
         freq_obs(j, o) = real(j) / (dt * int_len_dim_obs(o))
      enddo
   enddo
   
 end subroutine make_freq_array_obs
