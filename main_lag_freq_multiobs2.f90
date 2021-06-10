include 'header.h'

program cross
  use dyn_lc
  use rebin
  implicit none 

  integer              :: i, j, jj, reb_dim, o, tot_intervals
  character (len=200)  :: name_base
  logical              :: yes_no

  double precision, parameter   :: pi = acos(-1.0)
  double precision              :: bias2

  double precision, allocatable :: rc_freq_obs(:,:,:), ic_freq_obs(:,:,:)
  double precision, allocatable :: pw_freq1_obs(:,:,:), pw_freq2_obs(:,:,:)
  double precision, allocatable :: pw_freq1_ave(:),pw_freq2_ave(:),&
                                   pw2_freq1_ave(:),pw2_freq2_ave(:)
  double precision, allocatable :: rc_freq_ave(:), ic_freq_ave(:), &
                                   rc2_freq_ave(:),ic2_freq_ave(:),&
                                   rc_ic_freq_ave(:)

  double precision, allocatable :: lag_freq(:), std_rc_freq_ave(:),&
                                   std_ic_freq_ave(:), deriv_rc(:),&
                                   deriv_ic(:), error_prop(:), &
                                   coher2_freq(:), err_cohe_lag(:),&
                                   covariance(:), errA_rc_ic(:),&
                                   errA_lag(:)
  
!NOT power of 2 segments
  double precision, allocatable :: re_lc1(:), im_lc1(:), re_lc2(:), im_lc2(:)

  ! call set_seed(seed)

!******************************************************************
    
  ! call execute_command_line('ls')
          
!GET LIGHT CURVEs 

  verbose_merge = .false.
  obs_num = 7
  max_gap_sec_init = 1
  int_len_dim = 1024
  
  call load_lc_lag_freq_multi2()

  ! call run_test()
  
  write(*,*)
  write(*,*)
  write(*,*) '-------------------------------------'
  write(*,*) 'Start the Fourier analysis'
  write(*,*)
  if(.not. allocated(rc_freq_obs )) allocate(rc_freq_obs (int_number_max, int_len_dim / 2, obs_num))
  if(.not. allocated(ic_freq_obs )) allocate(ic_freq_obs (int_number_max, int_len_dim / 2, obs_num))
  if(.not. allocated(pw_freq1_obs)) allocate(pw_freq1_obs(int_number_max, int_len_dim / 2, obs_num))
  if(.not. allocated(pw_freq2_obs)) allocate(pw_freq2_obs(int_number_max, int_len_dim / 2, obs_num))

  do o = 1, obs_num
     ! write(*,*) 'ciao ciao ', check_power2
     if (.not. check_power2) then
       if (.not. allocated(re_lc1)) allocate(re_lc1(int_len_dim / 2))
       if (.not. allocated(im_lc1)) allocate(im_lc1(int_len_dim / 2))
       if (.not. allocated(re_lc2)) allocate(re_lc2(int_len_dim / 2))
       if (.not. allocated(im_lc2)) allocate(im_lc2(int_len_dim / 2))
     endif
      
     do i = 1, int_number_obs(o)
        if (check_power2) then 
           ! write(*,*) 'start FFT of interval ', i 
           
           call ncperiodogram_frac_rms(lc_freq2_obs(i, :, o), lc_freq1_obs(i, :, o) , rc_freq_obs(i, :, o), ic_freq_obs(i, :, o), dt, int_len_dim) !Cross spectrum  
           call   periodogram_frac_rms(lc_freq1_obs(i, :, o), pw_freq1_obs(i, :, o), dt, int_len_dim)
           call   periodogram_frac_rms(lc_freq2_obs(i, :, o), pw_freq2_obs(i, :, o), dt, int_len_dim)


           write(*,*) 'FFT done of interval ', i 
        else
           call FT_not_fast(lc_freq1_obs(i, :, o), re_lc1, im_lc1, int_len_dim)
           call FT_not_fast(lc_freq2_obs(i, :, o), re_lc2, im_lc2, int_len_dim)

           !Calculate Power spectra and cross spectrum of light curves and ref band
           do jj = 1, int_len_dim / 2    
              rc_freq_obs (i, jj, o) = (re_lc2(jj) * re_lc1(jj)) + (im_lc2(jj) * im_lc1(jj))
              ic_freq_obs (i, jj, o) = (im_lc2(jj) * re_lc1(jj)) - (re_lc2(jj) * im_lc1(jj))
              pw_freq1_obs(i, jj, o) = (re_lc1(jj) * re_lc1(jj)) + (im_lc1(jj) * im_lc1(jj)) 
              pw_freq2_obs(i, jj, o) = (re_lc2(jj) * re_lc2(jj)) + (im_lc2(jj) * im_lc2(jj))  
           enddo

        endif
     enddo

     if (allocated(re_lc1)) deallocate(re_lc1)
     if (allocated(im_lc1)) deallocate(im_lc1)
     if (allocated(re_lc2)) deallocate(re_lc2)
     if (allocated(im_lc2)) deallocate(im_lc2)

  enddo

 
  call make_freq_array()
    
  write(*,*) '   Completed Fourier transforms of the segmets', int_number
  write(*,*)
      
  if(.not. allocated(rc_freq_ave)   ) allocate(rc_freq_ave   (int_len_dim / 2))
  if(.not. allocated(ic_freq_ave)   ) allocate(ic_freq_ave   (int_len_dim / 2))
  if(.not. allocated(rc2_freq_ave)  ) allocate(rc2_freq_ave  (int_len_dim / 2))
  if(.not. allocated(ic2_freq_ave)  ) allocate(ic2_freq_ave  (int_len_dim / 2))
  if(.not. allocated(rc_ic_freq_ave)) allocate(rc_ic_freq_ave(int_len_dim / 2))
  if(.not. allocated(pw_freq1_ave)  ) allocate(pw_freq1_ave  (int_len_dim / 2))
  if(.not. allocated(pw_freq2_ave)  ) allocate(pw_freq2_ave  (int_len_dim / 2))
  if(.not. allocated(pw2_freq1_ave) ) allocate(pw2_freq1_ave (int_len_dim / 2))
  if(.not. allocated(pw2_freq2_ave) ) allocate(pw2_freq2_ave (int_len_dim / 2))

  rc_freq_ave    = 0.d0
  ic_freq_ave    = 0.d0
  rc2_freq_ave   = 0.d0
  ic2_freq_ave   = 0.d0
  rc_ic_freq_ave = 0.d0
  pw_freq1_ave   = 0.d0
  pw_freq2_ave   = 0.d0
  pw2_freq1_ave  = 0.d0
  pw2_freq2_ave  = 0.d0

  tot_intervals  = 0

  do o = 1, obs_num

     do jj = 1, int_len_dim / 2
        do i = 1, int_number_obs(o)
           rc_freq_ave   (jj) = rc_freq_ave   (jj) + rc_freq_obs(i, jj, o)
           ic_freq_ave   (jj) = ic_freq_ave   (jj) + ic_freq_obs(i, jj, o)
           rc2_freq_ave  (jj) = rc2_freq_ave  (jj) + (rc_freq_obs(i, jj, o) * rc_freq_obs(i, jj, o))
           ic2_freq_ave  (jj) = ic2_freq_ave  (jj) + (ic_freq_obs(i, jj, o) * ic_freq_obs(i, jj, o))
           rc_ic_freq_ave(jj) = rc_ic_freq_ave(jj) + (rc_freq_obs(i, jj, o) * ic_freq_obs(i, jj, o))
           pw_freq1_ave  (jj) = pw_freq1_ave  (jj) +  pw_freq1_obs(i, jj, o)
           pw2_freq1_ave (jj) = pw2_freq1_ave (jj) + (pw_freq1_obs(i, jj, o) * pw_freq1_obs(i, jj, o))
           pw_freq2_ave  (jj) = pw_freq2_ave  (jj) +  pw_freq2_obs(i, jj, o)
           pw2_freq2_ave (jj) = pw2_freq2_ave (jj) + (pw_freq2_obs(i, jj, o) * pw_freq2_obs(i, jj, o))
        enddo
     enddo

     tot_intervals = tot_intervals + int_number_obs(o)  
  enddo

  rc_freq_ave    = rc_freq_ave    / dble(tot_intervals)
  ic_freq_ave    = ic_freq_ave    / dble(tot_intervals)
  rc2_freq_ave   = rc2_freq_ave   / dble(tot_intervals)
  ic2_freq_ave   = ic2_freq_ave   / dble(tot_intervals)
  rc_ic_freq_ave = rc_ic_freq_ave / dble(tot_intervals)
  pw_freq1_ave   = pw_freq1_ave   / dble(tot_intervals)
  pw_freq2_ave   = pw_freq2_ave   / dble(tot_intervals)
  pw2_freq1_ave  = pw2_freq1_ave  / dble(tot_intervals)
  pw2_freq2_ave  = pw2_freq2_ave  / dble(tot_intervals)

  write(*,*) '   Averaging over segmets completed'
  write(*,*) 

  ! do i = 1, int_number
  !    do jj = 1, int_len_dim / 2
  !       write(10,*)  freq(jj), df, rc_freq(i, jj)
  !       write(11,*)  freq(jj), df, ic_freq(i, jj)
  !    enddo
  !    write(10,*)'no no'
  !    write(11,*)'no no'
  ! enddo

  if(.not. allocated(lag_freq    )) allocate(lag_freq    (int_len_dim / 2)) 
  if(.not. allocated(std_rc_freq_ave  )) allocate(std_rc_freq_ave  (int_len_dim / 2)) 
  if(.not. allocated(std_ic_freq_ave  )) allocate(std_ic_freq_ave  (int_len_dim / 2)) 
  if(.not. allocated(deriv_rc    )) allocate(deriv_rc    (int_len_dim / 2))
  if(.not. allocated(deriv_ic    )) allocate(deriv_ic    (int_len_dim / 2))
  if(.not. allocated(error_prop  )) allocate(error_prop  (int_len_dim / 2))
  if(.not. allocated(coher2_freq )) allocate(coher2_freq (int_len_dim / 2))
  if(.not. allocated(err_cohe_lag)) allocate(err_cohe_lag(int_len_dim / 2))

  if(.not. allocated(covariance  )) allocate(covariance  (int_len_dim / 2))
  if(.not. allocated(errA_rc_ic  )) allocate(errA_rc_ic  (int_len_dim / 2))
  if(.not. allocated(errA_lag    )) allocate(errA_lag    (int_len_dim / 2))


   write(*,*)
   name_base = 'freq_cross_spec_real.dat'  
   write(*,*) ' Name file real vs frequency: ', trim(name_base)
   write(*,*) '   The errors are calculated with the standard error on the mean '
   write(*,*)
   open(11, file=trim(name_base))
   name_base = 'freq_cross_spec_imaginary.dat'  
   write(*,*) ' Name file imaginary vs frequency: ', trim(name_base)
   write(*,*) '   The errors are calculated with the standard error on the mean '
   write(*,*)
   open(12, file=trim(name_base))

   name_base = 'lag_freq_err_prop.dat'  
   write(*,*) ' Name file lag vs frequency: ', trim(name_base)
   write(*,*) '    The errors are calculated with propagation error formula on the standard deviation of real and imaginary part'
   write(*,*)
   open(13, file=trim(name_base))
   name_base = 'lag_freq_err_cohe.dat'  
   write(*,*) ' Name file lag vs frequency: ', trim(name_base)
   write(*,*) '   The errors are calculated with coherence formula'
   write(*,*)
   open(14, file=trim(name_base))
   
   write(11, *) 'skip on'
   write(11, *) 'read serr 1 2'
   write(12, *) 'skip on'
   write(12, *) 'read serr 1 2'
   write(13, *) 'skip on'
   write(13, *) 'read serr 1 2'
   write(14, *) 'skip on'
   write(14, *) 'read serr 1 2'
   
!lag calculation and errors
   ! int_len_dim = int_len_dim_obs(1)
   
   do jj = 1, int_len_dim / 2 

!Errors through standard deviation on real and imaginary part      
      !lag with propagation formula
      std_rc_freq_ave(jj) = sqrt((rc2_freq_ave(jj) - (rc_freq_ave(jj) * rc_freq_ave(jj)) ))  / (sqrt(dble(tot_intervals)))
      std_ic_freq_ave(jj) = sqrt((ic2_freq_ave(jj) - (ic_freq_ave(jj) * ic_freq_ave(jj)) ))  / (sqrt(dble(tot_intervals)))
      covariance(jj) = (rc_ic_freq_ave(jj) - ((rc_freq_ave(jj) * ic_freq_ave(jj))))  / dble(tot_intervals)
      ! covariance(jj) = 0.0 
      deriv_rc(jj) =  (-1.d0 *  ic_freq_ave(jj) / rc_freq_ave(jj)**2) / (1 + (ic_freq_ave(jj) / rc_freq_ave(jj))**2 )
      deriv_ic(jj) =                      (1.d0 / rc_freq_ave(jj))    / (1 + (ic_freq_ave(jj) / rc_freq_ave(jj))**2 )
      error_prop(jj) = sqrt( (deriv_rc(jj)**2 * std_rc_freq_ave(jj)**2) + (deriv_ic(jj)**2 * std_ic_freq_ave(jj)**2) + (2.d0 * deriv_rc(jj) * deriv_ic(jj) * covariance(jj)) ) / (2.d0 * pi * freq(jj))

      
!bias term
      ! if ((obs_freq_bins(jj) * int_number) .lt. 500 ) then 
      !    bias2 = ((pw_fq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(obs_freq_bins(jj) * int_number))
      ! else
      !    bias2 = 0.0 
      ! endif

      bias2 = 0.0 

!Coherence and error on the lag with the coherence formula
      coher2_freq(jj) = (rc_freq_ave(jj)**2 + ic_freq_ave(jj)**2 - bias2) / (pw_freq1_ave(jj) * pw_freq2_ave(jj))

      err_cohe_lag(jj) = sqrt( (1.d0 - coher2_freq(jj)) / (2.d0 * coher2_freq(jj) * dble(tot_intervals) ) ) / (2.d0 * pi * freq(jj)) 

      ! write(10,*) jj, freq(jj), coher2_freq(jj), rc_freq_ave(jj)**2, ic_freq_ave(jj)**2, pw_freq1_ave(jj),  pw_freq2_ave(jj)

!Adam's formula real and imaginary part
         ! errA_rc_ic(jj) = sqrt (pw_freq_ave2(jj) * (pw_freq_ave1(jj, k) - ( (rc_freq_ave(jj)**2 + ic_freq_ave(jj)**2 - bias2) / (pw_freq_ave2(jj, k) - P_noise_ext_ref(k)) ) ) / (2 * real( int_number) ) )

!Adam's formula lag
         ! errA_lag(jj) = sqrt( pw_fq_en_ref(jj, k) * ( (pw_fq_en(jj, k) / (rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2 - bias2) ) - ( 1 / (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  real(obs_freq_bins(jj) * int_number)) ) / (2 * pi * relevant_freq(jj)) 

      write(11,*) freq(jj), df, rc_freq_ave(jj), std_rc_freq_ave(jj) !/ sqrt(real(int_number))
      write(12,*) freq(jj), df, ic_freq_ave(jj), std_ic_freq_ave(jj) !/ sqrt(real(int_number))
      
      lag_freq(jj) = atan2(ic_freq_ave(jj), rc_freq_ave(jj)) / (2.d0 * pi * freq(jj))

      write(13, *) freq(jj), df, lag_freq(jj), error_prop(jj) 
      write(14, *) freq(jj), df, lag_freq(jj), err_cohe_lag(jj)
   enddo
   write(11, *) 'log x on'
   write(12, *) 'log x on'
   write(13, *) 'log x on'
   write(14, *) 'log x on'

   write(11, *) 'no no'
   write(12, *) 'no no'
   write(13, *) 'no no'
   write(14, *) 'no no' 
   close(11)
   close(12)
   close(13)
   close(14)

    write(*,*) '   Lag and error calculation completed'
    write(*,*) 


! -------------- REBIN --------------      
666   continue 
      if (yes_no('   Do you want to rebin the lag freq spectrum?')) then

         call rebin_lag_freq_obs(freq, rc_freq_obs, ic_freq_obs, pw_freq1_obs, pw_freq2_obs, int_number_max, int_number_obs, int_len_dim, obs_num, reb_dim)
         name_base = 'lag_freq_err_reb.qdp'  
         write(*,*) '   The name the rebinned lag vs freq is ', trim(name_base)
         open(71, file = trim(name_base))
         write(71, *) 'skip on'
         write(71, *) 'read serr 1 2 3'
         do j = 1, reb_dim - 1            
            write(71, *) (reb_freq(j + 1) + reb_freq(j)) * 0.5, (reb_freq(j + 1) - reb_freq(j)) * 0.5, reb_lag_freq(j), reb_err_cohe_lag(j), reb_lag_freq(j), reb_err_prop_lag(j)
         enddo
         write(71, *) 'no no'
         write(71, *) 'log x on'
         write(71, *) 'scr white'
         write(71, *) 'lw 5'
         write(71, *) 'la x Frequency [Hz]'
         write(71, *) 'la y Lag [s]'
         write(71, *) 't off'
         close(71)
         
      endif

      write(*,*)
      if (yes_no('   Are you satified with the rebin?')) then
      else
         goto 666
      endif

      write(*,*)
      write(*,*) '   Rebin done'
      write(*,*)
      write(*,*) '   Analysis concluded'

      
    
  end program cross
