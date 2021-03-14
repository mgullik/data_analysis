include 'header.h'

program cross
  use dyn_lc
  use rebin
  implicit none 

  integer              :: i, j, jj, reb_dim
  character (len=200)  :: filename1, filename2, name_base
  logical              :: yes_no

  double precision, parameter   :: pi = acos(-1.0)
  double precision              :: bias2

  double precision, allocatable :: rc_freq(:,:), ic_freq(:,:)
  double precision, allocatable :: pw_freq1(:,:), pw_freq2(:,:)
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
  
  ! call set_seed(seed)

!******************************************************************
    
  ! call execute_command_line('ls')
          
!GET LIGHT CURVEs 
  filename1 = '/Users/gullo/Work/BHB_project/CygX1_Ole/ni2636010101_0.01s_1000eV-2000eV.lc'
  filename2 = '/Users/gullo/Work/BHB_project/CygX1_Ole/ni2636010101_0.01s_5000eV-8000eV.lc'

  call load_lc_lag_freq(filename1, filename2)

  ! call run_test()

  
  write(*,*)
  write(*,*)
  write(*,*) '-------------------------------------'
  write(*,*) 'Start the Fourier analysis'
  write(*,*)
      if(.not. allocated(rc_freq )) allocate(rc_freq (int_number, int_len_dim / 2))
      if(.not. allocated(ic_freq )) allocate(ic_freq (int_number, int_len_dim / 2))
      if(.not. allocated(pw_freq1)) allocate(pw_freq1(int_number, int_len_dim / 2))
      if(.not. allocated(pw_freq2)) allocate(pw_freq2(int_number, int_len_dim / 2))
       
    do i = 1, int_number
       if (check_power2) then 
          ! call ncperiodogram(lc_freq2(i, :), lc_freq1(i, :) , rc_freq(i, :), ic_freq(i, :), int_len_dim) !Cross spectrum 
          ! call   periodogram(lc_freq1(i, :), pw_freq1(i, :), int_len_dim)
          ! call   periodogram(lc_freq2(i, :), pw_freq2(i, :), int_len_dim)


          call ncperiodogram_frac_rms(lc_freq2(i, :), lc_freq1(i, :) , rc_freq(i, :), ic_freq(i, :), dt, int_len_dim) !Cross spectrum  
          call   periodogram_frac_rms(lc_freq1(i, :), pw_freq1(i, :), dt, int_len_dim)
          call   periodogram_frac_rms(lc_freq2(i, :), pw_freq2(i, :), dt, int_len_dim)
          ! call ncperiodogram_no_norm(lc_freq2(i, :), lc_freq1(i, :) , rc_freq(i, :), ic_freq(i, :), dt, int_len_dim) !Cross spectrum 
          ! call   periodogram_no_norm(lc_freq1(i, :), pw_freq1(i, :), dt, int_len_dim)
          ! call   periodogram_no_norm(lc_freq2(i, :), pw_freq2(i, :), dt, int_len_dim)
         
          ! write(*,*) 'FFT done of interval ', i 
       else
          write(*,*) 'The interval length is not a power of 2'
          exit
       endif
    enddo

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
    
    do jj = 1, int_len_dim / 2
       do i = 1, int_number
          rc_freq_ave   (jj) = rc_freq_ave   (jj) + rc_freq(i, jj)
          ic_freq_ave   (jj) = ic_freq_ave   (jj) + ic_freq(i, jj)
          rc2_freq_ave  (jj) = rc2_freq_ave  (jj) + (rc_freq(i, jj) * rc_freq(i, jj))
          ic2_freq_ave  (jj) = ic2_freq_ave  (jj) + (ic_freq(i, jj) * ic_freq(i, jj))
          rc_ic_freq_ave(jj) = rc_ic_freq_ave(jj) + (rc_freq(i, jj) * ic_freq(i, jj))
          pw_freq1_ave (jj)  = pw_freq1_ave  (jj) + pw_freq1(i, jj)
          pw2_freq1_ave(jj)  = pw2_freq1_ave (jj) + (pw_freq1(i, jj) * pw_freq1(i, jj))
          pw_freq2_ave (jj)  = pw_freq2_ave  (jj) + pw_freq2(i, jj)
          pw2_freq2_ave(jj)  = pw2_freq2_ave (jj) + (pw_freq2(i, jj) * pw_freq2(i, jj))
       enddo
    enddo
    
    rc_freq_ave    = rc_freq_ave    / real(int_number)
    ic_freq_ave    = ic_freq_ave    / real(int_number)
    rc2_freq_ave   = rc2_freq_ave   / real(int_number)
    ic2_freq_ave   = ic2_freq_ave   / real(int_number)
    rc_ic_freq_ave = rc_ic_freq_ave / real(int_number)
    pw_freq1_ave   = pw_freq1_ave   / real(int_number)
    pw_freq2_ave   = pw_freq2_ave   / real(int_number)
    pw2_freq1_ave  = pw2_freq1_ave  / real(int_number)
    pw2_freq2_ave  = pw2_freq2_ave  / real(int_number)

    write(*,*) '   Averaging over segmets completed'
    write(*,*) 

    call make_freq_array()
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
   do jj = 1, int_len_dim / 2 

!Errors through standard deviation on real and imaginary part      
      !lag with propagation formula
      std_rc_freq_ave(jj) = sqrt((rc2_freq_ave(jj) - (rc_freq_ave(jj) * rc_freq_ave(jj)) ))  / (sqrt(dble(int_number)))
      std_ic_freq_ave(jj) = sqrt((ic2_freq_ave(jj) - (ic_freq_ave(jj) * ic_freq_ave(jj)) ))  / (sqrt(dble(int_number)))
      covariance(jj) = (rc_ic_freq_ave(jj) - ((rc_freq_ave(jj) * ic_freq_ave(jj))))  / dble(int_number)
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

      err_cohe_lag(jj) = sqrt( (1.d0 - coher2_freq(jj)) / (2.d0 * coher2_freq(jj) * dble(int_number) ) ) / (2.d0 * pi * freq(jj)) 

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

         call rebin_lag_freq(freq, rc_freq, ic_freq, pw_freq1, pw_freq2, int_number, int_len_dim, reb_dim)
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
