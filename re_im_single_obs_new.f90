MODULE dyn_lc
!---------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!---------------------------------------------------------------------
! int_len_dim: length of the interval (in units of element) -> this is decided by the user (must be a power of 2)
! int_number: number of interval in the light curve -> this is calculated automatically
  implicit none 
  integer              :: dim_lc, dim_GTI, int_number = -1, int_len_dim, gap = -1, check_gap_num = -1, en_num
  ! integer, parameter   :: int_len_dim_max = 2000000
  double precision      :: dt = -1.d0
  logical              :: check_power2
  double precision   , allocatable :: time(:), start_GTI(:), end_GTI(:)
  real               , allocatable :: lc(:), err_rate(:), bkg(:)
  double precision   , allocatable :: time_int(:,:)
  real               , allocatable :: lc_int(:,:), bkg_int(:,:)
  real   , allocatable ::  lc_en(:,:,:), bkg_en(:,:,:), lc_ref(:,:)
  ! real   , allocatable :: time_int_o(:,:,:), lc_en_o(:,:,:,:), bkg_en_o(:,:,:,:)
  integer, allocatable :: split_ind(:)
END MODULE dyn_lc
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
!-----------------------------------------------------------------------
module rand_fun

  implicit none 
  private
  integer  :: idum 

  public   :: set_seed, poidev, ran1

  contains 
    
    subroutine set_seed(seed)
      implicit none 
      integer :: seed
      idum = seed
    end subroutine set_seed

!-----------------------------------------------------------------------
   function poidev(xm)
     implicit none
     double precision  poidev, xm, pi
     parameter (pi=3.141592654)
     real alxm, em, g, oldm, sq, t, y!, gammln, ran1
     save alxm, g, oldm, sq
      data oldm /-1./

      if (xm .lt. 12) then
        if(xm .ne. oldm) then
          oldm = xm
          g    = exp(-xm)
        end if
        em = -1.d0
        t  = 1.d0
 2      em = em + 1.d0
        t  = t * ran1()
        if( t .gt. g) goto 2
      else
        if(xm .ne. oldm)then
          oldm = xm
          sq   = sqrt(2. * xm)
          alxm = log(xm)
          g    = xm * alxm - gammln( xm + 1.d0)
        end if
 1      y = tan(pi * ran1())
        em = sq * y + xm
        if(em .lt. 0.) goto 1
        em = int(em)
        t = 0.9d0*(1. + y**2.) * exp(em * alxm - gammln(em + 1.d0) - g)
        if(ran1() .gt. t) goto 1
      end if
      poidev = em
      return
    end function poidev
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      FUNCTION gammln(xx)
        implicit none 
        double precision  gammln, xx
        INTEGER j
        DOUBLE PRECISION ser, stp, tmp, x, y, cof(6)
        SAVE cof,stp
        DATA cof, stp / 76.18009172947146d0, -86.50532032941677d0, 24.01409824083091d0, &
             -1.231739572450155d0, 0.1208650973866179d-2, -0.5395239384953d-5, 2.5066282746310005d0/
        x   = xx
        y   = x
        tmp = x + 5.5d0
        tmp = (x + 0.5d0) * log(tmp) - tmp
        ser = 1.000000000190015d0
        do 11 j = 1, 6
           y   = y + 1.d0
           ser = ser +cof(j) / y
11         continue
           gammln = real(tmp + log(stp * ser / x))
           return
        END function gammln
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        FUNCTION ran1()
          implicit none 
          INTEGER IA,IM,IQ,IR,NTAB,NDIV
          REAL ran1,AM,EPS,RNMX
          PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
          NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
          INTEGER j,k,iv(NTAB),iy
          SAVE iv,iy
          DATA iv /NTAB*0/, iy /0/
          if (idum.le.0.or.iy.eq.0) then
             idum=max(-idum,1)
             do 11 j=NTAB+8,1,-1
                k=idum/IQ
                idum=IA*(idum-k*IQ)-IR*k
                if (idum.lt.0) idum=idum+IM
                if (j.le.NTAB) iv(j)=idum
11              continue
                iy=iv(1)
             endif
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum.lt.0) idum=idum+IM
             j=1+iy/NDIV
             iy=iv(j)
             iv(j)=idum
             ran1=min(AM*iy,RNMX)
             return
          END function ran1
!-----------------------------------------------------------------------

  end module rand_fun
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
program cross
  use dyn_lc
  use rand_fun
  implicit none 

  integer               :: i, j, jj, k, kk, file_line_num
  character (len = 200) :: filename,  name_suffix, name_base, name_base_2,&
                           name_path, name_extension, filename_resp
  character (len = 400) :: filename_ch
  integer, parameter    :: seed = -376736
  real   , parameter    :: pi = acos(-1.0)

!General spectral timing analysis
  real                 :: p_noise_pw
  real   , allocatable :: pw(:), tot_rate_int(:), tot_re(:), tot_im(:)
  logical              :: yes_no

!Rebinning 
  ! integer              :: dim_rebin
  ! real                 :: bin_factor
  ! integer, allocatable :: rebin(:)
  ! real   , allocatable :: freq_rebin(:), pw_rebin(:), err_pw_rebin(:)
  

!LAG vs ENERGY
  real                 :: channel_energy_conversion
  integer, allocatable :: l_bin_ch(:), r_bin_ch(:)
  real   , allocatable :: l_bin(:), r_bin(:),  &
                          re_lc(:), im_lc(:), re_ref(:), im_ref(:), &
                          rc(:), ic(:), pw_ref(:), lc_ref_cross(:), lag(:), &
                          pw_ref_unic(:), pw2_ref_unic(:)
  real   , allocatable :: ave_pw_int(:,:), ave_pw_ref(:,:), pw_ref_int(:,:)

  real                 :: c, ave_rate, df
  real   , allocatable :: ave_rate_en(:), ave_bkg_en(:), freq(:)
  real   , allocatable :: rc_en(:,:,:), ic_en(:,:,:), pw_en(:,:,:), &
                          pw_ref_en(:,:,:)
!P_noise
  real                 :: freq_limit
  real   , allocatable :: P_noise_ext(:), P_noise_ext_ref(:), pw_tot(:)

!multiple freq
  integer              :: freq_num
  integer, allocatable :: lower_fq(:), upper_fq(:), obs_freq_bins(:)
  real                 :: min_freq, max_freq
  real   , allocatable :: fq_min_arr(:), fq_max_arr(:), fq(:),  relevant_freq(:)
  real   , allocatable :: rc_fq_en(:,:), ic_fq_en(:,:), rc2_fq_en(:,:), &
                          ic2_fq_en(:,:), rc_ic_fq_en(:,:), pw_fq_en(:,:), &
                          pw_fq_en_ref(:,:)
!Error lag vs energy
  real                 :: bias2, err_pw_ref_unic
  real   , allocatable :: var_rc(:,:), var_ic(:,:), deriv_rc(:,:), &
                          deriv_ic(:,:), covariance(:,:), &
                          coher2_freq(:,:), error_Aform_rc_ic(:,:), &
                          error_Aform_lag(:,:), error_prop(:,:), &
                          error_cohe_lag(:,:), lag_freq(:,:)

    call set_seed(seed)

!*******************************************************************************************************************************************************************************************************************************************************
!*******************************************!
!*********** LAG-ENERGY SPECTRUM ***********!       
!*******************************************!

      print *, ' '
      print *, ' '
      print *, ' ---------------------------------------------------------------------------'
      print *, ' '
      print *, ' '
   ! if (yes_no('   Do you want to create LAG vs ENERGY spectrum? ')) then       
      call execute_command_line('ls')
      
      print *, ' '
      print *, '   Name of the path!'
      name_path = '/Users/gullo/Work/BHB_project/maxi_j1820/Gullo/ni1200120105'
      write(*,*) name_path
      
      print *, ' '
      print *, '   Name of the channel/energy bin file with path!'
      ! read(*,*)  filename_ch
      filename_ch = '/Users/gullo/Work/BHB_project/maxi_j1820/Gullo/ni1200120105/bins.txt'
      write(*,*) filename_ch
      print *, ' '

      en_num = file_line_num(filename_ch)
      write(*,*) en_num


      if (.not. allocated(l_bin)    ) allocate(l_bin   (en_num))
      if (.not. allocated(r_bin)    ) allocate(r_bin   (en_num))
      if (.not. allocated(l_bin_ch) ) allocate(l_bin_ch(en_num))
      if (.not. allocated(r_bin_ch) ) allocate(r_bin_ch(en_num))

      open(55, file = filename_ch)
         do k = 1, en_num
            read(55, *) l_bin_ch(k), r_bin_ch(k)
            l_bin(k) = real(l_bin_ch(k)) / 1000.0
            r_bin(k) = real(r_bin_ch(k)) / 1000.0
         enddo
      close(55)

      if (.not. allocated(ave_rate_en) ) allocate(ave_rate_en(en_num))
      ! if (.not. allocated(ave_bkg_en)  ) allocate(ave_bkg_en (en_num))
      ave_rate_en = 0.0
      ! ave_bkg_en  = 0.0

!GET ALL THE LIGHT CURVES  

      write(name_base, '(A, A)')  trim(name_path), '/ni1200120105_0mpu7_cl_barycorr_'  
      write(*,*)  name_base
      write(name_extension, '(A)') '.lc'
!      name_extension = '_log_en.lc'

!--------------------------------------------------------------------!
      ! filename = '/Users/gullo/Work/BHB_project/maxi_j1820/Gullo/ni1200120105/ni1200120105_0mpu7_cl_barycorr_50_999.lc'
      ! call extract_lc(filename)
      
      !    if (allocated(lc)) then 
      !       call split_lc()
      !    else
      !       write(*,*) '    !!! NO LIGHT CURVE !!! '
      !       write(*,*) '        STOP HERE   '
      !       stop
      !    endif

      !    do i = 1, int_number
      !       ave_rate = 0.0
      !       do j = 1, int_len_dim
      !          ave_rate = ave_rate + lc_int(i, j)
      !       enddo
      !       write(10,*) time_int(i, 1), ave_rate / real(int_len_dim)
      !    enddo
         
      ! stop
!--------------------------------------------------------------------!
      
      do k = 1, en_num 

         !Work out the name of the light curve created with SAS
         write(filename, '(A, I0, A, I0, A)') trim(name_base), l_bin_ch(k), '_', r_bin_ch(k), trim(name_extension)
         write(*,*)
         write(*,*)
         write(*,*)
         write(*,'(A, A, A)') '    Cross spectrum calculation of  ', trim(filename), '  with the reference band'
         ! write(*,'(A, A, A)') trim(filename)

         call extract_lc(filename)
         if (allocated(lc)) then 
            call split_lc()
         else
            write(*,*) '    !!! NO LIGHT CURVE !!! '
            write(*,*) '        STOP HERE   '
            stop
         endif

         if (.not. allocated(lc_en)   ) allocate(lc_en   (int_number, int_len_dim, en_num))
         ! if (.not. allocated(bkg_en)  ) allocate(bkg_en  (int_number, int_len_dim, en_num))
         if (.not. allocated(time_int)) allocate(time_int(int_number, int_len_dim))

         do j = 1, int_len_dim
            do i = 1, int_number
               lc_en (i, j, k) = lc_int(i, j)
               ! bkg_en(i, j, k) = bkg_int(i, j)
               !For the Poisson noise calculation                
               ! We consider light curves for different energies in a single interval 
               ! On this long single light curve we calculate the Poisson noise
               ave_rate_en(k) = ave_rate_en(k) + lc_int(i, j)
               ! ave_bkg_en(k)  = ave_bkg_en(k)  + bkg_int(i, j)
            enddo
         enddo
         ave_rate_en(k) = ave_rate_en(k)  / real(int_len_dim * int_number)
         ! ave_bkg_en (k) = ave_bkg_en (k)  / real(int_len_dim * int_number)

         write(*,*) 
         write(*,*) '*************************************'
         write(*,*) '   Light curve NAME    '          , trim(filename)
         write(*,*) '   Number of intervals '          , int_number
         write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
         write(*,*) '   Actual exposure (sec) '        , int_len_dim * dt * int_number
         write(*,*) '   Average count rate (count/s)'  , ave_rate_en(k)
         write(*,*) '*************************************'
         write(*,*) 
      enddo

!reference band
      if(.not. allocated(lc_ref)) allocate(lc_ref(int_number, int_len_dim))
      write(filename, '(A, A, A)') trim(name_base), '50_999', trim(name_extension)
      call extract_lc(filename)
      call split_lc()
      do j = 1, int_len_dim
         do i = 1, int_number
            lc_ref(i, j) = lc_int(i, j)
         enddo
      enddo

      
! !Frequency arrays
      df = 0.5 / (dt * int_len_dim)
      if (.not. allocated(freq)) allocate(freq(int_len_dim / 2))
      do j = 1, int_len_dim / 2 
         freq(j) = real(j) / (real(dt * int_len_dim))
         ! write(*,*) freq(j)
      enddo

!Calculating the PDS of the reference band and then print it
      if (.not. allocated(pw_ref_int)) allocate(pw_ref_int(int_number, int_len_dim / 2))
      do i = 1, int_number
         call periodogram_frac_rms(lc_ref(i,:), pw_ref_int(i, :), int_len_dim)
      enddo

    if (.not. allocated(pw_ref_unic)) allocate(pw_ref_unic(int_len_dim/2))
    if (.not. allocated(pw2_ref_unic)) allocate(pw2_ref_unic(int_len_dim/2))
      pw_ref_unic  = 0.0
      pw2_ref_unic = 0.0
      do j = 1, int_len_dim / 2 
         do i = 1, int_number
            pw_ref_unic (j) = pw_ref_unic (j) + pw_ref_int(i, j)
            pw2_ref_unic (j) = pw2_ref_unic (j) + pw_ref_int(i, j)**2
         enddo
      enddo
      pw_ref_unic  = pw_ref_unic  / real(int_number)
      pw2_ref_unic = pw2_ref_unic / real(int_number)
      
      ! name_base_2 = 'ener_PDS_ref_band.dat'  
      ! write(*,*) 'name file of the reference band PDS ', trim(name_base_2)
      ! open(71, file = trim(name_base_2))
      write(71, *) 'skip on'
      write(71, *) 'read serr 1 2 3'
      do j = 1, int_len_dim / 2 - 1
         err_pw_ref_unic = sqrt((pw2_ref_unic(j) - pw_ref_unic(j)**2) / real(int_number))
         write(71, *) (freq(j + 1) + freq(j)) * 0.5 , df, pw_ref_unic(j), err_pw_ref_unic, pw_ref_unic(j), pw_ref_unic(j) / sqrt(real(int_number))  
      enddo
      write(71, *) 'no no'
      write(71, *) 'log x y on'

      freq_limit = 100.0 
      write(*,*) 'Poisson noise reference band', p_noise_pw(pw_ref_unic, int_len_dim / 2, freq, freq_limit )

      
!deallocation for the new observation 
      if(allocated(lc)   ) deallocate(lc  )
      if(allocated(time) ) deallocate(time)
      if(allocated(bkg)  ) deallocate(bkg )

      if(allocated(start_GTI)) deallocate(start_GTI)
      if(allocated(end_GTI)  ) deallocate(end_GTI  )
      if(allocated(split_ind)) deallocate(split_ind)

      if(allocated(lc_int)   ) deallocate(lc_int   )
      ! if(allocated(time_int) ) deallocate(time_int )
      if(allocated(bkg_int)  ) deallocate(bkg_int  )
      gap = -1 
      check_gap_num = -1



!PRINT the total light curve (for each obs if more than one) and the corresponding power spectrum 
!***************************************************************************************************************************!
!       if(yes_no('    Do you want to print the total light curve divided in intervals and its power spectrum? ')) then 

!             allocate(tot_rate_int( int_len_dim * int_number))

!             name_base = 'energy_tot.lc'  
!             write(*,*) 'name file light curve ', trim(name_base)
!             open(70, file = trim(name_base))
            
!             write(70,*) 'skip on'
!             jj = 0
!             do i = 1, int_number
!                do j = 1, int_len_dim
!                   jj = jj + 1 
!                   tot_rate_int(jj) = 0.0
!                   do k = 1, en_num 
!                      tot_rate_int(jj) = tot_rate_int(jj) + lc_en(i, j, k)
! !                     write(*,*)  j, k, lc_en(i, j, k)
!                   enddo
!                   write(70, *) time_int(i, j), tot_rate_int(jj)  
!                enddo
!                write(70, *) 'no no'
!             enddo
!             close(70)
            
!             deallocate(tot_rate_int)

!             ! write(*,*) 'name of the PDS ', trim(name_base)


!       endif
!***************************************************************************************************************************!

      write(*,*) 


!***************************************************************************************************************************!
!Fourier analysis

!Arrays to save the fourier trasform analysis for every energy       
      if(.not. allocated(rc_en)) allocate(rc_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(ic_en)) allocate(ic_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(pw_en)) allocate(pw_en(int_number, int_len_dim / 2, en_num))
      if(.not. allocated(pw_ref_en)) allocate(pw_ref_en(int_number, int_len_dim / 2, en_num))

      print *, ' '
      print *, '   Starting the cross spectrum analysis'

!temporaty arrays for the fourier transform (both FFT and slow one)      
      if (.not. allocated(rc)    ) allocate(rc    (int_len_dim / 2))
      if (.not. allocated(ic)    ) allocate(ic    (int_len_dim / 2))
      if (.not. allocated(pw)    ) allocate(pw    (int_len_dim / 2))
      if (.not. allocated(pw_ref)) allocate(pw_ref(int_len_dim / 2))
      if (.not. allocated(re_lc )) allocate(re_lc (int_len_dim / 2))
      if (.not. allocated(im_lc )) allocate(im_lc (int_len_dim / 2))
      if (.not. allocated(re_ref)) allocate(re_ref(int_len_dim / 2))
      if (.not. allocated(im_ref)) allocate(im_ref(int_len_dim / 2))

      if (.not. allocated(lc_ref_cross)) allocate(lc_ref_cross(int_len_dim))
      
      do k = 1, en_num             
         do i = 1, int_number
            lc_ref_cross = 0.0

!Calculate the reference band for every light curve 
            do kk = 1, en_num 
               do j = 1, int_len_dim
                  if (kk .ne. k) then 
                     lc_ref_cross(j) = lc_ref_cross(j) + lc_en(i, j, kk)
! For the Poisson noise calculation 
                     ! ave_rate_ref(k) =  ave_rate_ref(k) + lc_en(i, j, kk)
                     ! ave_bkg_ref(k)  =  ave_bkg_ref(k)  + bkg_en(i, j, kk)
                  endif
               enddo
            enddo

            
!Print just the refecernce band light curve REMEMBER TO SET en_num = 1
            ! do j = 1, int_len_dim
            !    write(13, *) lc_en(i, j, k)
            ! enddo
            ! write(13, *) 'no no'
            ! do j = 1, int_len_dim
            !    write(14, *) lc_ref_cross(j)
            ! enddo
            ! write(14, *) 'no no'

            if (check_power2) then 
               ! call ncperiodogram(lc_en(i, :, k), lc_ref_cross, rc_en(i, :, k), ic_en(i, :, k), int_len_dim) 
               ! call periodogram(lc_en(i, :, k), pw_en(i, :, k), int_len_dim)
               ! call periodogram(lc_ref_cross, pw_ref_en(i, :, k), int_len_dim)
               call ncperiodogram_frac_rms(lc_en(i, :, k), lc_ref_cross, rc_en(i, :, k), ic_en(i, :, k), int_len_dim) 
               call periodogram_frac_rms(lc_en(i, :, k), pw_en(i, :, k), int_len_dim)
               call periodogram_frac_rms(lc_ref_cross, pw_ref_en(i, :, k), int_len_dim)

               ! call ncperiodogram_frac_rms(lc_en(i, :, k), lc_ref(i,:), rc_en(i, :, k), ic_en(i, :, k), int_len_dim) 
               ! call periodogram_frac_rms(lc_en(i, :, k), pw_en(i, :, k), int_len_dim)
               ! call periodogram_frac_rms(lc_ref(i,:), pw_ref_en(i, :, k), int_len_dim)

               ! write(*,*) 'FFT done', i
            else

!Calculate the Fourier transform of the light curves and reference bands
               call FT_not_fast(lc_en(i, :, k), re_lc, im_lc, int_len_dim)
               call FT_not_fast(lc_ref_cross, re_ref, im_ref, int_len_dim)

!Calculate Power spectra and cross spectrum of light curves and ref band
               do jj = 1, int_len_dim / 2    
                  rc_en    (i, jj, k) = (re_lc(jj) * re_ref(jj)) + (im_lc(jj) * im_ref(jj))
                  ic_en    (i, jj, k) = (im_lc(jj) * re_ref(jj)) - (re_lc(jj) * im_ref(jj))
                  pw_en    (i, jj, k) = (re_lc(jj) * re_lc(jj)) + (im_lc(jj) * im_lc(jj)) 
                  pw_ref_en(i, jj, k) = (re_ref(jj) * re_ref(jj)) + (im_ref(jj) * im_ref(jj))  
               enddo               
               write(*,*) 'Slow Fourier transform done', k, i  
            endif
         enddo

         write(*,'(A, I3)') '   Finished to calculate the cross spectrum number ', k
      enddo
!Deallocate temporary arrays 
      if (allocated(rc)    ) deallocate(rc)
      if (allocated(ic)    ) deallocate(ic)
      if (allocated(pw)    ) deallocate(pw)
      if (allocated(pw_ref)) deallocate(pw_ref)
      if (allocated(re_lc )) deallocate(re_lc )
      if (allocated(im_lc )) deallocate(im_lc )
      if (allocated(re_ref)) deallocate(re_ref)
      if (allocated(im_ref)) deallocate(im_ref)
      if (allocated(lc_ref_cross)) deallocate(lc_ref_cross)     
         
! Print the power spectra for different energy bands for ONE obs (the first one)
   ! write(50,*) 'skip on'
   ! do k = 1, en_num
   !    do j = 1, int_len_dim(1) / 2
   !       write(50,*) freq(j, 1), pw_en(1, j, k, 1)
   !    enddo
   !    write(50, *) 'no no'
   ! enddo
         

!Frequency intervals
   min_freq = 0.1
   ! ! min_freq = 2.9
   max_freq = 30.0
   freq_num = 6
   
   if(.not. allocated(upper_fq)) allocate(upper_fq(freq_num))
   if(.not. allocated(lower_fq)) allocate(lower_fq(freq_num))
      
   if(.not. allocated(fq_min_arr)    ) allocate(fq_min_arr    (freq_num))
   if(.not. allocated(fq_max_arr)    ) allocate(fq_max_arr    (freq_num))
   if(.not. allocated(relevant_freq) ) allocate(relevant_freq (freq_num))
   if(.not. allocated(fq)            ) allocate(fq            (freq_num + 1))
   
!log shape for the frequency ranges 
   call Rlog_scale(min_freq, max_freq, freq_num + 1, fq)

   do jj = 1, freq_num
      fq_min_arr(jj) = fq(jj) 
      fq_max_arr(jj) = fq(jj + 1)
      write(*,*) fq_min_arr(jj), fq_max_arr(jj)
   enddo

   ! fq_min_arr(1) = 0.1
   ! fq_max_arr(1) = 1.0
   ! fq_min_arr(2) = 2.0
   ! fq_max_arr(2) = 15.0
   
    
   do jj = 1, freq_num 
      lower_fq(jj) = ceiling(fq_min_arr(jj) * real(int_len_dim) * dt)
      upper_fq(jj) = floor  (fq_max_arr(jj) * real(int_len_dim) * dt)
      write(*,*)
      write(*,*) "Min frequency bin=", lower_fq(jj)
      write(*,*) "Min true frequency  (Hz)=", lower_fq(jj) / (real(int_len_dim) * dt) 
      write(*,*) "Max frequency bin=", upper_fq(jj)
      write(*,*) "Max true frequency (Hz)=", upper_fq(jj) / (real(int_len_dim) * dt)
      write(*,*) 'Number of frequency bins ', real(upper_fq(jj) - lower_fq(jj) + 1)
      write(*,*) '------------------------------'
      write(*,*)         
   enddo

   if (.not. allocated(obs_freq_bins)) allocate(obs_freq_bins(freq_num))
   if (.not. allocated(relevant_freq)) allocate(relevant_freq(freq_num))
   
   do jj = 1, freq_num       
      obs_freq_bins(jj) = upper_fq(jj) - lower_fq(jj) + 1
      relevant_freq(jj) = (fq_max_arr(jj) + fq_min_arr(jj)) * 0.5 
   enddo
   
   if(.not. allocated(rc_fq_en    )) allocate(rc_fq_en    (freq_num, en_num))
   if(.not. allocated(ic_fq_en    )) allocate(ic_fq_en    (freq_num, en_num))
   if(.not. allocated(rc2_fq_en   )) allocate(rc2_fq_en   (freq_num, en_num))
   if(.not. allocated(ic2_fq_en   )) allocate(ic2_fq_en   (freq_num, en_num))
   if(.not. allocated(rc_ic_fq_en )) allocate(rc_ic_fq_en (freq_num, en_num))
   if(.not. allocated(pw_fq_en    )) allocate(pw_fq_en    (freq_num, en_num))
   if(.not. allocated(pw_fq_en_ref)) allocate(pw_fq_en_ref(freq_num, en_num))
   rc_fq_en     = 0.0
   ic_fq_en     = 0.0
   rc2_fq_en    = 0.0
   ic2_fq_en    = 0.0
   rc_ic_fq_en  = 0.0
   pw_fq_en     = 0.0 
   pw_fq_en_ref = 0.0


   !loop on the observations 
   do k = 1, en_num 
      do jj = 1, freq_num
         do j = lower_fq(jj), upper_fq(jj)
            do i = 1, int_number             
               rc_fq_en    (jj, k) = rc_fq_en    (jj, k) + rc_en(i, j, k)
               ic_fq_en    (jj, k) = ic_fq_en    (jj, k) + ic_en(i, j, k)
               rc2_fq_en   (jj, k) = rc2_fq_en   (jj, k) + rc_en(i, j, k)**2
               ic2_fq_en   (jj, k) = ic2_fq_en   (jj, k) + ic_en(i, j, k)**2
               rc_ic_fq_en (jj, k) = rc_ic_fq_en (jj, k) + rc_en(i, j, k) * ic_en(i, j, k)
               pw_fq_en    (jj, k) = pw_fq_en    (jj, k) + pw_en(i, j, k)
               pw_fq_en_ref(jj, k) = pw_fq_en_ref(jj, k) + pw_ref_en(i, j, k)
            enddo
         enddo
      enddo
   enddo

   if (allocated(rc_en)) deallocate(rc_en)
   if (allocated(ic_en)) deallocate(ic_en)
   write(*,*) 'averaging on intervals and frequency bins; total ', real(obs_freq_bins(1) * int_number)    
   
   do k = 1, en_num 
      do jj = 1, freq_num
         rc_fq_en    (jj, k) = rc_fq_en    (jj, k) / real(obs_freq_bins(jj) * int_number)
         ic_fq_en    (jj, k) = ic_fq_en    (jj, k) / real(obs_freq_bins(jj) * int_number)
         rc2_fq_en   (jj, k) = rc2_fq_en   (jj, k) / real(obs_freq_bins(jj) * int_number)
         ic2_fq_en   (jj, k) = ic2_fq_en   (jj, k) / real(obs_freq_bins(jj) * int_number)
         rc_ic_fq_en (jj, k) = rc_ic_fq_en (jj, k) / real(obs_freq_bins(jj) * int_number)
         pw_fq_en    (jj, k) = pw_fq_en    (jj, k) / real(obs_freq_bins(jj) * int_number)
         pw_fq_en_ref(jj, k) = pw_fq_en_ref(jj, k) / real(obs_freq_bins(jj) * int_number)
      enddo
   enddo


!!! Error calculation !!!

!Compute the poisson noise for every pw of each energy band and ref band for every obs 
!NOTE: it works only for the first interval
   if( .not. allocated(P_noise_ext)    ) allocate(P_noise_ext    (en_num))
   if( .not. allocated(P_noise_ext_ref)) allocate(P_noise_ext_ref(en_num))

   freq_limit = 100
   write(*,*) 'Calculating the Poisson noise for every energy PDS. Above frequency', freq_limit

   do k = 1, en_num
      P_noise_ext    (k) = p_noise_pw(pw_en(1, :, k), int_len_dim / 2, freq, freq_limit )
      P_noise_ext_ref(k) = p_noise_pw(pw_ref_en(1, :, k), int_len_dim / 2, freq, freq_limit )         
   enddo
   if (allocated(pw_en)) deallocate(pw_en)
   if (allocated(pw_ref_en)) deallocate(pw_ref_en)

   
   if(.not. allocated(var_rc           )) allocate(var_rc           (freq_num, en_num)) 
   if(.not. allocated(var_ic           )) allocate(var_ic           (freq_num, en_num)) 
   if(.not. allocated(covariance       )) allocate(covariance       (freq_num, en_num))
   if(.not. allocated(deriv_rc         )) allocate(deriv_rc         (freq_num, en_num))
   if(.not. allocated(deriv_ic         )) allocate(deriv_ic         (freq_num, en_num))
   if(.not. allocated(coher2_freq      )) allocate(coher2_freq      (freq_num, en_num))
   if(.not. allocated(error_Aform_rc_ic)) allocate(error_Aform_rc_ic(freq_num, en_num))
   if(.not. allocated(lag_freq         )) allocate(lag_freq         (freq_num, en_num)) 
   if(.not. allocated(error_prop       )) allocate(error_prop       (freq_num, en_num))
   if(.not. allocated(error_Aform_lag  )) allocate(error_Aform_lag  (freq_num, en_num))
   if(.not. allocated(error_cohe_lag   )) allocate(error_cohe_lag   (freq_num, en_num))

   write(*,*) 'Start error calculation'
!Error std of real and im part, propagation for the lag, Adam's formula for rc and ic and lag  
   do k = 1, en_num 
      do jj = 1, freq_num 

!propagation error 
      var_rc(jj, k) = rc2_fq_en(jj, k) - rc_fq_en(jj, k)**2
      var_ic(jj, k) = ic2_fq_en(jj, k) - ic_fq_en(jj, k)**2

      covariance(jj, k) = rc_ic_fq_en(jj, k) - (rc_fq_en(jj, k) * ic_fq_en(jj, k))

      deriv_rc(jj, k) =  (-1. *  ic_fq_en(jj, k) / rc_fq_en(jj, k)**2) / (1 + (ic_fq_en(jj, k) / rc_fq_en(jj, k))**2 )
      deriv_ic(jj, k) =                      (1. / rc_fq_en(jj, k))    / (1 + (ic_fq_en(jj, k) / rc_fq_en(jj, k))**2 )


!error lag with propagation formula
      error_prop(jj, k) = sqrt( (deriv_rc(jj, k)**2 * var_rc(jj, k) + deriv_ic(jj, k)**2 * var_ic(jj, k) + 2. * deriv_rc(jj, k) * deriv_ic(jj, k) * covariance(jj, k)) / real(obs_freq_bins(jj) * int_number) ) 

!bias term
      if ((obs_freq_bins(jj) * int_number) .gt. 500 ) then 
         bias2 = ((pw_fq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(obs_freq_bins(jj) * int_number))
      else
         bias2 = 0.0 
      endif

      ! bias2 = 0.0 

!Adam's formula real and imaginary part
         error_Aform_rc_ic(jj, k) = sqrt (pw_fq_en_ref(jj, k) * (pw_fq_en(jj, k) - ( (rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2 - bias2) / (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k)) ) ) / (2 * real(obs_freq_bins(jj) * int_number) ) )

!Adam's formula lag
         error_Aform_lag(jj, k) = sqrt( pw_fq_en_ref(jj, k) * ( (pw_fq_en(jj, k) / (rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2 - bias2) ) - ( 1 / (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  real(obs_freq_bins(jj) * int_number)) ) / (2 * pi * relevant_freq(jj)) 

!Coherence error lag 
         coher2_freq(jj, k) = (rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2 - bias2) / (pw_fq_en(jj, k) * pw_fq_en_ref(jj, k))

         error_cohe_lag(jj, k) = sqrt( (1.0 - coher2_freq(jj, k)) / (2.0 * coher2_freq(jj, k) * real(obs_freq_bins(jj) * int_number) ) ) / (2 * pi * relevant_freq(jj)) 

!lag calculation 
         lag_freq(jj, k) = atan2(ic_fq_en(jj, k), rc_fq_en(jj, k)) / (2 * pi * relevant_freq(jj))

         write(*,*) '          k,          jj,     fq_min,      fq_max,       relevant freq,    coherence,       cross^2,              bias^2,            num_freq_average ' 
         write(*,*) k, jj, fq_min_arr(jj), fq_max_arr(jj), relevant_freq(jj), coher2_freq(jj, k), rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2, bias2, obs_freq_bins(jj)

      enddo
   enddo


!printing different stuff
   write(84,  *) 'skip on'
   write(840, *) 'skip on'
   write(841, *) 'skip on'
   write(85,  *) 'skip on'
   write(850, *) 'skip on'
   write(851, *) 'skip on'
   write(86,  *) 'skip on'
   write(860, *) 'skip on'
   write(861, *) 'skip on'
   write(862, *) 'skip on'
   write(84,  *) 'read serr  1, 2, 3'
   write(840, *) 'read serr  1, 2'
   write(841, *) 'read serr  1, 2'
   write(85,  *) 'read serr  1, 2, 3'
   write(850, *) 'read serr  1, 2'
   write(851, *) 'read serr  1, 2'
   write(86,  *) 'read serr  1, 2, 3, 4'
   write(860, *) 'read serr  1, 2, 3, 4'
   write(861, *) 'read serr  1, 2'
   write(862, *) 'read serr  1, 2'

   r_bin = r_bin * 10.
   l_bin = l_bin * 10.
   do jj = 1, freq_num 
      do k = 1, en_num 
         write(84, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, rc_fq_en(jj, k), error_Aform_rc_ic(jj, k), &
              rc_fq_en(jj, k), sqrt(var_rc(jj, k) / real(obs_freq_bins(jj) * int_number) )
         write(840, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, rc_fq_en(jj, k), error_Aform_rc_ic(jj, k)
         write(841, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, rc_fq_en(jj, k), sqrt(var_rc(jj, k) / real(obs_freq_bins(jj) * int_number) )
         
         write(85, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, ic_fq_en(jj, k), error_Aform_rc_ic(jj, k), &
              ic_fq_en(jj, k), sqrt(var_ic(jj, k) / real(obs_freq_bins(jj) * int_number) )
         write(850, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, ic_fq_en(jj, k), error_Aform_rc_ic(jj, k)
         write(851, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, ic_fq_en(jj, k), sqrt(var_ic(jj, k) / real(obs_freq_bins(jj) * int_number) )

         write(86, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, lag_freq(jj, k), error_Aform_lag(jj, k), lag_freq(jj, k), error_cohe_lag(jj, k), lag_freq(jj, k), error_prop(jj, k) / (2 * pi * relevant_freq(jj))
         write(860, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, lag_freq(jj, k), error_Aform_lag(jj, k)
         write(861, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, lag_freq(jj, k), error_cohe_lag(jj, k)
         write(862, *)  (r_bin(k) + l_bin(k)) * 0.5, (r_bin(k) - l_bin(k)) * 0.5, lag_freq(jj, k), error_prop(jj, k) / (2 * pi * relevant_freq(jj))

      enddo
      write(84,  *) 'no no no no'
      write(840, *) 'no no no no'
      write(841, *) 'no no no no'
      write(85,  *) 'no no no no'
      write(850, *) 'no no no no'
      write(851, *) 'no no no no'
      write(86,  *) 'no no no no'
      write(860, *) 'no no no no'
      write(861, *) 'no no no no'
      Write(862, *) 'no no no no'
   enddo


   write(84, *) 'scr white'
   write(84, *) 'log x on'
   ! write(84, *) 'r x 0.3 11.'
   write(84, *) 'la y Real Cross'
   write(84, *) 'la x Energy (KeV)'
   write(84, *) 'lw 5'
   write(84, *) 'tim off'
   write(85, *) 'scr white'
   write(85, *) 'log x on'
   ! write(85, *) 'r x 0.3 11.'
   write(85, *) 'la y Imaginary Cross'
   write(85, *) 'la x Energy (KeV)'
   write(85, *) 'lw 5'
   write(85, *) 'tim off'

   write(86, *) 'scr white'
   write(86, *) 'log x on'
   write(86, *) 'la y Lag (s)'
   write(86, *) 'la x Energy (KeV)'
   write(86, *) 'lw 5'
   write(86, *) 'tim off'



  write(*,*)"----------------------------------------------------------------"
  write(*,*)"Plots:"
  write(*,*)"fort.84 ..... Re part of the cross spectrum (1:Adam error; 2:prop error)"
  write(*,*)"fort.85 ..... Im part of the cross spectrum (1:Adam error; 2:prop error)"
  write(*,*)"fort.86 ..... lag spectrum (1:Adam error; 2: coherence errors; 3:prop errors)"
  write(*,*)"----------------------------------------------------------------"


!XSPEC stuff  (use flx2xsp to create the pha file)
   do jj = 1, freq_num 
      write(name_base  , '(A,I1,A)') 'rc',jj - 1,'.dat'  
      write(name_base_2, '(A,I1,A)') 'ic',jj - 1,'.dat'  
      open(11, file = trim(name_base))
      open(12, file = trim(name_base_2))
      write(11, *) '0.1 0.3 0.01 0.01 '
      write(12, *) '0.1 0.3 0.01 0.01 '
      do k = 1, en_num 
         write(11, *)  l_bin(k), r_bin(k),  (r_bin(k) - l_bin(k)) * rc_fq_en(jj, k), (r_bin(k) - l_bin(k)) * error_Aform_rc_ic(jj, k)
         write(12, *)  l_bin(k), r_bin(k),  (r_bin(k) - l_bin(k)) * ic_fq_en(jj, k), (r_bin(k) - l_bin(k)) * error_Aform_rc_ic(jj, k)
      enddo
      write(11, *) '10.0 20. 0.01 0.01 '
      write(12, *) '10.0 20. 0.01 0.01'
      close(11)
      close(12)
   enddo
   write(*,*) 'Use flx2xsp to create the pha files from rc and ic'

! !XSPEC lags 
   do jj = 1, freq_num 
      write(name_base  , '(A,I1,A)') 'lag',jj - 1,'.dat'  
      open(11, file = trim(name_base))
      write(11, *) '0.1 0.3 0.01 0.01 '
      do k = 1, en_num 
         write(11, *)  l_bin(k), r_bin(k),  (r_bin(k) - l_bin(k)) * lag_freq(jj, k), (r_bin(k) - l_bin(k)) *  error_Aform_lag(jj, k)
      enddo
      write(11, *) '10.0 20. 0.01 0.01 '
      close(11)
   enddo
   write(*,*) 'Use flx2xsp to create the pha files from the lag file'



   ! endif




end program cross


!-----------------------------------------------------------------------
      subroutine Rlog_scale(a_start,a_end,dim,out)
      implicit none
      integer :: i,dim
      real :: a_start,a_end,out(dim)
      real :: du,u1,u2,exp
      if (dim.gt.1) then
         u1=log10(a_start)
         u2=log10(a_end)      
         du=(u2-u1)/float(dim-1)
         exp=log10(a_start)
         do i=1, dim
            out(i)=10**exp
            exp=exp+du
         enddo
      else if(dim.eq.1) then
       out(1)=a_start
      else if(dim.lt.1) then
       write(*,'(A)',advance="no") 'Error! It is impossible create' 
       write(*,*) ' the log scaled array: dim is less than 1'
       stop
      endif
      return
      end
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
      subroutine Rlog_rebin_power(vector, rebin_vector, rebin_array, dim, nf, c)
!-----------------------------------------------------------------------
!     PURPOUSE:  Logarithmic rebin of a vector. Find how many bins (nf) 
!     the logarithmic rebinned vector has. c is the binning factor, 
!     dim is the dimension of the lin vector. It takes the linear x_axis 
!        
!     INPUTS:  vector -- linear vector
!              dim -- # bin of the linear vector
!              c -- rebinning factor
!      
!     OUTPUTS: rebin_vector -- logarithmic rebinned vector
!              rebin_array -- how many numbers for each rebin
!
!     ROUTINE CALLED: none
!     AUTHOR: Mastroserio 
!     DATE WRITTEN: Feb 2019
!     LAST CHANGE: 
!     NOTE: 
!           
!-----------------------------------------------------------------------
        implicit none

        integer, intent(IN)    :: dim
        integer, intent(INOUT) :: nf
        integer, intent(OUT)   :: rebin_array(dim)
        real   , intent(IN)    :: vector(dim), c 
        real   , intent(OUT)   :: rebin_vector(dim)

        integer                :: i, j, nnp, remain, cont
        integer, allocatable   :: iar(:), np(:)
 

        i = 0
        j = 0
        cont = 0 
        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i
           if (floor( c**j ).eq.remain) cont = 1
           nnp  = min( floor( c**j ) , remain )
           nnp  = max( 1 , nnp )
           i   = i + nnp
        end do
        if (cont.eq.1) nf = j
        if (cont.ne.1) nf = j - 1

        allocate(iar(0:nf))
        allocate(np(nf))

        np = 0
! First calculate binning scheme
        i = 0
        j = 0
        iar(0) = 0

        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i

           if (floor( c**j ).gt.remain) then
              np(nf)  = np(nf) + remain
              iar(nf) = dim
              i      = i + np(nf)
           else
              np(j)  = min( floor( c**j ) , remain )
              np(j)  = max( 1 , np(j) )
              i      = i + np(j)
              iar(j) = i
           endif

        end do
        

! Now do the binning
        do j = 0, nf - 1
           rebin_vector(j + 1)  = 0.0
           cont = 0
           do i = iar(j) + 1,iar(j + 1)
              cont = cont + 1
              rebin_vector(j + 1)  = rebin_vector(j + 1) + vector(i)
           end do
           rebin_array (j + 1) = cont
           rebin_vector(j + 1) = rebin_vector(j + 1) / float( np(j + 1) )
        end do

        if(allocated(iar)) deallocate(iar)
        if(allocated(np)) deallocate(np)
        
        return
      end subroutine Rlog_rebin_power
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
      subroutine Rlog_rebin(vector, rebin_vector, dim, nf, c)
!-----------------------------------------------------------------------
!     PURPOUSE:  Logarithmic rebin of a vector. Find how many bins (nf) 
!     the logarithmic rebinned vector has. c is the binning factor, 
!     dim is the dimension of the lin vector. It takes the linear x_axis 
!        
!     INPUTS:  vector -- linear vector
!              dim -- # bin of the linear vector
!              c -- rebinning factor
!      
!     OUTPUTS: rebin_vector -- logarithmic rebinned vector 
!
!     ROUTINE CALLED: none
!     AUTHOR: Mastroserio 
!     DATE WRITTEN: Feb 2019
!     LAST CHANGE: 
!     NOTE: 
!           
!-----------------------------------------------------------------------
        implicit none

        integer, intent(IN)    :: dim
        integer, intent(INOUT) :: nf
        real   , intent(IN)    :: vector(dim), c 
        real   , intent(OUT)   :: rebin_vector(dim)

        integer                :: i, j, nnp, remain, cont
        integer, allocatable   :: iar(:), np(:)

        

        i = 0
        j = 0
        cont = 0 
        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i
           if (floor( c**j ).eq.remain) cont = 1
           nnp  = min( floor( c**j ) , remain )
           nnp  = max( 1 , nnp )
           i   = i + nnp
        end do
        if (cont.eq.1) nf = j
        if (cont.ne.1) nf = j - 1

        allocate(iar(0:nf))
        allocate(np(nf))


        np = 0
! First calculate binning scheme
        i = 0
        j = 0
        iar(0) = 0

        do while( i .lt. dim)
           j      = j + 1
           remain = dim - i

           if (floor( c**j ).gt.remain) then
              np(nf)  = np(nf) + remain
              iar(nf) = dim
              i      = i + np(nf)
           else
              np(j)  = min( floor( c**j ) , remain )
              np(j)  = max( 1 , np(j) )
              i      = i + np(j)
              iar(j) = i
           endif

        end do


! Now do the binning
        do j = 0, nf - 1
           rebin_vector(j + 1)  = 0.0

           do i = iar(j) + 1, iar(j + 1)
              rebin_vector(j + 1)  = rebin_vector(j + 1) + vector(i)
           end do
           rebin_vector(j + 1)  = rebin_vector(j + 1) / float( np(j + 1) )
        end do

        if(allocated(iar)) deallocate(iar)
        if(allocated(np)) deallocate(np)
        return
      end subroutine Rlog_rebin
!-----------------------------------------------------------------------


! !---------------------------------------------------------------------!
!     subroutine split_lc()
!       use dyn_lc
!       implicit none

!       integer               :: i, new_dim_GTI
!       double precision      :: tot_time_lc
!       logical               :: yes_no, power2, check_interval
!       double precision   , allocatable  :: temp_array(:), temp_GTI1(:), temp_GTI2(:)

! ! The first time (when int_number is -1 and split_ind is not allocated) we call the lc_split_first to work out split_index array
! ! Then we call only lc_split_index

! !At this point there is the possibility to print the complete light curve with the intxoerpolated gaps            
!       ! do i=1, nrow
!       !    write(98,*) time(i), lc(i)
!       ! enddo


!       if (.not. allocated(split_ind)) then
!          write(*,*) 
!          write(*,*) 'TOTAL length of the light curve (sec): ', time(dim_lc) - time(1)
!          write(*,*) 
!          write(*,*) 'Number of gaps in the light curve: ', dim_GTI - 1
!          write(*,*) 

!          check_gap_num = dim_GTI

!          write(*,*) 'GTI invervals in seconds:'
!          do i = 1, dim_GTI
!             write(*,*) start_GTI(i), end_GTI(i)
!          enddo
!          write(*,*) 
         
         
!          if (yes_no('    Do you want to interpolate the light curve gaps?')) then 

! ! This is to save the rate before the interpolation 
!             if(.not. allocated(temp_array)) allocate(temp_array(dim_lc))
!             if(.not. allocated(temp_GTI1) ) allocate(temp_GTI1 (dim_GTI))
!             if(.not. allocated(temp_GTI2) ) allocate(temp_GTI2 (dim_GTI))
!             temp_array = lc 
!             temp_GTI1  = start_GTI
!             temp_GTI2  = end_GTI
            
! ! This is while because after the interpolation the user might want to interpolate more 
!             do             
!                ! write(*,*) 'start do loop'
! ! CALL THE INTERPOLATION SUBROUTINE TO FILL THE SMALL GAPS      
!                call interpol_split_silent(new_dim_GTI)
!                write(*,*)
!                write(*,*) 'After the interpolation the number of gaps is: ', new_dim_GTI - 1
!                write(*,*)
!                if (yes_no('   Do you want to interpolate differently?')) then
!                   lc      = temp_array
!                   start_GTI = temp_GTI1
!                   end_GTI   = temp_GTI2
!                   gap = -1
!                else
!                   exit
!                endif

!                ! write(*,*) 'end do loop'
!             enddo
            
!             if(allocated(temp_array)) deallocate(temp_array)
!             if(allocated(temp_GTI1) ) deallocate(temp_GTI1)
!             if(allocated(temp_GTI2) ) deallocate(temp_GTI2)
            
! ! RE-SET THE GTI INTERVALS BASED ON THE INTERPOLATION
!             allocate(temp_array(new_dim_GTI))
!             do i=1, new_dim_GTI
!                temp_array(i) = start_GTI(i)
!             enddo
!             deallocate(start_GTI)      
!             allocate(start_GTI(new_dim_GTI))
!             do i = 1, new_dim_GTI
!                start_GTI(i) = temp_array(i)  ! write the new stat_GTI
!                temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
!             enddo
!             deallocate(end_GTI)
!             allocate(end_GTI(new_dim_GTI))
!             do i = 1, new_dim_GTI
!                end_GTI(i) = temp_array(i)  ! write the new end_GTI
!             enddo
!             deallocate(temp_array)

!             dim_GTI = new_dim_GTI

!          else
!             check_gap_num = dim_GTI
!          endif
         
!          write(*,*)
!          tot_time_lc = 0.d0
!          do i = 1, dim_GTI
!             tot_time_lc = tot_time_lc +  end_GTI(i) - start_GTI(i)
!          enddo
!          write(*,*) '   Maximum length of the light curve (sec)', tot_time_lc
!          write(*,*) '   It might not been continuus'
!          write(*,*)

!          do 
! ! !length of the intervals
!             write(*,*) '   resolution of the light curve', dt
!             write(*,*) '   Length of the interval in steps: '
!             read(*,*) int_len_dim

! !Compute the split in the light curve and create the split_ind array
!             call lc_split_first(check_interval)
!             if (check_interval) then 
!                write(*,*) 'int_number, int_len_dim', int_number, int_len_dim
               
!                if(.not. allocated(lc_int)  ) allocate(lc_int  (int_number, int_len_dim) )
!                if(.not. allocated(time_int)) allocate(time_int(int_number, int_len_dim) )
!                ! if(.not. allocated(bkg_int) ) allocate(bkg_int (int_number, int_len_dim) )
!             ! if(.not. allocated(err_int) ) allocate(err_int (int_number, int_len_dim) )

!             !Fill the time_int, lc_int and bkg_int which are the light curves separated in intervals
!                call lc_split_index_time()
! ! Check if int_len_dim is a power of 2
!                if (power2(int_len_dim) .eqv. .false.)then
!                   write(*,*) '   Precedure without the FFT (it is not a power of 2)'
!                   check_power2 = .false.
!                else
!                   check_power2 = .true.
!                endif

!                exit
!             endif
!             write(*,*) '  '
!             write(*,*) '  Try again...'
!          enddo

!       else
! !THIS is called if it's not the first time

! !Call interpol_split_silent because the new lc needs interpolation. The number of gaps is saved from the first time          
! !There is no need to adjust the GTI arrays because the split is already done. 
! ! The IF statement before calling the subroutine is because if it is false the first time it should be false every time (no gap interpolation!!)
!          ! write(*,*) 'gaaaaaaaaaaaaaaaaaap', gap

!          if (gap .ne. -1) then 
!             call interpol_split_silent(new_dim_GTI)
!          endif
!             allocate(temp_array(new_dim_GTI))
!             do i=1, new_dim_GTI
!                temp_array(i) = start_GTI(i)
!             enddo
!             deallocate(start_GTI)      
!             allocate(start_GTI(new_dim_GTI))
!             do i = 1, new_dim_GTI
!                start_GTI(i) = temp_array(i)  ! write the new stat_GTI
!                temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
!             enddo
!             deallocate(end_GTI)
!             allocate(end_GTI(new_dim_GTI))
!             do i = 1, new_dim_GTI
!                end_GTI(i) = temp_array(i)  ! write the new end_GTI
!             enddo
!             deallocate(temp_array)

!             dim_GTI = new_dim_GTI
         
!          print *, ' Number of gaps in the total light curve: : ', dim_GTI - 1

!          call lc_split_index_time()
!       endif


     
!     end subroutine split_lc
! !---------------------------------------------------------------------!
!---------------------------------------------------------------------!
    subroutine split_lc()
      use dyn_lc
      implicit none

      integer               :: i, new_dim_GTI
      real                  :: tot_time_lc
      logical               :: yes_no, power2, check_interval, check 
      real   , allocatable  :: temp_array(:), temp_GTI1(:), temp_GTI2(:)

! The first time (when int_number is -1 and split_ind is not allocated) we call the lc_split_first to work out split_index array
! Then we call only lc_split_index

!At this point there is the possibility to print the complete light curve without the interpolated gaps            
      ! do i=1, nrow
      !    write(98,*) time(i), lc(i)
      ! enddo


if (.not. allocated(split_ind)) then
   write(*,*) 
   write(*,*) 'TOTAL length of the light curve (sec) considering the gaps (final_time - starting time): ', time(dim_lc) - time(1)
   write(*,*) 
   write(*,*) 'Number of gaps in the light curve: ', dim_GTI - 1
   write(*,*) 

   check_gap_num = dim_GTI


! Check if you can interpolate with the current routine: the light curve needs no jumps in time 
   check = .true.
   do i = 2 , dim_lc
      !         write(*,*) 'time index', i
      ! write(*,*) time(i) - time(i - 1)
      if ((time(i) - time(i - 1)) .gt. 10 * dt  ) then
         check = .false.
      endif
   enddo

   if (check) then
      if (yes_no('    Do you want to interpolate the light curve gaps?')) then 

! This is to save the rate before the interpolation 
         if(.not. allocated(temp_array)) allocate(temp_array(dim_lc))
         if(.not. allocated(temp_GTI1) ) allocate(temp_GTI1 (dim_GTI))
         if(.not. allocated(temp_GTI2) ) allocate(temp_GTI2 (dim_GTI))
         temp_array = lc 
         temp_GTI1  = start_GTI
         temp_GTI2  = end_GTI

! This is a while loop because after the interpolation the user might want to interpolate differently
         do             
               ! write(*,*) 'start do loop'
! CALL THE INTERPOLATION SUBROUTINE TO FILL THE SMALL GAPS      
            call interpol_split_silent(new_dim_GTI)
            write(*,*)
            write(*,*) 'After the interpolation the number of gaps is: ', new_dim_GTI - 1
            write(*,*)
            if (yes_no('   Do you want to interpolate differently?')) then
               lc      = temp_array
               start_GTI = temp_GTI1
               end_GTI   = temp_GTI2
               gap = -1
            else
               exit
            endif

            ! write(*,*) 'end do loop'
         enddo

         if(allocated(temp_array)) deallocate(temp_array)
         if(allocated(temp_GTI1) ) deallocate(temp_GTI1)
         if(allocated(temp_GTI2) ) deallocate(temp_GTI2)
            
! RE-SET THE GTI INTERVALS BASED ON THE INTERPOLATION
         allocate(temp_array(new_dim_GTI))
         do i=1, new_dim_GTI
            temp_array(i) = start_GTI(i)
         enddo
         deallocate(start_GTI)      
         allocate(start_GTI(new_dim_GTI))
         do i = 1, new_dim_GTI
            start_GTI(i) = temp_array(i)  ! write the new stat_GTI
            temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
         enddo
         deallocate(end_GTI)
         allocate(end_GTI(new_dim_GTI))
         do i = 1, new_dim_GTI
            end_GTI(i) = temp_array(i)  ! write the new end_GTI
         enddo
         deallocate(temp_array)

         dim_GTI = new_dim_GTI

      else
         check_gap_num = dim_GTI
      endif

   else
      write(*,*) 'No interpolation of the light curve is possible with the current routine'
   endif
         
   write(*,*)
   tot_time_lc = 0.0
   do i = 1, dim_GTI
      tot_time_lc = tot_time_lc +  end_GTI(i) - start_GTI(i)
   enddo
   write(*,*) '   Length of the light curve, summing all the pieces', tot_time_lc, '(sec)'
   write(*,*)


   do 
! !length of the intervals
      write(*,*) '   Length of the interval in steps: '
      read(*,*) int_len_dim

!Compute the split in the light curve and create the split_ind array
      call lc_split_first(check_interval)
      if (check_interval) then 

         if(.not. allocated(lc_int)  ) allocate(lc_int  (int_number, int_len_dim) )
         if(.not. allocated(time_int)) allocate(time_int(int_number, int_len_dim) )
         if(.not. allocated(bkg_int) ) allocate(bkg_int (int_number, int_len_dim) )
         ! if(.not. allocated(err_int) ) allocate(err_int (int_number, int_len_dim) )

               
         !Fill the time_int, lc_int and bkg_int which are the light curves separated in intervals
         call lc_split_index_time()

! Check if int_len_dim is a power of 2
         if (power2(int_len_dim) .eqv. .false.)then
            write(*,*) '   Precedure without the FFT (it is not a power of 2)'
            check_power2 = .false.
         else
            check_power2 = .true.
         endif
         
         exit
      endif
      write(*,*) '  '
      write(*,*) '  Try again...'
   enddo

else !THIS is called if it's not the first time

!Call interpol_split_silent because the new lc needs interpolation. The number of gaps is saved from the first time          
!There is no need to adjust the GTI arrays because the split is already done. 
! The IF statement before calling the subroutine is because if it is false the first time it should be false every time (no gap interpolation!!)
   ! write(*,*) 'gaaaaaaaaaaaaaaaaaap', gap, new_dim_GTI
         
   if (gap .ne. -1) then 
      call interpol_split_silent(new_dim_GTI)

      allocate(temp_array(new_dim_GTI))
      do i=1, new_dim_GTI
         temp_array(i) = start_GTI(i)
      enddo
      deallocate(start_GTI)      
      allocate(start_GTI(new_dim_GTI))
      do i = 1, new_dim_GTI
         start_GTI(i) = temp_array(i)  ! write the new stat_GTI
         temp_array(i) = end_GTI(i)    ! re-write the temporary array with end_GTI
      enddo
      deallocate(end_GTI)
      allocate(end_GTI(new_dim_GTI))
      do i = 1, new_dim_GTI
         end_GTI(i) = temp_array(i)  ! write the new end_GTI
      enddo
      deallocate(temp_array)
   endif
         


   print *, ' Number of gaps in the total light curve: : ', dim_GTI - 1
   
   call lc_split_index_time()
endif


     
    end subroutine split_lc
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
  subroutine interpol_split_silent(new_dim_GTI)
    !This subroutine differs from interpol_split only because when it asks to  
    !set the maximum number of bins that you want to interpolate, it has the 
    !option to skip this if it has been already asked. This is why is called 
    !SILENT, it has the option to be silent. 

    !This subroutine fills the gaps in the light curve according to the maximum number
    ! of bins that the user chooses.
    !Then the GTI arrays (start and end) are modified and a new dimension is set.
    !Basically the become shorter.
    !NOTE: after this function the GTI arrays are re-set (deallocate and reallocate)
    !      with the correct number of bins
    
    use dyn_lc
    use rand_fun
    implicit none
    integer, intent(OUT)   :: new_dim_GTI

    integer    :: j, i, k, w, ee,  max_gap
    real       :: m, q, lc_interpol, m_b, q_b, bkg_interpol
    double precision, allocatable   :: new_start_GTI(:), new_end_GTI(:), max_gap_sec 

    ! do j = 1, dim_GTI
    !    write(*,*) start_GTI(j), end_GTI(j)
    ! enddo
    

    allocate(new_start_GTI(dim_GTI))
    allocate(new_end_GTI(dim_GTI))
    
    w = 1 
    new_start_GTI(w) = start_GTI(w)

    if (gap .eq. -1) then 

       write(*,*)
       write(*,*) '  The gap situation is: '
       do j = 1, dim_GTI - 1
          ee = int( ( start_GTI(j + 1) - end_GTI(j) ) / dt)
          write(*,*) 'Gap number: ', j, real(ee) * real(dt), 'sec', ee, 'bins'
       enddo

       write(*,*)
       write(*,*) '   Set the maximum length that you want to '
       write(*,*) '   interpolate. Remember that dt is ', real(dt)
       read(*,*) max_gap_sec
       ! read(*,*) max_gap
       gap = int(max_gap_sec / dt)
       max_gap = gap
    else
       max_gap = gap
    endif 
    
    do j = 1, dim_GTI - 1
       ee = int( ( start_GTI(j + 1) - end_GTI(j) ) / dt)


       if (ee .ge. 1 .and. ee .le. max_gap ) then
          i = 1
          do while(time(i) .lt. end_GTI(j))
             i = i + 1
             if (i .gt. dim_lc) then
                exit
             endif
          enddo

          ! write(*,*)  i, ee + i,  time(ee+i), lc(ee+i), time(i-1), lc(i-1)

          m = (lc(ee + i) - lc(i - 1)) / real(time(ee + i) - time(i - 1))
          q = lc(i - 1) - m * time(i - 1)
          if (allocated(bkg)) then
             write(*,*) 'bkg active'
             m_b = (bkg(ee + i) - bkg(i - 1)) / real(time(ee + i) - time(i - 1))
             q_b = bkg(i - 1) - m_b * real(time(i - 1))
          endif
          
          do k = i, i + ee - 1
             lc_interpol = m * real(time(k)) + q
             lc_interpol = lc_interpol * real(dt)
             ! write(*,*) 'time interpol', time(k), lc_interpol
             
             if (allocated(bkg)) then 
                bkg_interpol  = m_b * real(time(k)) + q_b
                bkg_interpol = bkg_interpol * real(dt)
             endif
!Extracting the interpolation value from a Poisson distribution with the same mean
             lc(k) = real(poidev(dble(lc_interpol)) / dt)
             ! write(*,*) 'extract!', lc(k)
             ! if (lc(k) .lt. 0.0) lc(k) = 0.0

             ! lc(k) = m * time(k) + q
             if (allocated(bkg)) then  
                bkg(k)  = real(poidev(dble(bkg_interpol)) / dt)
             endif 
             ! write(*,*) 'bkg', bkg(k)
             ! if (bkg(k) .lt. 0.0) bkg(k) = 0.0
             
          enddo
          new_end_GTI(w) = end_GTI(j+1)
       else
          new_end_GTI(w) = end_GTI(j)
          w = w + 1
          new_start_GTI(w) = start_GTI(j+1)
       endif
       
    enddo
    new_end_GTI(w) = end_GTI(j)

    ! write(*,*) "W" ,w 
    do j = 1, dim_GTI
       start_GTI(j) = new_start_GTI(j)
       end_GTI(j) = new_end_GTI(j)
       ! write(*,*)   new_start_GTI(j),new_end_GTI(j)
    enddo
    
    new_dim_GTI = w

    
  end subroutine interpol_split_silent
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
    subroutine lc_split_first(check_interval)
!! Working out the intervals in the light curve
! The best way to do that is to check with the GTI
      use dyn_lc
      implicit none
      logical, intent(OUT) :: check_interval
      
      integer :: i, j, count, count2
      logical :: check


      if(.not. allocated(split_ind)) allocate(split_ind(dim_lc))
      split_ind = -1

      
! Counter: this first part work out how many intervals there are in the complete light curve
!      the calculation is based on the time difference
      count  = 0   !count how many elements to form an interval  
      count2 = 1  !count how many intervals
      split_ind(1) = 1
      
!First step
      write(*,*) 'First time of the light curve ', time(1)
      write(*,*) 'First time of the first GTI', start_GTI(1)
      if( (time(1) .ge. start_GTI(1))  .and. (time(1) .le. end_GTI(1)) ) then
         count = count + 1
      endif

      
      do i = 2 , dim_lc
!         write(*,*) 'time index', i
         check = .false.
         if ((time(i) - time(i - 1)) .gt. dt * 10 ) then
            write(*,*) 'jump in time of the light curve'
            write(*,*) 'from ', time(i - 1), 'to ', time(i)
            count = 0
            
         else 
            
            do j = 1, dim_GTI
            ! write(*,*) time(i), start_GTI(j), end_GTI(j) 
               if( (time(i) .ge. start_GTI(j))  .and. (time(i) .le. end_GTI(j)) ) then
                  check = .true.
                  ! write(*,*) 'exit', check
                  exit
               endif
            enddo
         endif
         
        if (check) then
           count = count + 1
           ! write(11,*) check, count
        else 
           count = 0
           ! write(22,*) check, time(i), count
           split_ind((count2 * 2) - 1) = i + 1
        end if
        
        if (count .eq. int_len_dim) then
           split_ind( count2 * 2 ) = i  
           count2 = count2 + 1
           count = 0
           split_ind((count2 * 2) - 1) = i + 1
        end if        
     enddo

     int_number = count2 - 1 ! There is always one more than what we need 
! Check if there are at least one interval in the light curve (it can happen that the light curve has too many holes or the int_len_dim is too big)
     check_interval = .true.
     if (int_number .lt. 1 ) then
        write(*,*) '   The interval chosen is too long for the light curve.'
        if (allocated(split_ind)) deallocate(split_ind)
        check_interval = .false.
     endif

     
   end subroutine lc_split_first
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
    subroutine lc_split_index_time()
! This suroutine writes in lc matrix all the intervals based on the split_ind array      
      use dyn_lc
      implicit none
      integer   :: i, j, count

      do i = 1, int_number
         ! write(*,*) i, split_ind(2 * i - 1), split_ind(2 * i)
         count = 1 
         do j = split_ind(2 * i - 1), split_ind(2 * i)
            lc_int  (i, count) = lc(j)
            time_int(i, count) = time(j)
            ! if(allocated(bkg)) bkg_int (i, count) = bkg(j)
           ! if (lc(j) .lt. 0.0) lc(j) = 0.0
           ! if (bkg(j) .lt. 0.0) bkg(j) = 0.0

           ! if (lc(j) .lt. 0.0) write(*,*) 'ATTENTION!! (in lc_split_index_time) negative lc (time, lc)', time(j), lc(j)
           ! if (bkg(j) .lt. 0.0) write(*,*) 'ATTENTION!! (in lc_split_index_time) negative bkg (time, bkg)', time(j), bkg(j)
           count = count + 1
        enddo
     enddo

   end subroutine lc_split_index_time
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
    subroutine extract_lc(filename)
      use dyn_lc
      implicit none

! CFITSIO variable
      character (len=200) :: filename

      character (len=30)  :: error_description, hdu_name &
                             ,col_name_temp, keyword
      logical             :: anynul
      integer             :: i, status, readwrite, blocksize, unit, chdu  &
                             ,nrow, nrow_GTI, ncol, ncol_GTI, colnum & 
                             ,felem, nelem, datacode, repeat, width, frow
      integer             :: colnum_f
      real                :: dt_f, dt_temp
      real                :: nullval
      double precision    :: first_time_bin

      real  , dimension(:), allocatable :: time_e, rate_e, bkg_e, start_GTI_e, end_GTI_e, err_rate_e
      double precision, dimension(:), allocatable :: time_d, rate_d, bkg_d, start_GTI_d, end_GTI_d, err_rate_d


! Cross-spectum variables
!      integer :: i,count,count2,dim
      

      status = 0
      readwrite = 0  !file access mode: readonly = 0 , read and write = 1
! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)
! Open the FITS file
      call ftopen(unit,filename,readwrite,blocksize,status)
!      call ftnopn(unit,filename,readwrite,status)  !thisi is used with the extension in the name (e.g. namefile.fits+2)
      
      call ftgerr(status,error_description)
      if (status .ne. 0) then
         write(*,*) '!! ATTENTION !! ', error_description
         write(*,*) '    Exit...'
         stop
      endif
      
!Let's move in the correct HDU (with the name in hdu_name)      
      hdu_name = 'RATE'
      call hdu_move(unit,hdu_name,status)
      call ftghdn(unit,chdu)
      status = 0 
!      write(*,*) '   Current hdu',chdu

! Get the number of rows and columns in the CHDU
      call ftgnrw(unit,nrow,status)
      call ftgncl(unit,ncol,status)
!      write(*,*) "number of rows", nrow


!**************************************************************!      
! Let's get the column number for a particolar name with colnum_f (function)
      col_name_temp = '*TIME*'
      write(*,*) 'Look at the *TIME* column'
      colnum = colnum_f(unit,col_name_temp,status)
      write(*,*) 'TIME found'

! Get the datatype of a column and read that column
!   repeat tells you if there is more than one element in very space of the table  

! Get the datatype
      call ftgtcl(unit,colnum,datacode,repeat,width,status)
       if( status .gt. 0 )call printerror(status)

! Get the the column values      
      frow = 1 !starting row
      nelem = nrow ! last row to read
      felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         write(*,*) "time is a real"
         if(.not.allocated(time_e)) allocate(time_e(nrow))
         call ftgcve(unit,colnum,frow,felem,nelem,nullval,time_e,anynul,status)         
      else if (datacode .eq. 82) then
         write(*,*) "time is a double"
         if(.not.allocated(time_d)) allocate(time_d(nrow))
         call ftgcvd(unit,colnum,frow,felem,nelem,nullval,time_d,anynul,status)
      else
         write(*,*) "   TIME column is niether a real nor a double"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      
     
!**************************************************************!      
      col_name_temp = '*RATE*'
      write(*,*) 'Look at the *RATE* column'
      colnum = colnum_f(unit,col_name_temp,status)
      write(*,*) 'RATE found'
      call ftgtcl(unit,colnum,datacode,repeat,width,status)
       if( status .gt. 0 )call printerror(status)

      frow = 1 !starting row
      nelem = nrow ! last row to read
      felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         write(*,*) 'rate real'
         if(.not.allocated(rate_e)) allocate(rate_e(nrow))
         call ftgcve(unit,colnum,frow,felem,nelem,nullval,rate_e,anynul,status)         
      else if (datacode .eq. 82) then
         write(*,*) 'rate double'
         if(.not.allocated(rate_d)) allocate(rate_d(nrow))
         call ftgcvd(unit,colnum,frow,felem,nelem,nullval,rate_d,anynul,status)         
      else
         write(*,*) "   RATE column is niether an real nor a double"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      

! !**************************************************************!      
      col_name_temp = '*ERROR*'
      write(*,*) 'Look at the *ERROR* column'
      colnum = colnum_f(unit, col_name_temp, status)
      if (colnum .gt. 0.0) then 
         write(*,*) 'ERROR found'
         call ftgtcl(unit, colnum, datacode, repeat, width, status)
         if( status .gt. 0 )call printerror(status)

         frow    = 1    !starting row
         nelem   = nrow ! last row to read
         felem   = 1    ! first pixel of the element vector (ignored for ASCII tables)
         nullval = -1
         if(datacode .eq. 42) then
            if(.not.allocated(err_rate_e)) allocate(err_rate_e(nrow))
            call ftgcve(unit, colnum, frow, felem, nelem, nullval, err_rate_e, anynul, status)
         else if (datacode .eq. 82) then
            if(.not.allocated(err_rate_d)) allocate(err_rate_d(nrow))
            call ftgcvd(unit, colnum, frow, felem, nelem, nullval, err_rate_d, anynul, status) 
         else
            write(*,*) "   ERROR column is niether an real nor a double"
            write(*,*) "   Modify the source code "
         end if
      endif
! !**************************************************************!      


!**************************************************************!      
      ! col_name_temp = '*BACKV*'
      ! write(*,*) 'Look at the *BACKV* column'
      ! colnum        = colnum_f(unit, col_name_temp, status)
      ! if (colnum .gt. 0.0) then 
      !    write(*,*) 'BACKV found'
      !    call ftgtcl(unit, colnum, datacode, repeat, width, status)
      !    if( status .gt. 0 )call printerror(status)

      !    frow    = 1    !starting row
      !    nelem   = nrow ! last row to read
      !    felem   = 1    ! first pixel of the element vector (ignored for ASCII tables)
      !    nullval = -1
      !    if(datacode .eq. 42) then
      !       if(.not.allocated(bkg_e)) allocate(bkg_e(nrow))
      !       call ftgcve(unit, colnum, frow, felem, nelem, nullval, bkg_e, anynul, status)
      !    else if (datacode .eq. 82) then
      !       if(.not.allocated(bkg_d)) allocate(bkg_d(nrow))
      !       call ftgcvd(unit, colnum, frow, felem, nelem, nullval, bkg_d, anynul, status) 
      !    else
      !       write(*,*) "   BACKV column is niether an real nor a double"
      !       write(*,*) "   Modify the source code "
      !    end if
      ! endif
!**************************************************************!      


            
! extract the dt from fits file with dt_f (function)
      keyword = 'TIMEDEL'
      dt_temp = dt
      dt = dt_f(unit,keyword,status)

! Check if the dt is equal to the previous one, besides the first time the dt is calculated (starting dt=-1)
      if(dt_temp .ne. -1) then 
         if (dt .ne. dt_temp) then
            write(*,*) '   dt of this new light curve is different from the previous one'
            write(*,*) '   to make the cross spectrum you need to have the same dt. EXIT...'
            stop
         endif
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      ! ! GET THE GTI ARRAY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!Let's move in the correct HDU (with the name in hdu_name)      
      hdu_name = 'GTI'
      call hdu_move(unit, hdu_name,status)
      call ftghdn(unit, chdu)
      status = 0 
! Get the number of rows and columns in the CHDU
! and check if the number of columns is equal to 2      
      call ftgnrw(unit, nrow_GTI, status)
      call ftgncl(unit, ncol_GTI, status)
      if (ncol_GTI .ne. 2) then
         write(*,*) '    The selected GTI extension has not 2 columns! Exit...'
         stop
      endif
      
!**************************************************************!      
      col_name_temp = '*START*'
      colnum = colnum_f(unit, col_name_temp, status)
      call ftgtcl(unit, colnum, datacode, repeat, width, status)
       if( status .gt. 0 )call printerror(status)

      frow    = 1 !starting row
      nelem   = nrow_GTI ! last row to read
      felem   = 1 ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         if(.not.allocated(start_GTI_e)) allocate(start_GTI_e(nrow_GTI))
         call ftgcve(unit, colnum, frow, felem, nelem, nullval, start_GTI_e, anynul, status)         
      else if (datacode .eq. 82) then
         if(.not.allocated(start_GTI_d)) allocate(start_GTI_d(nrow_GTI))
         call ftgcvd(unit, colnum, frow, felem, nelem, nullval, start_GTI_d, anynul, status)         
      else
         write(*,*) "   START GTI column is niether an real nor a double"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      

!**************************************************************!      
      col_name_temp = '*STOP*'
      colnum = colnum_f(unit, col_name_temp, status)
      call ftgtcl(unit, colnum, datacode, repeat, width, status)
      if( status .gt. 0 )call printerror(status)

      frow    = 1        !starting row
      nelem   = nrow_GTI ! last row to read
      felem   = 1        ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         if(.not.allocated(end_GTI_e)) allocate(end_GTI_e(nrow_GTI))
         call ftgcve(unit, colnum, frow, felem, nelem, nullval, end_GTI_e, anynul, status)         
      else if (datacode .eq. 82) then
         if(.not.allocated(end_GTI_d)) allocate(end_GTI_d(nrow_GTI))
         call ftgcvd(unit, colnum, frow, felem, nelem, nullval, end_GTI_d, anynul, status)         
      else
         write(*,*) "   END GTI column is niether an real nor a double"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      

! !CONVERSION FROM DOUBLE PRECISION TO REAL OF TIME_D, START
!       first_time_bin = time_d(1)
!       if (allocated(time_d)) then
!          if(.not. allocated(time_e)) allocate(time_e(nrow))
!          do i = 1, nrow
!             time_e(i) = real(time_d(i) - first_time_bin)
!          enddo
!          deallocate(time_d)
!       endif

! !CONVERSION FROM DOUBLE PRECISION TO REAL OF RATE_D, START
      if (allocated(rate_d)) then
         if(.not. allocated(rate_e)) allocate(rate_e(nrow))
         do i = 1, nrow
            rate_e(i) = real(rate_d(i))
         enddo
         deallocate(rate_d)
      endif

! !CONVERSION FROM DOUBLE PRECISION TO REAL OF BKG_D, START
!       if (allocated(bkg_d)) then
!          if(.not. allocated(bkg_e)) allocate(bkg_e(nrow))
!          do i = 1, nrow
!             bkg_e(i) = real(bkg_d(i))
!          enddo
!          deallocate(bkg_d)
!       endif
      
!       if (allocated(start_GTI_d)) then
!          if(.not. allocated(start_GTI_e)) allocate(start_GTI_e(nrow_GTI))
!          do i = 1, nrow_GTI
!             start_GTI_e(i) = real(start_GTI_d(i) - first_time_bin)
!          enddo
!          deallocate(start_GTI_d)
!       endif

!        if (allocated(end_GTI_d)) then
!          if(.not. allocated(end_GTI_e)) allocate(end_GTI_e(nrow_GTI))
!          do i = 1, nrow_GTI
!             end_GTI_e(i) = real(end_GTI_d(i) - first_time_bin)
!          enddo
!          deallocate(end_GTI_d)
!       endif

! HERE it is possible to print the GTI and the gaps
!       write(*,*) '-----------------------------------------'
!       write(*,*) 'GTI: ', nrow_GTI
!       do i = 1, nrow_GTI - 1
!          write(*,*) start_GTI_e(i), end_GTI_e(i)
!          write(*,*) 'ok time ',  end_GTI_e(i) -  start_GTI_e(i)
!          write(*,*) 'gap ', start_GTI_e(i + 1) - end_GTI_e(i)
!       enddo
!       write(*,*) start_GTI_e(nrow_GTI), end_GTI_e(nrow_GTI)
!       write(*,*) '-----------------------------------------'

! HERE it is possible to print the light curve before we fill the gaps

!       do i=1, nrow
!          write(99,*) time_e(i),rate_e(i)
!       enddo

!       do i=1, nrow_GTI
!          write(80,*) start_GTI_e(i),end_GTI_e(i)
!       enddo
      
!Allocation of the general array
      if(.not. allocated(lc) ) then
         write(*,*) 'allocation lc', nrow
         allocate(lc (nrow))
      else 
         if (nrow .ne. dim_lc) then 
            write(*,*) '   ATTENTION!! The light curves do not have the same length!'
            stop 
         endif
      endif
      
      if(.not. allocated(time) ) allocate(time (nrow))
      if(.not. allocated(err_rate)) allocate(err_rate(nrow))
      ! if(.not. allocated(bkg)  ) allocate(bkg  (nrow))

!Set the dimension of the lc and write the lc, time, and bkg in the common arrays
      dim_lc = nrow
      
      ! do i = 1, nrow
      !    lc(i)   = rate_e(i) 
      !    time(i) = time_e(i)
      !    bkg(i)  = bkg_e(i)
      ! enddo

      if (allocated(rate_e)) then 
         do i = 1, nrow
            lc(i)   = rate_e(i)
         enddo
      else
         deallocate(lc)
         write(*,*) '!! ATTENTION !! No RATE column'
      endif
      ! if (allocated(rate_e)) then 
      !    do i = 1, nrow
      !       lc(i)   = rate_e(i)
      !    enddo
      ! else
      !    deallocate(lc)
      !    write(*,*) '!! ATTENTION !! No RATE column'
      ! endif

      if (allocated(time_d)) then 
         do i = 1, nrow
            time(i) = time_d(i)
         enddo
      else 
         deallocate(time)
         write(*,*) '!! ATTENTION !! No TIME column'
      endif
      ! if (allocated(time_e)) then 
      !    do i = 1, nrow
      !       time(i) = time_e(i)
      !    enddo
      ! else 
      !    deallocate(time)
      !    write(*,*) '!! ATTENTION !! No TIME column'
      ! endif

      ! if (allocated(err_rate_e)) then 
      !    do i = 1, nrow
      !       err_rate(i)  = err_rate_e(i)
      !    enddo
      ! else 
      !    deallocate(err_rate)
      !    write(*,*) '!! ATTENTION !! No ERROR column'
      ! endif

      ! if (allocated(bkg_e)) then 
      !    do i = 1, nrow
      !       bkg(i)  = bkg_e(i)
      !    enddo
      ! else 
      !    deallocate(bkg)
      !    write(*,*) '!! ATTENTION !! No BKG column'
      ! endif


      
!Fill the GTI 
!It is complicated because we have to distinguish between the first call and the others
! and check if the GTIs are the same in all the calls (checking if the GTIs are identical to the previous call)      

      first_time_bin = start_GTI_d(1)
      ! write(*,*) 'first_time_bin', first_time_bin

      if(.not. allocated(start_GTI)) then  
         write(*,*) 'First call for filling the GTI'
         allocate(start_GTI(nrow_GTI))
         allocate(end_GTI(nrow_GTI))
!Set the dimension of the GTI and write them in the common arrays      
         dim_GTI = nrow_GTI
!Fill the actual GTI         
         do i = 1, dim_GTI
            ! start_GTI(i) = start_GTI_e(i)
            ! end_GTI  (i) = end_GTI_e(i)
            start_GTI(i) = start_GTI_d(i) - first_time_bin
            end_GTI  (i) = end_GTI_d(i)   - first_time_bin
            write(*,*) start_GTI(i), end_GTI  (i) 
         enddo
      else
         write(*,*) 'This is not the first time the GTI filling is called'
         if (check_gap_num .ne. nrow_GTI) then 
            write(*,*) '  ATTENTION!! The GTIs are not the same in all the light curves'
            stop 
         else 
!deallocate the GTIs because they have been modified by the interpolation
            write(*,*) 'GTI have the same number of rows of the first light curve'
            deallocate(start_GTI)
            deallocate(end_GTI)
!allocate the new GTIs (they will be modify by interpolation)
            allocate(start_GTI(nrow_GTI))
            allocate(end_GTI(nrow_GTI))
            dim_GTI = nrow_GTI
         do i = 1, dim_GTI
            ! start_GTI(i) = start_GTI_e(i)
            ! end_GTI  (i) = end_GTI_e(i)
            start_GTI(i) = start_GTI_d(i) - first_time_bin
            end_GTI  (i) = end_GTI_d(i)   - first_time_bin
         enddo
         endif
      endif

! Close and show if there are errors      
      call ftclos(unit,status)
      call ftgerr(status,error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description

      if(allocated(time_e)) deallocate(time_e)
      if(allocated(time_d)) deallocate(time_d)
      if(allocated(rate_e)) deallocate(rate_e)
      if(allocated(rate_d)) deallocate(rate_d)
      if(allocated(bkg_e) ) deallocate(bkg_e )
      if(allocated(err_rate_e) ) deallocate(err_rate_e )
      
      if(allocated(start_GTI_e)) deallocate(start_GTI_e)
      if(allocated(end_GTI_e)) deallocate(end_GTI_e)
      if(allocated(start_GTI_d)) deallocate(start_GTI_d)
      if(allocated(end_GTI_d)) deallocate(end_GTI_d)

    end subroutine extract_lc
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
  function channel_energy_conversion(filename, channel)
    implicit none 
    integer, intent(IN)    :: channel 
    character (len = 200)  :: filename 

    character (len = 30)   :: error_description, hdu_name, col_name_temp
    integer                :: status, readwrite, blocksize, chdu, colnum, unit, &
                              datacode, repeat, width, frow, nelem, felem
    logical                :: anynul
    real                   :: nullval
    integer                :: colnum_f
    integer                :: nrow, ncol, index_energy
    real   , allocatable   :: energy1(:), energy2(:), channels(:)
    
    real                   :: channel_energy_conversion

      status = 0
      readwrite = 0  !file access mode: readonly = 0 , read and write = 1
! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit, status)
! Open the FITS file
      call ftopen(unit, filename, readwrite, blocksize, status)
      
      call ftgerr(status, error_description)
      if (status .ne. 0) then
         write(*,*) '!! ATTENTION !! ', error_description
         write(*,*) '    Exit...'
         stop
      endif
      
!Let's move in the correct HDU (with the name in hdu_name)      
      hdu_name = 'EBOUNDS'
      call hdu_move(unit, hdu_name, status)
      call ftghdn(unit, chdu)
      status = 0 
!      write(*,*) '   Current hdu',chdu

! Get the number of rows and columns in the CHDU
      call ftgnrw(unit, nrow, status)
      call ftgncl(unit, ncol, status)

      allocate(energy1 (nrow))
      allocate(energy2 (nrow))
      allocate(channels(nrow))

!**************************************************************!      
! Let's get the column number for a particolar name with colnum_f (function)
      col_name_temp = '*E_MIN*'
      colnum = colnum_f(unit, col_name_temp, status)

!   repeat tells you if there is more than one element in very space of the table  

! Get the datatype
      call ftgtcl(unit, colnum, datacode, repeat, width, status)
       if( status .gt. 0 )call printerror(status)

! Get the the column values      
      frow = 1 !starting row
      nelem = nrow ! last row to read
      felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         call ftgcve(unit, colnum, frow, felem, nelem, nullval, energy1, anynul, status)         
      else
         write(*,*) "   ENERGY columns in the response matrix is not a real"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      

!**************************************************************!      
! Let's get the column number for a particolar name with colnum_f (function)
      col_name_temp = '*E_MAX*'
      colnum = colnum_f(unit, col_name_temp, status)

!   repeat tells you if there is more than one element in very space of the table  

! Get the datatype
      call ftgtcl(unit, colnum, datacode, repeat, width, status)
       if( status .gt. 0 )call printerror(status)

! Get the the column values      
      frow = 1 !starting row
      nelem = nrow ! last row to read
      felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
      nullval = -1
      if(datacode .eq. 42) then
         call ftgcve(unit, colnum, frow, felem, nelem, nullval, energy2, anynul, status)         
      else
         write(*,*) "   ENERGY columns in the response matrix is not a real"
         write(*,*) "   Modify the source code"
      end if
!**************************************************************!      



! !**************************************************************!      
! ! Let's get the column number for a particolar name with colnum_f (function)
!       col_name_temp = '*CHANNEL*'
!       colnum = colnum_f(unit, col_name_temp, status)

! !   repeat tells you if there is more than one element in very space of the table  

! ! Get the datatype
!       call ftgtcl(unit, colnum, datacode, repeat, width, status)
!        if( status .gt. 0 )call printerror(status)

! ! Get the the column values      
!       frow = 1 !starting row
!       nelem = nrow ! last row to read
!       felem = 1 ! first pixel of the element vector (ignored for ASCII tables)
!       nullval = -1
!       if(datacode .eq. 42) then
!          call ftgcve(unit, colnum, frow, felem, nelem, nullval, channels, anynul, status)         
!       else
!          write(*,*) "   ENERGY columns in the response matrix is not a real"
!          write(*,*) "   Modify the source code"
!       end if
! !**************************************************************!      
 
      call ftclos(unit, status)
      call ftgerr(status, error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
      
      index_energy = channel + 1 

      channel_energy_conversion = (energy2(index_energy) + energy1(index_energy)) * 0.5   

      deallocate(energy1 )
      deallocate(energy2 )
      deallocate(channels)

  end function channel_energy_conversion  
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
    subroutine hdu_move(unit,hdu_name,status)
!PURPOSE:  Move to the HDU based on the name. 
      implicit none
      integer :: unit,hdutype,hdu_extver,status
      character (len=30) hdu_name,error_description

! The hdutype parameter may have a value of IMAGE HDU(0), ASCII TBL (1),
! BINARY TBL (2), or ANY HDU (-1) where ANY HDU means that
! only the extname and extver values will be used to locate the correct extension.
! If the input value of extver is 0 then the EXTVER keyword is ignored and the first HDU with a matching
! EXTNAME (or HDUNAME) keyword will be found. If no matching HDU is found in the file
! then the current HDU will remain unchanged and a status = BAD HDU NUM (301) will be returned.     
      hdutype = -1
      hdu_extver = 0
      call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      call ftgerr(status,error_description)
      if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
      
! if the extension is not RATE the code asks for a new name until it gets a right one      
      do while(status .eq. 301)
         write(*,*) "   No extension called: ",hdu_name
         write(*,*) "   Please specify the correct name or type 'no' to quit"
         read(*,*) hdu_name
         if (hdu_name .eq. 'no') stop
         hdutype = -1
         hdu_extver = 0
         status = 0
         call ftmnhd(unit,hdutype,hdu_name,hdu_extver,status)
      enddo
     
    end subroutine  hdu_move
!---------------------------------------------------------------------!
 
  
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
      
!---------------------------------------------------------------------!
    function dt_f(unit,keyword,status)
      implicit none
      integer :: unit,status
      character (len=30) keyword,comment,error_description
      real dt_f

      call ftgkye(unit,keyword,dt_f,comment,status)

      if(status .ne. 0) then 
         call ftgerr(status,error_description)
         if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
         write(*,*) '   The keyword does not match, try another one'
         read(*,*) keyword
         status = 0 
         call ftgkye(unit,keyword,dt_f,comment,status)
         if(status .ne. 0) then 
            call ftgerr(status,error_description)
            if (status .ne. 0) write(*,*) '!! ATTENTION !! ', error_description
            write(*,*) '   Neither this one, EXIT...'
         endif
      endif
    end function dt_f
!---------------------------------------------------------------------!

!---------------------------------------------------------------------!
      subroutine deletefile(filename,status)
!  A simple little routine to delete a FITS file
      implicit none
      integer status,unit,blocksize
      character*(*) filename
!  Simply return if status is greater than zero
      if (status .gt. 0)return

!  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

!  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
!         file was opened;  so now delete it 
          call ftdelt(unit,status)
          !write(*,*)"Deleted a file"
      else if (status .eq. 103)then
!         file doesn't exist, so just reset status to zero and clear errors
          status=0
          !write(*,*)"Didn't delete a file"
          call ftcmsg
      else
!         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

!  Free the unit number for later reuse
      call ftfiou(unit, status)
    end subroutine deletefile
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      subroutine printerror(status)
!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.
      integer status
      character errtext*30,errmessage*80
!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return
!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext
!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
    end subroutine printerror
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------                            
     function file_line_num(filename)
       implicit none

       character (*), intent(IN) :: filename 
       integer     :: count, file_line_num

       open(1, file = trim(filename))
       count = 0
       do
          read(1, *, END=10)
          count = count + 1
       enddo
10     close(1)
       
       file_line_num = count
     end function file_line_num
!-----------------------------------------------------------------------      


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

!-----------------------------------------------------------------------
      function power2(number)
! Returns TRUE if the number is a power of 2, if not it returns FALSE 
        implicit none
        integer :: number,x
        logical :: power2
        x = number
           ! write(*,*) x
        if (x .eq. 0)then 
           power2 = .false.
           return
        endif
        
        do while(MOD(x,2) .eq. 0)
           x = x / 2
        end do
        if (x .gt. 1) then
           power2 =.false.
           return
        else
           power2 = .true.
        endif
        
      end function power2
!-----------------------------------------------------------------------


!---------------------------------------------------------------------!
     function p_noise_pw(pw, dim, freq, freq_limit)
       implicit none 
       integer, intent(IN) :: dim
       real   , intent(IN) :: pw(dim), freq(dim), freq_limit
       real                :: p_noise_pw

       integer             :: i, count
       real                :: mean 

       count = 0
       mean = 0.0

       do i = 1, dim
          if (freq(i) .ge. freq_limit) then 
             mean = mean + pw(i)
             count = count + 1
          endif 
       enddo
       ! write(*,*) 'p noise: sum, count', mean, count
       
       p_noise_pw = mean / real(count)
     end function p_noise_pw
!---------------------------------------------------------------------!


!------------------------------------------------------------------------
      subroutine ncperiodogram(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim), st(dim) 
        real   , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        real                 :: meanh, means
        real   , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.0
        end do
        call four1(datah, dim, 1)
        call four1(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
!           rc(j) = rc(j) * 2 * dt / (real(dim) * meanh * means) 
           rc(j) = rc(j) * 2 * dt / (real(dim)) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
!           ic(j) = ic(j) *  2 * dt / (real(dim) * meanh * means) 
           ic(j) = ic(j) *  2 * dt / (real(dim)) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
    end subroutine ncperiodogram
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine periodogram(ht, pw, dim)
! Calculates power spectrum between the time series ht(int_len_dim)
! In absolute  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim)
        real   , intent(OUT) :: pw(dim / 2)

        integer              :: j
        real                 :: mean !sum, var
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           ! sum = sum + ht(j)**2
           datah(2 * j - 1) = ht(j)
           datah(2 * j)   = 0.0
        end do

        call four1(datah, dim, 1)
        mean = datah(1) / dim

        do j = 1, dim / 2
!           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * dt / (float(dim) * mean**2) 
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * dt / (float(dim)) 
        end do

!Check the Parcival theorem
!         var = 0.0
!         do j = 0, dim - 1
! !           write(*,*) 2 * j + 1, 2 * j + 2
!            var = var + datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
!            write(111,*) j, datah(2 * j + 1)**2 + datah(2 * j + 2 )**2
!         enddo
!         var = var / dim
!         write(*,*)  'Parcival theorem'
!         write(*,*)  sum, var

!         var = 0.0
!         do j = 1, dim / 2 
!            var = var + pw(j)
!         enddo
!         var = (2 * var + datah(1)**2)  / dim
!         write(*,*)  'Parcival theorem'
!         write(*,*)  sum, var
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine ncperiodogram_frac_rms(ht, st, rc, ic, dim)
! Calculates a cross spectrum between the time series ht(int_len_dim) and st(int_len_dim)
! In fractional  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition, this is S^*(\nu)H(\nu)
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim), st(dim) 
        real   , intent(OUT) :: rc(dim / 2), ic(dim / 2) 
        integer              :: j
        real                 :: meanh, means
        real   , allocatable :: datah(:), datas(:)
        if (.not. allocated(datah)) allocate(datah(2 * dim))
        if (.not. allocated(datas)) allocate(datas(2 * dim))

        do j = 1, dim         
           datah(2 * j - 1) = ht(j)
           datah(2 * j)     = 0.0
           datas(2 * j - 1) = st(j)
           datas(2 * j)     = 0.0
        end do
        call four1(datah, dim, 1)
        call four1(datas, dim, 1)

        meanh = datah(1) / dim
        means = datas(1) / dim
        write(10,*) dim, meanh, means 
        
        do j = 1, dim / 2
           rc(j) = datah(2 * j + 1) * datas(2 * j + 1) + datah(2 * j + 2) * datas(2 * j + 2)
           rc(j) = rc(j) * 2 * dt / (real(dim) * meanh * means) 
           ic(j) = datah(2 * j + 2) * datas(2 * j + 1) - datah(2 * j + 1) * datas(2 * j + 2)
           ic(j) = ic(j) *  2 * dt / (real(dim) * meanh * means) 

        end do
        if (allocated(datah)) deallocate(datah)
        if (allocated(datas)) deallocate(datas)
        return
      end subroutine ncperiodogram_frac_rms
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine periodogram_frac_rms(ht, pw, dim)
! Calculates power spectrum between the time series ht(int_len_dim)
! In fractional  rms normalisation
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        integer, intent(IN)  :: dim
        real   , intent(IN)  :: ht(dim)
        real   , intent(OUT) :: pw(dim / 2)

        integer              :: j
        real                 :: mean !sum, var
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * dim))
        ! sum = 0.0
        

        do j = 1, dim
           datah(2 * j - 1) = ht(j)
           datah(2 * j)   = 0.0
        end do

        call four1(datah, dim, 1)
        mean = datah(1) / dim

        do j = 1, dim / 2
          pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2. * dt / (float(dim) * mean**2) 
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_frac_rms
!------------------------------------------------------------------------


!------------------------------------------------------------------------
      subroutine periodogram_leahy(ht, pw)
! Calculates power spectrum between the time series ht(int_len_dim)
! In leahy normalisation 
! Phase is such that +ve lag corresponds to ht lagging st
! ***MODIFIED*** from Press et al (1992) DFT definition
        use dyn_lc
        implicit none
        real   , intent(IN)  :: ht(int_len_dim)
        real   , intent(OUT) :: pw(int_len_dim / 2)

        integer              :: j
        real                 :: mean
        real   ,allocatable  :: datah(:)

        if (.not. allocated(datah)) allocate(datah(2 * int_len_dim))
        ! sum = 0.0
        do j = 1, int_len_dim
           ! sum = sum + ht(j)**2
           datah(2 * j - 1) = ht(j)
           datah(2 * j)   = 0.0
        end do
        call four1(datah, int_len_dim, 1)
        mean = datah(1) / int_len_dim
        do j = 1, int_len_dim / 2
           pw(j) = (datah(2 * j + 1)**2 + datah(2 * j + 2 )**2) * 2.0 / datah(1)  
        end do
        
        if (allocated(datah)) deallocate(datah)
        return
      end subroutine periodogram_leahy
!------------------------------------------------------------------------


!-----------------------------------------------------------------------
      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      REAL data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end
!-----------------------------------------------------------------------




!---------------------------------------------------------------------!
subroutine cross_FT(lc1, lc2, rc, ic, NN, dt)
! This subroutine calculates the cross-spectrum between lc1 and lc2
! lc2 is the complex cojugate light curve (reference light curve)

! re1, im1, re2, im2 are the real and imaginary part to calculare the 
! Fourier transform 
  implicit none 
  integer, intent(IN)    :: NN
  real   , intent(IN)    :: lc1(NN), lc2(NN), dt
  real   , intent(OUT)   :: rc(NN / 2), ic(NN / 2)
  
  integer                :: i
  real   , allocatable   :: re1(:), im1(:), re2(:), im2(:)

  allocate (re1(NN / 2))
  allocate (im1(NN / 2))
  allocate (re2(NN / 2))
  allocate (im2(NN / 2))
 
  call FT_not_fast(lc1, re1, im1, NN)
  call FT_not_fast(lc2, re2, im2, NN)

  do i = 1, NN / 2    
     rc(i) = (re1(i) * re2(i) + im1(i) * im2(i)) *  2 * dt / (real(NN))
     ic(i) = (im1(i) * re2(i) - re1(i) * im2(i)) *  2 * dt / (real(NN))
  enddo

   end subroutine Cross_FT
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
subroutine power_FT(lc1, pw, NN, dt)
! This subroutine calculates the power-spectrum of lc1
! re1, im1 are the real and imaginary part to calculare the Fourier transform 
  implicit none 
  integer, intent(IN)    :: NN
  real   , intent(IN)    :: lc1(NN), dt 
  real   , intent(OUT)   :: pw(NN / 2)
  
  integer                :: i
  real   , allocatable   :: re1(:), im1(:)

  allocate (re1(NN / 2))
  allocate (im1(NN / 2))
 
  call FT_not_fast(lc1, re1, im1, NN)

  do i = 1, NN / 2    
     pw(i) = (re1(i) * re1(i) + im1(i) * im1(i) ) *  2 * dt / (real(NN))
  enddo

end subroutine Power_FT
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
subroutine FT_not_fast(lc, re, im, NN)
! This subroutine calculates the Fourier transfor of a series lc with NN data points 
! re and im are respectively the real and the imaginary part of the FT
! re and im don't store the average of the light curve, so they have NN/2
! The first frequency is nu = 1/(NN*dt) the last is the Nyquist frequency nuNy = 1/(2*dt)
  implicit none 
  integer, intent(IN)    :: NN
  real   , intent(IN)    :: lc(NN)
  real   , intent(INOUT) :: re(NN / 2), im(NN / 2)
  
  integer                :: n, k
  real                   :: arg
  real   , parameter     :: pi = acos(-1.0)
  


  do n = 1, NN / 2  
     re(n) = 0.0 
     im(n) = 0.0 
     
     do k = 1, NN 

        arg = 2.0 * pi * real(k * n) / real(NN)
        re(n) = re(n) + lc(k) * cos(arg)
        im(n) = im(n) + lc(k) * sin(arg)

     enddo
     re(n) = re(n) / real(NN)
     im(n) = im(n) / real(NN)
  enddo

end subroutine FT_not_fast
!---------------------------------------------------------------------!


!-----------------------------------------------------------------------
subroutine SFT(at,nnmax,nn,Af)
  implicit none
  integer nnmax,nn,j,k
  real at(nnmax)
  complex Af(nnmax/2)
  real pi,arg
  pi = acos(-1.0)
  do j = 1,nn/2
     Af(j) = 0.0
     do k = 1,nn
        arg   = 2.0 * pi * real( j * k ) / real(nn)
        Af(j) = Af(j) + at(k) * complex( cos(arg) , sin(arg) )
     end do
     Af(j) = Af(j) / real(nn)
  end do
  return
end subroutine SFT
!-----------------------------------------------------------------------
