MODULE dyn_lc
!---------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!---------------------------------------------------------------------
! int_len_dim: length of the interval (in units of element) -> this is decided by the user (must be a power of 2)
! int_number: number of interval in the light curve -> this is calculated automatically
  implicit none 
  integer              :: dim_lc, dim_GTI, int_number = -1, int_len_dim, gap = -1, check_gap_num = -1, en_num
  integer, parameter   :: int_len_dim_max = 2000000
  real                 :: dt = -1
  logical              :: check_power2
  real   , allocatable :: lc(:), time(:), err_rate(:), bkg(:), start_GTI(:), end_GTI(:)
  real   , allocatable :: lc_int(:,:), time_int(:,:), bkg_int(:,:)
  real   , allocatable :: lc_en(:,:,:), bkg_en(:,:,:)
  real   , allocatable :: time_int_o(:,:,:), lc_en_o(:,:,:,:), bkg_en_o(:,:,:,:)
  integer, allocatable :: split_ind(:), int_len_dim_o(:)
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
     real poidev, xm, pi
     parameter (pi=3.141592654)
     real alxm, em, g, oldm, sq, t, y!, gammln, ran1
     save alxm, g, oldm, sq
      data oldm /-1./

      if (xm .lt. 12) then
        if(xm .ne. oldm) then
          oldm = xm
          g    = exp(-xm)
        end if
        em = -1.
        t  = 1.
 2      em = em + 1.
        t  = t * ran1()
        if( t .gt. g) goto 2
      else
        if(xm .ne. oldm)then
          oldm = xm
          sq   = sqrt(2. * xm)
          alxm = log(xm)
          g    = xm * alxm - gammln( xm + 1.)
        end if
 1      y = tan(pi * ran1())
        em = sq * y + xm
        if(em .lt. 0.) goto 1
        em = int(em)
        t = 0.9*(1. + y**2.) * exp(em * alxm - gammln(em + 1.) - g)
        if(ran1() .gt. t) goto 1
      end if
      poidev = em
      return
    end function poidev
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
      FUNCTION gammln(xx)
        implicit none 
        REAL gammln, xx
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
  character (len = 500) :: filename,  name_base
  integer, parameter    :: seed = -376736
  real   , parameter    :: pi = acos(-1.0)

!General spectral timing analysis
  real                 :: p_noise_pw
  logical              :: yes_no

!Rebinning 
  integer              :: reb_freq_dim, start, end_reb 
  real                 :: reb_fac
  integer, allocatable :: reb_array(:)
  real   , allocatable :: reb_freq(:), lag_reb(:)
  real   , allocatable :: rc_reb(:), ic_reb(:), rc2_reb(:), ic2_reb(:), rc_ic_reb(:), var_rc_reb(:), var_ic_reb(:), deriv_rc_reb(:), deriv_ic_reb(:), error_prop_reb(:)
  

!LAG vs Freq
  real   , allocatable :: lc_freq1(:,:), lc_freq2(:,:), time_lc(:,:)
  real   , allocatable :: freq(:), pw1(:), pw2(:), rc(:), ic(:), &
       re1(:), im1(:), re2(:), im2(:), rc_freq(:,:), ic_freq(:,:), &
       pw_freq1(:,:), pw_freq2(:,:)
  real   , allocatable :: pw_tot1(:), pw_tot2(:), pw2_tot1(:), pw2_tot2(:), &
       rc_tot(:), ic_tot(:), rc2_tot(:), ic2_tot(:), rc_ic_tot(:)

  real                 :: ave_rate_freq1, ave_rate_freq2, df

  !errors
  real   , allocatable :: lag_freq(:), var_rc_tot(:), var_ic_tot(:), &
     deriv_rc(:), deriv_ic(:), error_prop(:), coher2_freq(:), err_cohe_lag(:)
  real   , allocatable :: covariance(:), errA_rc_ic(:), errA_lag(:)
  real                 :: bias2
  
    call set_seed(seed)

!*******************************************************************************************************************************************************************************************************************************************************
!*******************************************!
!*********** LAG-FREQUENCY SPECTRUM ***********!       
!*******************************************!
    
      ! call execute_command_line('ls')
      
      
!GET FIRST LIGHT CURVE  

    filename = '/Users/gullo/Work/BHB_project/CygX1_Ole/ni2636010101_0.01s_1000eV-2000eV.lc'
    write(*,*) 'First light curve upload'
    
    call extract_lc(filename)
    if (allocated(lc)) then 
       call split_lc()
    else
       write(*,*) '    !!! NO LIGHT CURVE !!! '
       write(*,*) '        STOP HERE   '
       stop
    endif

    if (.not. allocated(lc_freq1)   ) allocate(lc_freq1   (int_number, int_len_dim))
    if (.not. allocated(time_lc)) allocate(time_lc(int_number, int_len_dim))
    ave_rate_freq1 = 0.0
    ! ave_bkg_freq1  = 0.0
    do j = 1, int_len_dim
       do i = 1, int_number
          time_lc  (i, j) = time_int(i, j)
          lc_freq1 (i, j) = lc_int  (i, j)
          ! bkg_freq1(i, j) = bkg_int(i, j)
          ave_rate_freq1 = ave_rate_freq1 + lc_int(i, j)
          ! ave_bkg_freq1  = ave_bkg_freq1  + bkg_int(i, j)
       enddo
    enddo
    ave_rate_freq1 = ave_rate_freq1  / real(int_len_dim * int_number)
    ! ave_bkg_freq1  = ave_bkg_freq1   / real(int_len_dim * int_number)

    write(*,*) 
    write(*,*) '*************************************'
    write(*,*) '   Light curve NAME    '          , trim(filename)
    write(*,*) '   Number of intervals '          , int_number
    write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
    write(*,*) '   Actual exposure (sec) '        , int_len_dim * dt * int_number
    write(*,*) '   Average count rate (count/s)'  , ave_rate_freq1
    write(*,*) '*************************************'
    write(*,*) 



!GET SECOND LIGHT CURVE

    filename = '/Users/gullo/Work/BHB_project/CygX1_Ole/ni2636010101_0.01s_5000eV-8000eV.lc'
    write(*,*) 'Second light curve upload'

    call extract_lc(filename)
    write(*,*) 'ciao ciao ciao '
    if (allocated(lc)) then 
       call split_lc()
    else
       write(*,*) '    !!! NO LIGHT CURVE !!! '
       write(*,*) '        STOP HERE   '
       stop
    endif

    if (.not. allocated(lc_freq2)   ) allocate(lc_freq2   (int_number, int_len_dim))

    
    ave_rate_freq2 = 0.0
    ! ave_bkg_freq2  = 0.0
    do j = 1, int_len_dim
       do i = 1, int_number
          lc_freq2 (i, j) = lc_int(i, j)
          ! bkg_freq2(i, j) = bkg_int(i, j)
          ave_rate_freq2 = ave_rate_freq2 + lc_int(i, j)
          ! ave_bkg_freq2  = ave_bkg_freq2  + bkg_int(i, j)
       enddo
    enddo
    ave_rate_freq2 = ave_rate_freq2  / real(int_len_dim * int_number)
    ! ave_bkg_freq2  = ave_bkg_freq2   / real(int_len_dim * int_number)

    write(*,*) 
    write(*,*) '*************************************'
    write(*,*) '   Light curve NAME    '          , trim(filename)
    write(*,*) '   Number of intervals '          , int_number
    write(*,*) '   Length of the intervals (sec) ', int_len_dim * dt
    write(*,*) '   Actual exposure (sec) '        , int_len_dim * dt * int_number
    write(*,*) '   Average count rate (count/s)'  , ave_rate_freq2
    write(*,*) '*************************************'
    write(*,*) 

    !deallocation
    if(allocated(lc)   ) deallocate(lc  )
    if(allocated(time) ) deallocate(time)
    if(allocated(bkg)  ) deallocate(bkg )

    if(allocated(start_GTI)) deallocate(start_GTI)
    if(allocated(end_GTI)  ) deallocate(end_GTI  )
    if(allocated(split_ind)) deallocate(split_ind)

    if(allocated(lc_int)   ) deallocate(lc_int   )
    if(allocated(time_int) ) deallocate(time_int )
    if(allocated(bkg_int)  ) deallocate(bkg_int  )
    gap = -1 
    check_gap_num = -1

! !Frequency array
    df = 0.5 / (dt * int_len_dim)
    if (.not. allocated(freq)) allocate(freq(int_len_dim / 2))
    ! print *, dt, int_len_dim
    ! print *, 'Frequencies'
    do j = 1, int_len_dim / 2 
       freq(j) = real(j) / (dt * int_len_dim)
       ! write(*,*) freq(j)
    enddo


!PRINT the total light curve (for each obs if more than one) and the corresponding power spectrum 
!***************************************************************************************************************************!
 ! if(yes_no('    Do you want to print the total light curve divided in intervals and its power spectrum? ')) then 

    ! name_base = 'freq1_obs_lc.dat'  
    ! write(*,*) 'name file first light curve ', trim(name_base)
    ! open(70, file = trim(name_base))

    ! write(70,*) 'skip on'
    ! do i = 1, int_number
    !    do j = 1, int_len_dim
    !       write(70, *) time_lc(i, j), lc_freq1(i, j)  
    !    enddo
    !    write(70, *) 'no no'
    ! enddo
    ! close(70)

    ! name_base = 'freq2_obs_lc.dat'  
    ! write(*,*) 'name file second light curve ', trim(name_base)
    ! open(70, file = trim(name_base))

    ! write(70,*) 'skip on'
    ! do i = 1, int_number
    !    do j = 1, int_len_dim
    !       write(70, *) time_lc(i, j), lc_freq2(i, j)  
    !    enddo
    !    write(70, *) 'no no'
    ! enddo
    ! close(70)

 !------------------------------------------------------------------------!
    !FOURIER ANALISYS!
    write(*,*) 
    write(*,*) '-------------------------------------------------------'
    write(*,*) 'Start the Fourier analysis'
      if(.not. allocated(rc_freq )) allocate(rc_freq (int_number, int_len_dim / 2))
      if(.not. allocated(ic_freq )) allocate(ic_freq (int_number, int_len_dim / 2))
      if(.not. allocated(pw_freq1)) allocate(pw_freq1(int_number, int_len_dim / 2))
      if(.not. allocated(pw_freq2)) allocate(pw_freq2(int_number, int_len_dim / 2))
      
      if (.not. allocated(rc) ) allocate(rc (int_len_dim / 2))
      if (.not. allocated(ic) ) allocate(ic (int_len_dim / 2))
      if (.not. allocated(pw1)) allocate(pw1(int_len_dim / 2))
      if (.not. allocated(pw2)) allocate(pw2(int_len_dim / 2))

    do i = 1, int_number
       if (check_power2) then 
          call ncperiodogram(lc_freq2(i, :), lc_freq1(i, :) , rc, ic, int_len_dim) !Cross spectrum 
          call periodogram(lc_freq1(i, :), pw1, int_len_dim)
          call periodogram(lc_freq2(i, :), pw2, int_len_dim)

          do jj = 1, int_len_dim / 2
             rc_freq (i, jj) = rc (jj)
             ic_freq (i, jj) = ic (jj)
             pw_freq1(i, jj) = pw1(jj)
             pw_freq2(i, jj) = pw2(jj)
          enddo
          ! write(*,*) 'FFT done of interval ', i 
       else

          allocate(re1(int_len_dim / 2))
          allocate(im1(int_len_dim / 2))
          allocate(re2(int_len_dim / 2))
          allocate(im2(int_len_dim / 2))

          call FT_not_fast(lc_freq1(i, :), re1, im1, int_len_dim)
          call FT_not_fast(lc_freq2(i, :), re2, im2, int_len_dim)

          do jj = 1, int_len_dim / 2    
             rc_freq (i, jj) = (re1(jj) * re2(jj)) + (im1(jj) * im2(jj))
             ic_freq (i, jj) = (im1(jj) * re2(jj)) - (re1(jj) * im2(jj))
             pw_freq1(i, jj) = (re1(jj) * re1(jj)) + (im1(jj) * im1(jj)) 
             pw_freq2(i, jj) = (re2(jj) * re2(jj)) + (im2(jj) * im2(jj))  
          enddo
          ! write(*,*) 'Fourier transform (slow) done of interval', i
          
          deallocate(re1)
          deallocate(im1)
          deallocate(re2)
          deallocate(im2)

       endif

    enddo
    if (allocated(rc) ) deallocate(rc )
    if (allocated(ic) ) deallocate(ic )
    if (allocated(pw1)) deallocate(pw1)
    if (allocated(pw2)) deallocate(pw2)


!***************************************************************************************************************************!

    write(*,*) 

    if(.not. allocated(pw_tot1) ) allocate(pw_tot1 (int_len_dim / 2))
    if(.not. allocated(pw_tot2) ) allocate(pw_tot2 (int_len_dim / 2))
    if(.not. allocated(pw2_tot1)) allocate(pw2_tot1(int_len_dim / 2))
    if(.not. allocated(pw2_tot2)) allocate(pw2_tot2(int_len_dim / 2))

    pw_tot1  = 0.0
    pw_tot2  = 0.0 
    pw2_tot1 = 0.0
    pw2_tot2 = 0.0 
!Write the power spectra of the two light curve and the total
    do jj = 1, int_len_dim / 2
       do i = 1, int_number
          pw_tot1 (jj) = pw_tot1 (jj) + pw_freq1(i, jj)
          pw2_tot1(jj) = pw2_tot1(jj) + pw_freq1(i, jj) * pw_freq1(i, jj)
          pw_tot2 (jj) = pw_tot2 (jj) + pw_freq2(i, jj)
          pw2_tot2(jj) = pw2_tot2(jj) + pw_freq2(i, jj) * pw_freq2(i, jj)
       enddo
    enddo
    pw_tot1  = pw_tot1  / int_number     
    pw2_tot1 = pw2_tot1 / int_number     
    pw_tot2  = pw_tot2  / int_number     
    pw2_tot2 = pw2_tot2 / int_number     

!Print the first PDS 
    ! name_base = 'freq1_PDS.dat'  
    ! write(*,*) 'name file first PDS ', trim(name_base)
    ! open(71, file = trim(name_base))
    ! write(71, *) 'skip on'
    ! write(71, *) 'read serr 1 2 '
    ! do jj = 1, int_len_dim / 2 - 1 
    !    write(71, *) (freq(jj + 1) + freq(jj)) * 0.5 , df, &
    !         pw_tot1(jj), sqrt( (pw2_tot1(jj) - (pw_tot1(jj) * pw_tot1(jj))) / int_number)
    ! enddo
    ! write(71, *) 'no no'
    ! write(71, *) 'log x y on'
    ! close(71)

!Print the second PDS 
    ! name_base = 'freq2_PDS.dat'  
    ! write(*,*) 'name file second PDS ', trim(name_base)
    ! open(71, file = trim(name_base))
    ! write(71, *) 'skip on'
    ! write(71, *) 'read serr 1 2 '
    ! do jj = 1, int_len_dim / 2 - 1 
    !    write(71, *) (freq(jj + 1) + freq(jj)) * 0.5 , df, &
    !         pw_tot2(jj), sqrt( (pw2_tot2(jj) - (pw_tot2(jj) * pw_tot2(jj))) / int_number)
    ! enddo
    ! write(71, *) 'no no'
    ! write(71, *) 'log x y on'
    ! close(71)
!-----------------------------------------------------------------------!    
!Cross Spectrum

    if(.not. allocated(rc_tot)   ) allocate(rc_tot   (int_len_dim / 2))
    if(.not. allocated(ic_tot)   ) allocate(ic_tot   (int_len_dim / 2))
    if(.not. allocated(rc2_tot)  ) allocate(rc2_tot  (int_len_dim / 2))
    if(.not. allocated(ic2_tot)  ) allocate(ic2_tot  (int_len_dim / 2))
    if(.not. allocated(rc_ic_tot)) allocate(rc_ic_tot(int_len_dim / 2))

    rc_tot    = 0.0
    ic_tot    = 0.0
    rc2_tot   = 0.0
    ic2_tot   = 0.0
    rc_ic_tot = 0.0
    
    do jj = 1, int_len_dim / 2
       do i = 1, int_number
          rc_tot   (jj) = rc_tot   (jj) + rc_freq(i, jj)
          ic_tot   (jj) = ic_tot   (jj) + ic_freq(i, jj)
          rc2_tot  (jj) = rc2_tot  (jj) + rc_freq(i, jj) * rc_freq(i, jj)
          ic2_tot  (jj) = ic2_tot  (jj) + ic_freq(i, jj) * ic_freq(i, jj)
          rc_ic_tot(jj) = rc_ic_tot(jj) + rc_freq(i, jj) * ic_freq(i, jj)
       enddo
    enddo

    rc_tot    = rc_tot    / int_number
    ic_tot    = ic_tot    / int_number
    rc2_tot   = rc2_tot   / int_number
    ic2_tot   = ic2_tot   / int_number
    rc_ic_tot = rc_ic_tot / int_number
    
   if(.not. allocated(lag_freq    )) allocate(lag_freq    (int_len_dim / 2)) 
   if(.not. allocated(var_rc_tot  )) allocate(var_rc_tot  (int_len_dim / 2)) 
   if(.not. allocated(var_ic_tot  )) allocate(var_ic_tot  (int_len_dim / 2)) 
   if(.not. allocated(deriv_rc    )) allocate(deriv_rc    (int_len_dim / 2))
   if(.not. allocated(deriv_ic    )) allocate(deriv_ic    (int_len_dim / 2))
   if(.not. allocated(error_prop  )) allocate(error_prop  (int_len_dim / 2))
   if(.not. allocated(coher2_freq )) allocate(coher2_freq (int_len_dim / 2))
   if(.not. allocated(err_cohe_lag)) allocate(err_cohe_lag(int_len_dim / 2))

   if(.not. allocated(covariance  )) allocate(covariance  (int_len_dim / 2))
   if(.not. allocated(errA_rc_ic  )) allocate(errA_rc_ic  (int_len_dim / 2))
   if(.not. allocated(errA_lag    )) allocate(errA_lag    (int_len_dim / 2))

   name_base = 'lag_freq_err_prop_old.dat'  
   write(*,*) 'name file lag vs frequency, errors calculated with propagation error formula on the standard deviation of real and imaginary part ', trim(name_base)
   open(72, file=trim(name_base))
   name_base = 'lag_freq_err_cohe_old.dat'  
   write(*,*) 'name file lag vs frequency, errors calculated with coherence formula ', trim(name_base)
   open(73, file=trim(name_base))
   write(72, *) 'skip on'
   write(72, *) 'read serr 1 2'
   write(73, *) 'skip on'
   write(73, *) 'read serr 1 2'
   
!lag calculation and errors
   do jj = 1, int_len_dim / 2 - 1 
      ! write(*,*) ic_tot(jj), rc_tot(jj)
      lag_freq(jj) = atan2(ic_tot(jj), rc_tot(jj)) / (2 * pi * freq(jj))

!Errors through standard deviation on real and imaginary part      
!lag with propagation formula
      var_rc_tot(jj) = sqrt((rc2_tot(jj) - rc_tot(jj) * rc_tot(jj)) / int_number)
      var_ic_tot(jj) = sqrt((ic2_tot(jj) - ic_tot(jj) * ic_tot(jj)) / int_number)
      ! covariance(jj) = rc_ic_tot(jj) - (rc_tot(jj) * ic_tot(jj))
      covariance(jj) = 0.0 
      deriv_rc(jj) =  (-1. *  ic_tot(jj) / rc_tot(jj)**2) / (1 + (ic_tot(jj) / rc_tot(jj))**2 )
      deriv_ic(jj) =                      (1. / rc_tot(jj))    / (1 + (ic_tot(jj) / rc_tot(jj))**2 )
      error_prop(jj) = sqrt( (deriv_rc(jj)**2 * var_rc_tot(jj) + deriv_ic(jj)**2 * var_ic_tot(jj) + 2. * deriv_rc(jj) * deriv_ic(jj) * covariance(jj)) ) !/ real(int_number) ) 

!bias term
      ! if ((obs_freq_bins(jj) * int_number) .lt. 500 ) then 
      !    bias2 = ((pw_fq_en(jj, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (real(obs_freq_bins(jj) * int_number))
      ! else
      !    bias2 = 0.0 
      ! endif

      bias2 = 0.0 

!Coherence and error on the lag with the coherence formula
      coher2_freq(jj) = (rc_tot(jj)**2 + ic_tot(jj)**2 - bias2) / (pw_tot1(jj) * pw_tot2(jj))

      err_cohe_lag(jj) = sqrt( (1.0 - coher2_freq(jj)) / (2.0 * coher2_freq(jj) * real(int_number) ) ) / (2 * pi * freq(jj)) 

      ! write(10,*) jj, freq(jj), coher2_freq(jj), rc_tot(jj)**2 + ic_tot(jj)**2, bias2

!Adam's formula real and imaginary part
         ! errA_rc_ic(jj) = sqrt (pw_tot2(jj) * (pw_tot1(jj, k) - ( (rc_tot(jj)**2 + ic_tot(jj)**2 - bias2) / (pw_tot2(jj, k) - P_noise_ext_ref(k)) ) ) / (2 * real( int_number) ) )

!Adam's formula lag
         ! errA_lag(jj) = sqrt( pw_fq_en_ref(jj, k) * ( (pw_fq_en(jj, k) / (rc_fq_en(jj, k)**2 + ic_fq_en(jj, k)**2 - bias2) ) - ( 1 / (pw_fq_en_ref(jj, k) - P_noise_ext_ref(k))  ) ) / (2 *  real(obs_freq_bins(jj) * int_number)) ) / (2 * pi * relevant_freq(jj)) 


      write(72, *) (freq(jj + 1) + freq(jj)) * 0.5 ,  df, lag_freq(jj), error_prop(jj)
      write(73, *) (freq(jj + 1) + freq(jj)) * 0.5 ,  df, lag_freq(jj), err_cohe_lag(jj)
   enddo
   write(72, *) 'log x on'
   write(72, *) 'r y -0.1 0.1'
   write(72, *) 'no no'
   close(72)
   write(73, *) 'log x on'
   write(73, *) 'r y -0.1 0.1'
   write(73, *) 'no no' 
   close(73)

   
   
!Rebinning in freqeuncy

   if(.not. allocated(reb_freq) ) allocate(reb_freq (int_len_dim / 2))
   if(.not. allocated(reb_array)) allocate(reb_array(int_len_dim / 2))

   if(yes_no('Do you want to rebin?')) then
!Define the rebinning factor
      write(*,*) 'Choose exponensial rebinning factor'
      read(*,*) reb_fac
      ! reb_fac = 1.2
   
!Call rebinning function for the frequency to detirmine how many bins
! are going to be rebinned.    
! This info is in reb_array
      call Rlog_rebin_power(freq, reb_freq, reb_array, int_len_dim / 2, reb_freq_dim, reb_fac)

      if(.not. allocated(rc_reb) ) allocate(rc_reb  (reb_freq_dim))
      if(.not. allocated(ic_reb) ) allocate(ic_reb  (reb_freq_dim))
      if(.not. allocated(rc2_reb)) allocate(rc2_reb (reb_freq_dim))
      if(.not. allocated(ic2_reb)) allocate(ic2_reb (reb_freq_dim))
      if(.not. allocated(lag_reb)) allocate(lag_reb (reb_freq_dim))
      if(.not. allocated(rc_ic_reb) ) allocate(rc_ic_reb  (reb_freq_dim))

      ! write(*,*) "ciao", int_len_dim / 2, reb_freq_dim

      rc_reb  = 0.0
      ic_reb  = 0.0
      rc2_reb = 0.0
      ic2_reb = 0.0
      rc_ic_reb = 0.0

      do i = 1, int_number
         start   = 0
         end_reb = 0
         do jj = 1, reb_freq_dim
            start   = end_reb + 1 
            end_reb = end_reb + reb_array(jj)

            ! write(*,*) reb_array(jj)
            ! write(*,*) start, end_reb
            do k = start, end_reb
               rc_reb   (jj) = rc_reb   (jj) + rc_freq(i, k)
               ic_reb   (jj) = ic_reb   (jj) + ic_freq(i, k)
               rc2_reb  (jj) = rc2_reb  (jj) + rc_freq(i, k) * rc_freq(i, k)
               ic2_reb  (jj) = ic2_reb  (jj) + ic_freq(i, k) * ic_freq(i, k)
               rc_ic_reb(jj) = rc_ic_reb(jj) + rc_freq(i, k) * ic_freq(i, k)
            enddo

         enddo
      enddo


      if(.not. allocated(var_rc_reb)) allocate(var_rc_reb(reb_freq_dim))
      if(.not. allocated(var_ic_reb)) allocate(var_ic_reb(reb_freq_dim))
      if(.not. allocated(deriv_rc_reb)) allocate(deriv_rc_reb(reb_freq_dim))
      if(.not. allocated(deriv_ic_reb)) allocate(deriv_ic_reb(reb_freq_dim))
      if(.not. allocated(error_prop_reb)) allocate(error_prop_reb(reb_freq_dim))


      do jj = 1, reb_freq_dim
         rc_reb   (jj) = rc_reb   (jj) / real(int_number * reb_array(jj))
         ic_reb   (jj) = ic_reb   (jj) / real(int_number * reb_array(jj))
         rc2_reb  (jj) = rc2_reb  (jj) / real(int_number * reb_array(jj))
         ic2_reb  (jj) = ic2_reb  (jj) / real(int_number * reb_array(jj))
         rc_ic_reb(jj) = rc_ic_reb(jj) / real(int_number * reb_array(jj))

         lag_reb(jj) = atan2(ic_reb(jj), rc_reb(jj)) / (2 * pi * reb_freq(jj))
         ! write(*,*) 're', rc2_reb(jj), rc_reb(jj) * rc_reb(jj)
         ! write(*,*) 'im', ic2_reb(jj), ic_reb(jj) * ic_reb(jj)
         var_rc_reb(jj) = sqrt( (rc2_reb(jj) - (rc_reb(jj) * rc_reb(jj))) / real(int_number))
         var_ic_reb(jj) = sqrt( (ic2_reb(jj) - (ic_reb(jj) * ic_reb(jj))) / real(int_number))

         deriv_rc_reb(jj) =  (-1. *  ic_reb(jj) / rc_reb(jj)**2) / (1 + (ic_reb(jj) / rc_reb(jj))**2 )
         deriv_ic_reb(jj) =                 (1. / rc_reb(jj))    / (1 + (ic_reb(jj) / rc_reb(jj))**2 )
         error_prop_reb(jj) = sqrt( (deriv_rc_reb(jj)**2 * var_rc_reb(jj) + deriv_ic_reb(jj)**2 * var_ic_reb(jj)) / real(int_number) )
         ! write(*,*) deriv_rc_reb(jj),  var_rc_reb(jj),  deriv_ic_reb(jj),   var_ic_reb(jj),   real(int_number)
      enddo

      !print the rebin
      name_base = 'lag_freq_err_prop_rebin.dat'
      write(*,*) 'name file lag vs frequency REBINNED ', trim(name_base)
      open(72, file=trim(name_base))
      write(72, *) 'skip on'
      write(72, *) 'read serr 1 2'
      do jj = 1, reb_freq_dim - 1
         write(72, *) (reb_freq(jj) + reb_freq(jj+1)) * 0.5, (reb_freq(jj) - reb_freq(jj+1)) * 0.5, lag_reb(jj), error_prop_reb(jj)
      enddo
      write(72, *) 'log x on'
      write(72, *) 'r y -0.1 0.1'
      write(72, *) 'no no'
      close(72)

   endif

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
         check_gap_num = dim_GTI
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

      else
!THIS is called if it's not the first time

!Call interpol_split_silent because the new lc needs interpolation. The number of gaps is saved from the first time          
!There is no need to adjust the GTI arrays because the split is already done. 
! The IF statement before calling the subroutine is because if it is false the first time it should be false every time (no gap interpolation!!)
         ! write(*,*) 'gaaaaaaaaaaaaaaaaaap', gap

         if (gap .ne. -1) then 
            call interpol_split_silent(new_dim_GTI)
         endif
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

    integer                :: j, i, k, w, ee,  max_gap
    real                   :: m, q, lc_interpol, m_b, q_b, bkg_interpol, max_gap_sec 
    real   , allocatable   :: new_start_GTI(:), new_end_GTI(:)

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
       write(*,*) '   interpolate (in sec). Remember that dt is ', real(dt)
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

          m = (lc(ee + i) - lc(i - 1)) / (time(ee + i) - time(i - 1))
          q = lc(i - 1) - m * time(i - 1)
          if (allocated(bkg)) then
             write(*,*) 'bkg active'
             m_b = (bkg(ee + i) - bkg(i - 1)) / (time(ee + i) - time(i - 1))
             q_b = bkg(i - 1) - m_b * time(i - 1)
          endif
          
          do k = i, i + ee - 1
             lc_interpol = m * time(k) + q
             lc_interpol = lc_interpol * dt
             ! write(*,*) 'time interpol', time(k), lc_interpol
             
             if (allocated(bkg)) then 
                bkg_interpol  = m_b * time(k) + q_b
                bkg_interpol = bkg_interpol * dt
             endif
!Extracting the interpolation value from a Poisson distribution with the same mean
             lc(k) = poidev(lc_interpol) / dt
             ! write(*,*) 'extract!', lc(k)
             ! if (lc(k) .lt. 0.0) lc(k) = 0.0

             ! lc(k) = m * time(k) + q
             if (allocated(bkg)) then  
                bkg(k)  = poidev(bkg_interpol) / dt
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
         if ((time(i) - time(i - 1)) .gt. 2 * dt  ) then
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
     count = 1 
        do j = split_ind(2 * i - 1), split_ind(2 * i)
           lc_int  (i, count) = lc(j)
           time_int(i, count) = time(j)
           if(allocated(bkg)) bkg_int (i, count) = bkg(j)
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
         if(.not.allocated(time_e)) allocate(time_e(nrow))
         call ftgcve(unit,colnum,frow,felem,nelem,nullval,time_e,anynul,status)         
      else if (datacode .eq. 82) then
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
         if(.not.allocated(rate_e)) allocate(rate_e(nrow))
         call ftgcve(unit,colnum,frow,felem,nelem,nullval,rate_e,anynul,status)         
      else if (datacode .eq. 82) then
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
!      write(*,*) 'dt', dt

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

!CONVERSION FROM DOUBLE PRECISION TO REAL OF TIME_D, START
      first_time_bin = time_d(1)
      if (allocated(time_d)) then
         if(.not. allocated(time_e)) allocate(time_e(nrow))
         do i = 1, nrow
            time_e(i) = real(time_d(i) - first_time_bin)
         enddo
         deallocate(time_d)
      endif

!CONVERSION FROM DOUBLE PRECISION TO REAL OF RATE_D, START
      if (allocated(rate_d)) then
         if(.not. allocated(rate_e)) allocate(rate_e(nrow))
         do i = 1, nrow
            rate_e(i) = real(rate_d(i))
         enddo
         deallocate(rate_d)
      endif

!CONVERSION FROM DOUBLE PRECISION TO REAL OF BKG_D, START
      if (allocated(bkg_d)) then
         if(.not. allocated(bkg_e)) allocate(bkg_e(nrow))
         do i = 1, nrow
            bkg_e(i) = real(bkg_d(i))
         enddo
         deallocate(bkg_d)
      endif

      first_time_bin = start_GTI_d(1)
      write(*,*) 'first_time_bin', first_time_bin
      
      if (allocated(start_GTI_d)) then
         if(.not. allocated(start_GTI_e)) allocate(start_GTI_e(nrow_GTI))
         do i = 1, nrow_GTI
            start_GTI_e(i) = real(start_GTI_d(i) - first_time_bin)
         enddo
         deallocate(start_GTI_d)
      endif

       if (allocated(end_GTI_d)) then
         if(.not. allocated(end_GTI_e)) allocate(end_GTI_e(nrow_GTI))
         do i = 1, nrow_GTI
            end_GTI_e(i) = real(end_GTI_d(i) - first_time_bin)
         enddo
         deallocate(end_GTI_d)
      endif

!HERE it is possible to print the GTI and the gaps
      ! write(*,*) '-----------------------------------------'
      ! write(*,*) 'GTI: ', nrow_GTI
      ! do i = 1, nrow_GTI - 1
      !    write(*,*) start_GTI_e(i), end_GTI_e(i)
      !    write(*,*) 'ok time ',  end_GTI_e(i) -  start_GTI_e(i)
      !    write(*,*) 'gap ', start_GTI_e(i + 1) - end_GTI_e(i)
      ! enddo
      ! write(*,*) start_GTI_e(nrow_GTI), end_GTI_e(nrow_GTI)
      ! write(*,*) '-----------------------------------------'

!HERE it is possible to print the light curve before we fill the gaps

      ! do i=1, nrow
      !    write(99,*) time_e(i),rate_e(i)
      ! enddo

      ! do i=1, nrow_GTI
      !    write(80,*) start_GTI_e(i),end_GTI_e(i)
      ! enddo
      
!Allocation of the general array
      if(.not. allocated(lc) ) then 
         allocate(lc (nrow))
      else 
         if (nrow .ne. dim_lc) then 
            write(*,*) '   ATTENTION!! The light curves do not have the same length!'
            stop 
         endif
      endif
      
      if(.not. allocated(time) ) allocate(time (nrow))
      if(.not. allocated(err_rate)) allocate(err_rate(nrow))
      if(.not. allocated(bkg)  ) allocate(bkg  (nrow))

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

      if (allocated(time_e)) then 
         do i = 1, nrow
            time(i) = time_e(i)
         enddo
      else 
         deallocate(time)
         write(*,*) '!! ATTENTION !! No TIME column'
      endif

      if (allocated(err_rate_e)) then 
         do i = 1, nrow
            err_rate(i)  = err_rate_e(i)
         enddo
      else 
         deallocate(err_rate)
         write(*,*) '!! ATTENTION !! No ERROR column'
      endif

      if (allocated(bkg_e)) then 
         do i = 1, nrow
            bkg(i)  = bkg_e(i)
         enddo
      else 
         deallocate(bkg)
         write(*,*) '!! ATTENTION !! No BKG column'
      endif
      
!Fill the GTI 
!It is complicated because we have to distinguish between the first call and the others
! and check if the GTIs are the same in all the calls (checking if the GTIs are identical to the previous call)      
      if(.not. allocated(start_GTI)) then  
         allocate(start_GTI(nrow_GTI))
         allocate(end_GTI(nrow_GTI))
!Set the dimension of the GTI and write them in the common arrays      
         dim_GTI = nrow_GTI
!Fill the actual GTI         
         do i = 1, dim_GTI
            start_GTI(i) = start_GTI_e(i)
            end_GTI  (i) = end_GTI_e(i)
         enddo
      else 
         if (check_gap_num .ne. nrow_GTI) then 
            write(*,*) '  ATTENTION!! The GTIs are not the same in all the light curves'
            stop 
         else 
!deallocate the GTIs because they have been modified by the interpolation 
            deallocate(start_GTI)
            deallocate(end_GTI)
!allocate the new GTIs (they will be modify by interpolation)
            allocate(start_GTI(nrow_GTI))
            allocate(end_GTI(nrow_GTI))
            dim_GTI = nrow_GTI
         do i = 1, dim_GTI
            start_GTI(i) = start_GTI_e(i)
            end_GTI  (i) = end_GTI_e(i)
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
           rc(j) = rc(j) * 2 * dt / (real(dim) ) 
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
