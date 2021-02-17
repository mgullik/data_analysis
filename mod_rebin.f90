!------------------------------------------------------------------
module rebin 

  !rebin Power Spectrum variables 
  double precision, allocatable :: reb_freq(:), reb_pds(:), &
       reb_err_pds(:)

  !rebin lag vs frequency variables
  double precision, allocatable :: rc_freq_reb(:), ic_freq_reb(:)
  double precision, allocatable :: reb_lag_freq(:), reb_err_cohe_lag(:), reb_err_prop_lag(:)
  double precision, allocatable :: std_rc_freq_reb(:), std_ic_freq_reb(:)
  
contains 
  subroutine rebin_PDS(freq, pds_int, int_number, int_len_dim, reb_dim)
    implicit none

    integer         , intent(IN)  :: int_number, int_len_dim
    double precision, intent(IN)  :: freq(int_len_dim / 2), pds_int(int_number, int_len_dim / 2)

    integer              :: i, j, k, reb_dim, start, end_reb
    double precision     :: reb_fac
    integer, allocatable :: reb_array(:)
    double precision, allocatable :: reb_pds2(:)    


    if(allocated(reb_array  )) deallocate(reb_array)
    if(allocated(reb_freq   )) deallocate(reb_freq)
    if(allocated(reb_pds    )) deallocate(reb_pds )
    if(allocated(reb_err_pds)) deallocate(reb_err_pds)
    if(allocated(reb_pds2   )) deallocate(reb_pds2 )


    write(*,*) '   Choose exponensial rebinning factor'
    read(*,*) reb_fac
    print *, ' '

    allocate(reb_freq (int_len_dim /2))
    allocate(reb_array(int_len_dim /2))

!Finds out the dimension of the rebinned array and the rebin indexes
    call rebin_log(freq, reb_freq, reb_array, int_len_dim / 2, reb_dim, reb_fac)

    allocate(reb_pds    (reb_dim))
    allocate(reb_pds2   (reb_dim))
    allocate(reb_err_pds(reb_dim))

    reb_pds  = 0.d0
    reb_pds2 = 0.d0
       
    do i = 1, int_number
       start   = 0
       end_reb = 0
       do j = 1, reb_dim
          start   = end_reb + 1 
          end_reb = end_reb + reb_array(j)

          do k = start, end_reb
             reb_pds (j) = reb_pds (j) + pds_int(i, k)
             reb_pds2(j) = reb_pds2(j) + pds_int(i, k) * pds_int(i, k)
          enddo
       enddo
    enddo

    do j = 1, reb_dim - 1
       reb_pds (j) = reb_pds (j) / dble(int_number * reb_array(j))
       reb_pds2(j) = reb_pds2(j) / dble(int_number * reb_array(j))
       reb_err_pds(j) = sqrt((reb_pds2(j) - reb_pds(j)**2 ) / dble(int_number * reb_array(j)))
       
    enddo

    if (allocated(reb_array)) deallocate(reb_array)
    if (allocated(reb_pds2 )) deallocate(reb_pds2)

  end subroutine rebin_PDS



  subroutine rebin_lag_freq(freq, rc_freq, ic_freq, pw_freq1, pw_freq2, int_number, int_len_dim, reb_dim)
    implicit none

    integer         , intent(IN)  :: int_number, int_len_dim
    double precision, intent(IN)  :: freq(int_len_dim / 2), rc_freq(int_number, int_len_dim / 2), ic_freq(int_number, int_len_dim / 2), &
         pw_freq1(int_number, int_len_dim / 2), pw_freq2(int_number, int_len_dim / 2)

    integer              :: i, j, k, reb_dim, start, end_reb
    double precision     :: reb_fac, bias2
    integer, allocatable :: reb_array(:)

    double precision, parameter   :: pi = acos(-1.0) 
    double precision, allocatable :: rc2_freq_reb(:), &
         ic2_freq_reb(:), rc_ic_freq_reb(:), pw_freq1_reb(:), &
         pw_freq2_reb(:), pw2_freq1_reb(:), pw2_freq2_reb(:)
    double precision, allocatable :: deriv_rc(:), deriv_ic(:), &
         reb_coher2_freq(:), reb_covariance(:)

    if(allocated(reb_freq        )) deallocate(reb_freq        )
    if(allocated(reb_array       )) deallocate(reb_array       )

    if(allocated(reb_lag_freq    )) deallocate(reb_lag_freq    )
    if(allocated(reb_lag_freq    )) deallocate(reb_lag_freq    )
    if(allocated(reb_err_cohe_lag)) deallocate(reb_err_cohe_lag)
    if(allocated(reb_err_prop_lag)) deallocate(reb_err_prop_lag)
    if(allocated(std_rc_freq_reb )) deallocate(std_rc_freq_reb )
    if(allocated(std_ic_freq_reb )) deallocate(std_ic_freq_reb )
    if(allocated(rc_freq_reb     )) deallocate(rc_freq_reb     )
    if(allocated(ic_freq_reb     )) deallocate(ic_freq_reb     )

    write(*,*) '   Choose exponensial rebinning factor'
    read(*,*) reb_fac
    print *, ' '

    allocate(reb_freq (int_len_dim /2))
    allocate(reb_array(int_len_dim /2))

!Finds out the dimension of the rebinned array and the rebin indexes
    call rebin_log(freq, reb_freq, reb_array, int_len_dim / 2, reb_dim, reb_fac)

    if(.not. allocated(rc_freq_reb)   ) allocate(rc_freq_reb   (reb_dim))
    if(.not. allocated(ic_freq_reb)   ) allocate(ic_freq_reb   (reb_dim))
    if(.not. allocated(rc2_freq_reb)  ) allocate(rc2_freq_reb  (reb_dim))
    if(.not. allocated(ic2_freq_reb)  ) allocate(ic2_freq_reb  (reb_dim))
    if(.not. allocated(rc_ic_freq_reb)) allocate(rc_ic_freq_reb(reb_dim))
    if(.not. allocated(pw_freq1_reb)  ) allocate(pw_freq1_reb  (reb_dim))
    if(.not. allocated(pw_freq2_reb)  ) allocate(pw_freq2_reb  (reb_dim))
    if(.not. allocated(pw2_freq1_reb) ) allocate(pw2_freq1_reb (reb_dim))
    if(.not. allocated(pw2_freq2_reb) ) allocate(pw2_freq2_reb (reb_dim))       

    rc_freq_reb    = 0.d0
    ic_freq_reb    = 0.d0
    rc2_freq_reb   = 0.d0
    ic2_freq_reb   = 0.d0
    rc_ic_freq_reb = 0.d0
    pw_freq1_reb   = 0.d0
    pw_freq2_reb   = 0.d0
    pw2_freq1_reb  = 0.d0
    pw2_freq2_reb  = 0.d0
       
    do i = 1, int_number
       start   = 0
       end_reb = 0
       do j = 1, reb_dim
          start   = end_reb + 1 
          end_reb = end_reb + reb_array(j)
          ! write(*,*) j, start, end_reb
          
          do k = start, end_reb
             rc_freq_reb   (j) = rc_freq_reb   (j) + rc_freq(i, k)
             ic_freq_reb   (j) = ic_freq_reb   (j) + ic_freq(i, k)
             rc2_freq_reb  (j) = rc2_freq_reb  (j) + (rc_freq(i, k) * rc_freq(i, k))
             ic2_freq_reb  (j) = ic2_freq_reb  (j) + (ic_freq(i, k) * ic_freq(i, k))
             rc_ic_freq_reb(j) = rc_ic_freq_reb(j) + (rc_freq(i, k) * ic_freq(i, k))
             pw_freq1_reb  (j) = pw_freq1_reb  (j) + pw_freq1(i, k)
             pw2_freq1_reb (j) = pw2_freq1_reb (j) + (pw_freq1(i, k) * pw_freq1(i, k))
             pw_freq2_reb  (j) = pw_freq2_reb  (j) + pw_freq2(i, k)
             pw2_freq2_reb (j) = pw2_freq2_reb (j) + (pw_freq2(i, k) * pw_freq2(i, k))
          enddo
       enddo
    enddo

    do j = 1, reb_dim - 1
       rc_freq_reb   (j) = rc_freq_reb   (j) / dble(int_number * reb_array(j))
       ic_freq_reb   (j) = ic_freq_reb   (j) / dble(int_number * reb_array(j))
       rc2_freq_reb  (j) = rc2_freq_reb  (j) / dble(int_number * reb_array(j))
       ic2_freq_reb  (j) = ic2_freq_reb  (j) / dble(int_number * reb_array(j))
       rc_ic_freq_reb(j) = rc_ic_freq_reb(j) / dble(int_number * reb_array(j))
       pw_freq1_reb  (j) = pw_freq1_reb  (j) / dble(int_number * reb_array(j))
       pw_freq2_reb  (j) = pw_freq2_reb  (j) / dble(int_number * reb_array(j))
       pw2_freq1_reb (j) = pw2_freq1_reb (j) / dble(int_number * reb_array(j))
       pw2_freq2_reb (j) = pw2_freq2_reb (j) / dble(int_number * reb_array(j))
    enddo

    
    allocate(std_rc_freq_reb (reb_dim)) 
    allocate(std_ic_freq_reb (reb_dim)) 
    allocate(reb_err_prop_lag(reb_dim))
    allocate(reb_err_cohe_lag(reb_dim))
    allocate(reb_coher2_freq (reb_dim))
    allocate(reb_lag_freq    (reb_dim)) 

    allocate(deriv_rc        (reb_dim))
    allocate(deriv_ic        (reb_dim))
    allocate(reb_covariance  (reb_dim))    
   
!lag calculation and errors
   do j = 1, reb_dim - 1

!Errors through standard deviation on real and imaginary part      
      !lag with propagation formula
      std_rc_freq_reb(j) = sqrt((rc2_freq_reb(j) - (rc_freq_reb(j) * rc_freq_reb(j)) ))  / (sqrt(dble(int_number * reb_array(j))))
      std_ic_freq_reb(j) = sqrt((ic2_freq_reb(j) - (ic_freq_reb(j) * ic_freq_reb(j)) ))  / (sqrt(dble(int_number * reb_array(j))))
      reb_covariance(j)  = (rc_ic_freq_reb(j) - ((rc_freq_reb(j) * ic_freq_reb(j) ) ) )  /       dble(int_number * reb_array(j))
      ! covariance(j) = 0.0 
      deriv_rc(j) =  (-1.d0 *  ic_freq_reb(j) / rc_freq_reb(j)**2) / (1 + (ic_freq_reb(j) / rc_freq_reb(j))**2 )
      deriv_ic(j) =                      (1.d0 / rc_freq_reb(j))    / (1 + (ic_freq_reb(j) / rc_freq_reb(j))**2 )

      reb_err_prop_lag(j) = sqrt( (deriv_rc(j)**2 * std_rc_freq_reb(j)**2) + (deriv_ic(j)**2 * std_ic_freq_reb(j)**2) + (2.d0 * deriv_rc(j) * deriv_ic(j) * reb_covariance(j)) ) / (2.d0 * pi * reb_freq(j))

      
!bias term
      ! if ((obs_freq_bins(j) * int_number) .lt. 500 ) then 
      !    bias2 = ((pw_fq_en(j, k) - P_noise_ext(k)) * P_noise_ext_ref(k) + (pw_fq_en_ref(j, k) - P_noise_ext_ref(k)) * P_noise_ext(k) + P_noise_ext(k) * P_noise_ext_ref(k) ) / (dble(int_number))
      ! else
      !    bias2 = 0.0 
      ! endif

      bias2 = 0.0 

!Coherence and error on the lag with the coherence formula
      reb_coher2_freq(j) = (rc_freq_reb(j)**2 + ic_freq_reb(j)**2 - bias2) / (pw_freq1_reb(j) * pw_freq2_reb(j))
      ! write(*,*) 'iiii', pw_freq1_reb(j), pw_freq2_reb(j), reb_coher2_freq(j)

      
      reb_err_cohe_lag(j) = sqrt( (1.d0 - reb_coher2_freq(j)) / (2.d0 * reb_coher2_freq(j) * dble(int_number * reb_array(j)))) / (2.d0 * pi * reb_freq(j)) 

      ! write(*,*) 'iiiii', reb_err_cohe_lag(j)
      ! write(*,*) j, freq(j), reb_coher2_freq(j), rc_freq_reb(j)**2, ic_freq_reb(j)**2, pw_freq1_reb(j),  pw_freq2_reb(j)

!Adam's formula real and imaginary part
         ! errA_rc_ic(j) = sqrt (pw_freq_reb2(j) * (pw_freq_reb1(j, k) - ( (rc_freq_reb(j)**2 + ic_freq_reb(j)**2 - bias2) / (pw_freq_reb2(j, k) - P_noise_ext_ref(k)) ) ) / (2 * dble( int_number * reb_dim) ) )

!Adam's formula lag
         ! errA_lag(j) = sqrt( pw_fq_en_ref(j, k) * ( (pw_fq_en(j, k) / (rc_fq_en(j, k)**2 + ic_fq_en(j, k)**2 - bias2) ) - ( 1 / (pw_fq_en_ref(j, k) - P_noise_ext_ref(k))  ) ) / (2 *  dble(int_number * reb_dim)) ) / (2 * pi * relevant_freq(j)) 

      reb_lag_freq(j) = atan2(ic_freq_reb(j), rc_freq_reb(j)) / (2.d0 * pi * reb_freq(j))

   enddo
    
   deallocate(reb_array     )
   deallocate(deriv_rc      )
   deallocate(deriv_ic      )
   deallocate(reb_covariance)
   deallocate(rc_freq_reb   )
   deallocate(ic_freq_reb   )
   deallocate(rc2_freq_reb  )
   deallocate(ic2_freq_reb  )
   deallocate(rc_ic_freq_reb)   
   deallocate(pw_freq1_reb  )
   deallocate(pw_freq2_reb  )
   deallocate(pw2_freq1_reb )
   deallocate(pw2_freq2_reb )
   
    
  end subroutine rebin_lag_freq


  
  subroutine rebin_PDS_FPMA_B(freq, pds_intA, pds_intB, int_number, int_len_dim, reb_dim)
    implicit none

    integer         , intent(IN)  :: int_number, int_len_dim
    double precision, intent(IN)  :: freq(int_len_dim / 2), &
         pds_intA(int_number, int_len_dim / 2), &
         pds_intB(int_number, int_len_dim / 2)

    integer              :: i, j, k, reb_dim, start, end_reb
    double precision     :: reb_fac
    integer, allocatable :: reb_array(:)
    double precision, allocatable :: reb_pds2(:)    


    if(allocated(reb_array  )) deallocate(reb_array)
    if(allocated(reb_freq   )) deallocate(reb_freq)
    if(allocated(reb_pds    )) deallocate(reb_pds )
    if(allocated(reb_err_pds)) deallocate(reb_err_pds)
    if(allocated(reb_pds2   )) deallocate(reb_pds2 )


    write(*,*) '   Choose exponensial rebinning factor'
    read(*,*) reb_fac
    print *, ' '

    allocate(reb_freq (int_len_dim /2))
    allocate(reb_array(int_len_dim /2))

!Finds out the dimension of the rebinned array and the rebin indexes
    call rebin_log(freq, reb_freq, reb_array, int_len_dim / 2, reb_dim, reb_fac)

    
    allocate(reb_pds    (reb_dim))
    allocate(reb_pds2   (reb_dim))
    allocate(reb_err_pds(reb_dim))

    reb_pds  = 0.d0
    reb_pds2 = 0.d0

    do i = 1, int_number
       start   = 0
       end_reb = 0
       do j = 1, reb_dim
          start   = end_reb + 1 
          end_reb = end_reb + reb_array(j)

          do k = start, end_reb
             reb_pds (j) = reb_pds (j) + pds_intA(i, k)
             reb_pds (j) = reb_pds (j) + pds_intB(i, k)
             reb_pds2(j) = reb_pds2(j) + pds_intA(i, k) **2
             reb_pds2(j) = reb_pds2(j) + pds_intB(i, k) **2
          enddo
       enddo
    enddo

    do j = 1, reb_dim
       reb_pds (j) = reb_pds (j) / dble(2 * int_number * reb_array(j))
       reb_pds2(j) = reb_pds2(j) / dble(2 * int_number * reb_array(j))
       reb_err_pds(j) = sqrt((reb_pds2(j) - reb_pds(j)**2 ) / dble(2 * int_number * reb_array(j)))
       
    enddo

    if (allocated(reb_array)) deallocate(reb_array)
    if (allocated(reb_pds2 )) deallocate(reb_pds2)

  end subroutine rebin_PDS_FPMA_B

  
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
!----------------------------------------------------------------
!----------------------------------------------------------------
      subroutine rebin_log(vector, rebin_vector, rebin_array, dim, nf, c)
        implicit none

        integer         , intent(IN)    :: dim
        integer         , intent(INOUT) :: nf
        integer         , intent(OUT)   :: rebin_array(dim)
        double precision, intent(IN)    :: vector(dim), c 
        double precision, intent(OUT)   :: rebin_vector(dim)

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
           rebin_vector(j + 1)  = 0.d0
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
      end subroutine rebin_log
 
  
end module rebin
!------------------------------------------------------------------
