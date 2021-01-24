MODULE dyn_lc
!------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!-----------------------------------------------------------------
! int_len_dim: length of the interval (in units of element) -> this is decided by the user (must be a power of 2)
! int_number: number of interval in the light curve -> this is calculated automatically
  implicit none 
  integer     :: dim_lc, dim_GTI, int_number = -1
  integer     :: int_len_dim, gap = -1, check_gap_num = -1, en_num
  real                 :: dt = -1
  logical              :: check_power2

  real   , allocatable :: lc(:), time(:), err_rate(:), bkg(:)
  real   , allocatable :: start_GTI(:), end_GTI(:)

!split lc variables   
  real   , allocatable :: lc_int(:,:), time_int(:,:), bkg_int(:,:)

  real                 :: df
  real   , allocatable :: freq(:)
  
!lag vs freq variables
  real                 :: ave_rate_freq1, ave_rate_freq2
  real   , allocatable :: lc_freq1(:,:), lc_freq2(:,:)
  ! real   , allocatable :: pw_freq1_ave(:), pw_freq2_ave(:), &
  !                         pw2_freq1_ave(:), pw2_freq2_ave(:)


!lag vs energy variables
  real                 :: ave_rate_ref
  real   , allocatable :: l_en_bin(:), r_en_bin(:), ave_rate_en(:)
  real   , allocatable :: lc_en(:,:,:), bkg_en(:,:,:), lc_ref(:,:)
  real   , allocatable :: pw_ref(:), err_pw_ref(:)
  integer              :: freq_num
  integer, allocatable :: lower_fq(:), upper_fq(:)
  real   , allocatable :: relevant_freq(:)
  ! real   , allocatable :: time_int_o(:,:,:), lc_en_o(:,:,:,:), bkg_en_o(:,:,:,:)
  integer, allocatable :: split_ind(:), int_len_dim_o(:)
END MODULE dyn_lc
!------------------------------------------------------------------
