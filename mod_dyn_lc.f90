MODULE dyn_lc
!------------------------------------------------------------------
!  Module containing definitions needed to dynamically allocate 
!  the values of an array 
!-----------------------------------------------------------------
! int_len_dim: length of the interval (in units of element) -> this is decided by the user (must be a power of 2)
! int_number: number of interval in the light curve -> this is calculated automatically
  implicit none 
  integer     :: dim_lc, dim_GTI, int_number = -1
  integer     :: int_len_dim, gap = -1, check_gap_num = -1, en_num, obs_num
  double precision     :: dt = -1.d0
  logical              :: check_power2, check_merge, verbose_merge

  integer         , parameter   :: max_GTI_dim = 1000, int_number_max = 100, int_len_dim_max = 100000

  double precision, allocatable :: lc(:), time(:), err_rate(:), bkg(:)
  double precision, allocatable :: start_GTI(:), end_GTI(:)

!split lc variables
  integer                       :: max_gap_sec_init
  double precision, allocatable :: lc_int(:,:), time_int(:,:)

  double precision              :: df
  double precision, allocatable :: freq(:), freq_obs(:,:), df_obs(:)
  
!lag vs freq variables
  double precision              :: ave_rate_freq1, ave_rate_freq2
  double precision, allocatable :: lc_freq1(:,:), lc_freq2(:,:)
  ! real   , allocatable :: pw_freq1_ave(:), pw_freq2_ave(:), &
  !                         pw2_freq1_ave(:), pw2_freq2_ave(:)

!lag vs freq multi obs
  double precision, allocatable :: lc_freq1_obs(:,:,:), lc_freq2_obs(:,:,:), ave_rate_freq1_obs(:), ave_rate_freq2_obs(:)
  
  double precision              :: starting_telescope_time

!lag vs energy variables
  double precision              :: ave_rate_ref
  double precision, allocatable :: l_en_bin(:), r_en_bin(:),&
  ave_rate_en(:)
  double precision, allocatable :: lc_en(:,:,:),  lc_ref(:,:)
  double precision, allocatable :: pw_ref(:), err_pw_ref(:)
  integer                       :: freq_num
  integer         , allocatable :: lower_fq(:), upper_fq(:), lower_fq_obs(:,:), upper_fq_obs(:,:)
  double precision, allocatable :: relevant_freq(:)
  integer         , allocatable :: split_ind(:), int_len_dim_o(:)

!lag vs energy multi obs
  integer         , allocatable :: int_number_obs(:), int_len_dim_obs(:)
  double precision, allocatable :: lc_en_obs(:,:,:,:)
  
END MODULE dyn_lc
!------------------------------------------------------------------
