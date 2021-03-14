include 'header.h'
program ls_to_pds
  use dyn_lc
  implicit none

  integer                        :: i, j
  double precision               :: average_rate
  double precision, allocatable  :: pds(:), err_pds(:)
  
  character (len = 500) :: filename, lc_name 
  logical               :: yes_no
  
  ! filename = '/Users/gullo/Work/BHB_project/CygX1_Ole/sinusoid_gaps.fits'
  ! filename = '/Users/gullo/Work/StrayCats/maxij1535/80302312002/products/nu80302312002A01_full_FoV_3to80_01s_sr.lc'
  ! filename = '/Users/gullo/Work/StrayCats/maxij1535/80302312002/products/nu80302312002B01_full_FoV_3to80_01s_sr.lc'

  write(*,*) '-------------------------------------------------' 
  write(*,*) '    CALCULATES THE PDS OF A GIVEN LIGHT CURVE    '
  write(*,*) '-------------------------------------------------'
  write(*,*) 
  write(*,*) '   Enter the name of the light curve with the full path'
  read(*,'(A)') filename
  write(*,*) 
  write(*,*) '    Name: ',filename 
  write(*,*) 
  
  call load_single_lc(filename)


  ! open(10, file = 'Ole_lc_split.qdp')
  ! write(10,*) 'skip on'
  ! do i = 1, int_number
  !    do j = 1, int_len_dim
  !       write(10,*) lc_int(i, j)
  !    enddo
  !    write(10,*) 'no no'
  ! enddo


   
  average_rate = 0.d0
  do i = 1, int_number
     do j = 1, int_len_dim
        average_rate = average_rate + lc_int(i, j)
     enddo
  enddo
  average_rate = average_rate / dble(int_number * int_len_dim)

  write(*,'(A, F8.3)') '   Average rate: ', average_rate 
  
  write(*,*)
  if (yes_no('   Do you want to print the lc?')) then
     lc_name = 'lc_split.qdp'
     open(20, file = lc_name)
     write(20,*) 'skip on'
     do i = 1, int_number
        do j = 1, int_len_dim
           write(20,*) time_int(i, j), lc_int(i, j)
        enddo
        write(20,*) 'no no '
     enddo
     write(20,*) 'scr white'
     write(20,*)
     write(*,*) '   The file name of the light curve is ', trim(lc_name)
  endif
  write(*,*)

  call make_freq_array()

  if(.not. allocated(pds    )) allocate(pds    (int_len_dim/2))
  if(.not. allocated(err_pds)) allocate(err_pds(int_len_dim/2))
  
  call make_pds(freq, lc_int, int_number, dt, int_len_dim, pds, err_pds)

  
end program ls_to_pds
  
  
  
