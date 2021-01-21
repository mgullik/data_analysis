include 'header.h'
program fft_lc
  use dyn_lc
  use rand_fun
  implicit none

  integer               :: i
  real                  :: sum, mean 
  real   , allocatable  :: datah(:), re(:), im(:)
  
  character (len = 200) :: filename 

  filename = '/Users/gullo/Work/BHB_project/CygX1_Ole/fft_test_lc.txt'
  open(10, file = filename)

  dim_lc = 8
  allocate(lc(dim_lc))

  do i = 1, dim_lc
     read(10,*) lc(i)
  enddo

  write(*,*) '-----------------------------------'
  sum = 0.0
  do i = 1, dim_lc
     write(*,*) lc(i)
     sum = sum + lc(i)
  enddo
  mean = sum / real(dim_lc)
  write(*,*) 'sum of the time series', sum 
  write(*,*) 'mean of the time series', mean 
  write(*,*) '-----------------------------------'
  
  
  ! allocate(datah(2 * dim_lc))
  ! do i = 1, dim_lc
  !    datah(2 * i - 1) = lc(i)
  !    datah(2 * i)     = 0.0
  ! end do

  ! call four1(datah, dim_lc, 1)
  ! do i = 1, dim_lc
  !    write(*,*) i, datah(2 * i - 1), datah(2 * i)
  ! enddo

  allocate(re(dim_lc / 2))
  allocate(im(dim_lc / 2))
  call FT_not_fast(lc, re, im, dim_lc)

  do i = 1, dim_lc / 2
     write(*,*) i, re(i), im(i)
  enddo
  
  ! call periodogram_frac_rms(lc, pds, dim_lc)

  
  
end program fft_lc
