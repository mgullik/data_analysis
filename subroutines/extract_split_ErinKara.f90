!this subroutine substitutes the extract and the split routine in the specific case of Erin's light curves

subroutine extract_split_erin(filename, o)
  use dyn_lc
  implicit none

  character (len = 500)  :: filename
  integer, INTENT(IN)    :: o
  integer :: j,  read_file_length
  double precision :: temp, temp_time

  int_len_dim = read_file_length(filename) - 1
  int_number = 1

  ! write(*,*) 'length of the light curve', int_len_dim
  open(2, file = filename)

  read(2,*) !skip the first line which is the header 
  do j = 1, int_len_dim
     ! read(2, *) temp, time_obs(1, j, obs_index),
     read(2,*)  temp_time, lc_en_obs(1,j,1,o), temp, lc_en_obs(1,j,2,o), temp, lc_en_obs(1,j,3,o), temp, lc_en_obs(1,j,4,o), temp, lc_en_obs(1,j,5,o), temp, lc_en_obs(1,j,6,o), temp, lc_en_obs(1,j,7,o), temp, lc_en_obs(1,j,8,o), temp, lc_en_obs(1,j,9,o), temp, lc_en_obs(1,j,10,o), temp, lc_en_obs(1,j,11,o), temp, lc_en_obs(1,j,12,o), temp, lc_en_obs(1,j,13,o), temp, lc_en_obs(1,j,14,o), temp, lc_en_obs(1,j,15,o), temp, lc_en_obs(1,j,16,o), temp, lc_en_obs(1,j,17,o), temp, lc_en_obs(1,j,18,o), temp, lc_en_obs(1,j,19,o), temp, lc_en_obs(1,j,20,o), temp, lc_en_obs(1,j,21,o), temp, lc_en_obs(1,j,22,o), temp, lc_en_obs(1,j,23,o), temp, lc_en_obs(1,j,24,o), temp, lc_en_obs(1,j,25,o), temp, lc_en_obs(1,j,26,o), temp, lc_en_obs(1,j,27,o), temp, lc_en_obs(1,j,28,o), temp, lc_en_obs(1,j,29,o), temp, lc_en_obs(1,j,30,o), temp, lc_en_obs(1,j,31,o), temp, lc_en_obs(1,j,32,o), temp, lc_en_obs(1,j,33,o), temp, lc_en_obs(1,j,34,o), temp, lc_en_obs(1,j,35,o), temp
     ! write(*,*) j, temp_time

  enddo
  close(2)

  
 end subroutine extract_split_erin
