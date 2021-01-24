subroutine load_filename_lag_ene()
  implicit none

  character (len=200) :: name_path, filename_ch

  call execute_command_line('ls')
  write(*,*)
  
  if (yes_no('   Have you created the initial file init.txt? ')) then
     goto 666
  else
     write(*,*) '    Create the file init.txt with this structure:'
     write(*,*) '      line 1 path name: the full path of the folder where the light curves are (without / at the end)'
     write(*,*) '      line 2 filename of the energy bins: this is where the energy bins are '
     write(*,*)
     write(*,*)
     
     

  
666  continue
     

end subroutine load_filename_lag_ene
