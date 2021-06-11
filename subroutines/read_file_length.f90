  function read_file_length(filename)

  implicit none

  character (len = 500)  :: filename
  integer nlines, read_file_length
  
  nlines = 0 
  open (1, file = filename) 
  do 
     read (1,*, end=10) 
     nlines = nlines + 1 
  enddo 
10 close (1)  

  read_file_length = nlines

end function read_file_length
