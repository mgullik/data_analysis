!------------------------------------------------------------------
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
!----------------------------------------------------------------- 
