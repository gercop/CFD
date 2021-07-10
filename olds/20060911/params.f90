      module params

      implicit none

      save      
        integer, parameter, public      :: ghost_points=3

        !General MPI variables
        integer, dimension(4), public   :: local_sd
        integer :: local_imax,local_bx,global_x1,global_x2
        integer :: clock_start,clock_end,clock_rate

        integer, public                 :: my_rank,pro,ierror,status(MPI_STATUS_SIZE)
        

      end module params
