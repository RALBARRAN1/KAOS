module MPIReduceSum

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	9 MPI REDUCE SUM:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! ARTIFICIALLY COMPUTE MPI REDUCE SUM COMMAND:

	subroutine MPIReduceSumSub

    ! ----------------------------------------------------

		rrloopMPIRedSum: do rr= 0, ranksize(1)- 1, 1
			if ((rank == 0) .and. (rr == rank)) then
				MPIRedSumOut(1)= MPIRedSumIn(1)
				call mpi_send(MPIRedSumOut, 1, MPI_DOUBLE_PRECISION, rank+ 1, 0, MPI_COMM_WORLD, ierr)
			end if
			if ((rank > 0) .and. (rr == rank)) then
				call mpi_recv(MPIRedSumOut, 1, MPI_DOUBLE_PRECISION, rank- 1, 0, MPI_COMM_WORLD, status, ierr)
				MPIRedSumOut(1)= MPIRedSumOut(1)+ MPIRedSumIn(1)
				if (rr < ranksize(1)- 1) then
					call mpi_send(MPIRedSumOut, 1, MPI_DOUBLE_PRECISION, rank+ 1, 0, MPI_COMM_WORLD, ierr)
					cycle rrloopMPIRedSum
				end if
				if (rr == ranksize(1)- 1) then
					call mpi_send(MPIRedSumOut, 1, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)
					exit rrloopMPIRedSum
				end if
			end if
		end do rrloopMPIRedSum
		ranksizep= ranksize(1)- 1
		if (rank == 0) then
			call mpi_recv(MPIRedSumOut, 1, MPI_DOUBLE_PRECISION, ranksizep, 0, MPI_COMM_WORLD, status, ierr)
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine MPIReduceSumSub

end module MPIReduceSum
