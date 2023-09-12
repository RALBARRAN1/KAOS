module ExponentialRNG

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	3.2- EXPONENTIAL RANDOM NUMBER GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE EXPONENTIAL RANDOM VARIATES FROM UNIFORM VARIATES:

! ----------------------------------------------------

	subroutine ExponentialRNGSub

		! ----------------------------------------------------

		ExponentialRN(1)= -1d0
		do while ((ExponentialRN(1) < 0d0) .or. (ExponentialRN(1) > 1d0))
      call random_number(UniformRN)
      ExponentialRN(1)= log(1d0- UniformRN(1))/(-lambdaexp)
		end do

    ! ----------------------------------------------------

    ! DIAGNOSTIC FLAGS FOR EXPONENTIAL VARIATES WITHIN BOUNDS:

		if ((ExponentialRN(1) < 0d0) .or. (ExponentialRN(1) > 1d0)) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ExponentialRN= ', &
        ExponentialRN(1), ' IS OUTSIDE (0, 1) BOUNDS', &
        ' IN EXPONENTIAL RNG SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

	end subroutine ExponentialRNGSub

end module ExponentialRNG
