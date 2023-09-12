module GaussianRNG

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	3.1- GAUSSIAN RANDOM NUMBER GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE GAUSSIAN RANDOM VARIATES ON [-1 1] FROM UNIFORM VARIATES ON (0 1) BY BOX-MULLER TRANSFORM:

! ----------------------------------------------------

	subroutine GaussianRNGSub

		! ----------------------------------------------------

		call random_number(UniformRN1)
    call random_number(UniformRN2)
		GaussianRN(1)= sqrt(-2d0*log(UniformRN1(1)))*cos(2d0*pi*UniformRN2(1))

		! ----------------------------------------------------

	end subroutine GaussianRNGSub

end module GaussianRNG
