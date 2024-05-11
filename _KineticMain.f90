! ------------ KAOS (Kinetic model of Auroral ion OutflowS) ------------
! ------------ BY ROBERT M. ALBARRAN II, Ph.D. ------------
! ------------ contact: albarran1@atmos.ucla.edu

! Note: Without cross L-shell convection, solar forcing or other azimuthal asymmetries, particle motion is phi-invariant.
! Change this by not setting phi value in KinetiUpdateA KineticUpdateB KineticUpdateC and KineticRK4Update subroutines.

program KineticMain

use KineticMainParams
use SimParameterization
use DensityProfileA
use DensityProfileB
use DipolePolynomialSolver
use VelocityDistribution
use InitialConditions
use KineticSolver

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

! INITIALIZE MPI FUNCTION:

character(MPI_MAX_PROCESSOR_NAME) :: hostname
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

! ----------------------------------------------------

write(rankstring, '(i10)')  rank ! Convert current rank number to character

! ----------------------------------------------------

! INITIALIZE ALL SIMULATION SUBROUTINES:

call SimParameterizationSub
call DensityProfileASub
call DensityProfileBSub
call DipolePolynomialSolverSub
call VelocityDistributionSub
call InitialConditionsSub
call KineticSolverSub

! ----------------------------------------------------

if (rank == 0) then
	do s= 1, Stot, 1
		do f= 1, SpecieT(s)%NfT(1), 1
			write(*, *)
			call cpu_time(TotEnd) ! End Total CPU time
			write(ToTsGTring, '(i10)')  nint(TotEnd)
			write(*, *) trim('%%%%%% RANK= ' // adjustl(rankstring)) // &
				trim(', TOTAL CPU TIME= ' // adjustl(ToTsGTring)) // &
				trim(' s. KAOS SIMULATION COMPLETE %%%%%%')
		end do
	end do
end if

! ----------------------------------------------------

! FINALIZE MPI FUNCTION:

call MPI_FINALIZE(ierr)

! ----------------------------------------------------

end program KineticMain
