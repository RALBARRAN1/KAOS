module IonNeutralCollisions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.1.1 ION-NEUTRAL COLLISIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! PERFORM CHARGE EXCHANGE BY ION-NEUTRAL COLLISIONS:

	subroutine IonNeutralCollisionsSub

		! ----------------------------------------------------

		! COMPUTE ENA VELOCITY COMPONENTS:

		! ----------------------------------------------------

		ENAflag(j)= .true.
		ENAflagN1ind(j)= n

		! Compute Dipole ENA velocity components
		Vpp(1)= Vperp1N(j)
		Vphip(1)= Vperp2N(j)
		Vparp(1)= 3d0*cos(thetak4(j))*sin(thetak4(j))*(VxN(j)*cos(phik4(j))+ &
			VyN(j)*sin(phik4(j)))/sqrt(ellk4(j))+ VzN(j)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/ &
			sqrt(ellk4(j))

		! Compute Cartesian ENA velocity components
		VxN(j)= cos(phik4(j))*(3d0*Vparp(1)*cos(thetak4(j))*sin(thetak4(j))+ &
			Vpp(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))- &
			Vphip(1)*sin(phik4(j))
		VyN(j)= sin(phik4(j))*(3d0*Vparp(1)*cos(thetak4(j))*sin(thetak4(j))+ &
			Vpp(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))+ &
			Vphip(1)*cos(phik4(j))
		VzN(j)= Vparp(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))+ &
			3d0*Vpp(1)*cos(thetak4(j))*sin(thetak4(j))/sqrt(ellk4(j))
		Vperp1N(j)= 0d0
		Vperp2N(j)= 0d0
		VperpN(j)= 0d0

		! ----------------------------------------------------

	end subroutine IonNeutralCollisionsSub

end module IonNeutralCollisions
