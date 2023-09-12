module MBVelocityDistribution

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	6.1- MB VELOCITY DISTRIBUTION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use GaussianRNG

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE 3D MAXWELL-BOLTZMANN VELOCITY DISTRIBUTION IN CARTESIAN BY PROJECTING RANDOM NUMBERS FROM
! UNIFORM DISTRIBUTION ONTO NORMAL DISTRIBUTION BY BOX-MULLER TRANSFORM:

! ----------------------------------------------------

	subroutine MBVelocityDistributionSub

		! ----------------------------------------------------

		! Compute Gaussian random variates
		UniformRN1(1)= Vxrandn1(1)
		UniformRN2(1)= Vxrandn2(1)
		call GaussianRNGSub
		Vx0pp(1)= GaussianRN(1)

		UniformRN1(1)= Vyrandn1(1)
		UniformRN2(1)= Vyrandn2(1)
		call GaussianRNGSub
		Vy0pp(1)= GaussianRN(1)

		UniformRN1(1)= Vzrandn1(1)
		UniformRN2(1)= Vzrandn2(1)
		call GaussianRNGSub
		Vz0pp(1)= GaussianRN(1)

		! MB distribution standard deviations (Global Cartesian coords.)
		sigmaVx0(1)= sqrt(kB*(TsMB(1))/(SpecieT(s)%msT(1))) ! Where TparX= TparY= TparZ= Tpar= Ti
		sigmaVy0(1)= sigmaVx0(1)
		sigmaVz0(1)= sigmaVx0(1)

		! MB Cartesian velocity components
		Vx0p(1)= sigmaVx0(1)*Vx0pp(1)
		Vy0p(1)= sigmaVy0(1)*Vy0pp(1)
		Vz0p(1)= sigmaVz0(1)*Vz0pp(1)

		! ----------------------------------------------------

		! MB Dipole velocity components
		Vp0p(1)= ((1d0- 3d0*((cos(thetaMB(1)))**2d0))/sqrt(ellMB(1)))*(Vx0p(1)*cos(phiMB(1))+ Vy0p(1)*sin(phiMB(1)))+ &
			((3d0*Vz0p(1)*cos(thetaMB(1))*sin(thetaMB(1)))/sqrt(ellMB(1)))
		Vq0p(1)= ((3d0*cos(thetaMB(1))*sin(thetaMB(1)))/sqrt(ellMB(1)))*(Vx0p(1)*cos(phiMB(1))+ Vy0p(1)*sin(phiMB(1)))+ &
			((Vz0p(1)*(3d0*((cos(thetaMB(1)))**2d0)- 1d0))/sqrt(ellMB(1)))
		Vphi0p(1)= -Vx0p(1)*sin(phiMB(1))+ Vy0p(1)*cos(phiMB(1))

		! ----------------------------------------------------

		! SET INITIAL PERPENDICULAR VELOCITIES FOR ALL PARTICLES:
		! Note: Let Vperp1= Vp= Vy AND Vperp2= Vphid= Vx

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 0) then
			Vperp1MB(1)= Vp0p(1)
			Vperp2MB(1)= Vphi0p(1)
			VperpMB(1)= abs(sqrt(Vperp1MB(1)**2d0+ Vperp2MB(1)**2d0))
		end if

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 1) then
			Vperp1MB(1)= 0d0
			Vperp2MB(1)= 0d0
			VperpMB(1)= 0d0
		end if

		! ----------------------------------------------------

		! SET INITIAL PARALLEL VELOCITIES FOR ALL PARTICLES:
		! Note: Let Vpar= Vq= Vz

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 0) then
	    VparMB(1)= Vq0p(1)
		end if

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 1) then
			VparMB(1)= 0d0
		end if

		! ----------------------------------------------------

		! SET INITIAL CARTESIAN VELOCITIES FOR ALL PARTICLES:
		! Note: Translational velocities are from parallel velocity component only.

		! Cartesian velocities [m/s]

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 0) then
			VxMB(1)= cos(phiMB(1))*(3d0*VparMB(1)*cos(thetaMB(1))*sin(thetaMB(1)))/sqrt(ellMB(1))
			VyMB(1)= sin(phiMB(1))*(3d0*VparMB(1)*cos(thetaMB(1))*sin(thetaMB(1)))/sqrt(ellMB(1))
			VzMB(1)= VparMB(1)*(3d0*((cos(thetaMB(1)))**2d0)- 1d0)/sqrt(ellMB(1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1) == 1) then
			VxMB(1)= 0d0
			VyMB(1)= 0d0
			VzMB(1)= 0d0
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine MBVelocityDistributionSub

end module MBVelocityDistribution
