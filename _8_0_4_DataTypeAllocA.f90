module DataTypeAllocA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.4- DATA TYPE ALLOCATION A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! ALLOCATE KINETIC SOLVER DERIVED DATA TYPES:

	subroutine DataTypeAllocASub

    ! ----------------------------------------------------

    ! ALLOCATE PARTICLE-DEPENDENT DERIVED DATA TYPES:

    ! ----------------------------------------------------

    ! ALLOCATE KINETIC RK4 UPDATE VARIABLES:

		allocate(LBoutfluxIonMsk(NsTK(1)), UBoutfluxIonMsk(NsTK(1)), AEAmagN(NsTK(1)), AGmagN(NsTK(1)), AEPmagN(NsTK(1)))

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			allocate(LBoutfluxENAMsk(NsTK(1)), UBoutfluxENAMsk(NsTK(1)))
		end if

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			allocate(LBreplenishIonMsk(NsTK(1)))
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			allocate(UBreplenishIonMsk(NsTK(1)))
		end if

		allocate(Qindk1(NsTK(1)), Qindk0(NsTK(1)), Vparindk1(NsTK(1)), &
			Vparindk0(NsTK(1)), Vpindk1(NsTK(1)), Vqindk1(NsTK(1)), Vphiindk1(NsTK(1)))

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			allocate(Vperp1indk1(NsTK(1)), Vperp2indk1(NsTK(1)), &
				Vperp1indk0(NsTK(1)), Vperp2indk0(NsTK(1)))

		else

			allocate(Vperpindk1(NsTK(1)), Vperpindk0(NsTK(1)))

		end if

		allocate(VxN(NsTK(1)), VyN(NsTK(1)), VzN(NsTK(1)), xN(NsTK(1)), yN(NsTK(1)), zN(NsTK(1)), &
			Vperp1N(NsTK(1)), Vperp2N(NsTK(1)), VperpN(NsTK(1)))
		allocate(pk4(NsTK(1)), phik4(NsTK(1)), ellk4(NsTK(1)), qk4(NsTK(1)), rk4(NsTK(1)), thetak4(NsTK(1)))
		allocate(ENAflag(NsTK(1)), ENAflagN0ind(NsTK(1)), ENAflagN1ind(NsTK(1)))

    ! ----------------------------------------------------

    ! ALLOCATE KINETIC UPDATE A VARIABLES:

    allocate(AEAmag(NsTK(1)), AGmag(NsTK(1)), AEPmag(NsTK(1)), x(NsTK(1)), y(NsTK(1)), z(NsTK(1)), &
      Vperp1(NsTK(1)), Vperp2(NsTK(1)), Vperp(NsTK(1)), Vpar(NsTK(1)), &
      Vx(NsTK(1)), Vy(NsTK(1)), Vz(NsTK(1)), Vphi(NsTK(1)), Vp(NsTK(1)))

    ! ----------------------------------------------------

    ! ALLOCATE ION PARTICLE COUNTS VARIABLES IN DERIVED DATA TYPES:

    allocate(pdriftion(NsTK(1)), qdriftion(NsTK(1)), phidriftion(NsTK(1)), &
      rdriftENA(NsTK(1)), thetadriftENA(NsTK(1)), phidriftENA(NsTK(1)))

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DataTypeAllocASub

end module DataTypeAllocA
