module KineticRK4Update

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.2 KINETIC RK4 UPDATE:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use KineticUpdateA
use KineticUpdateB
use KineticUpdateC
use RK4DipolePolynomialSolver

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! UPDATE CARTESIAN KINETIC VARIABLES BY RK4 INTEGRATION:

	subroutine KineticRK4UpdateSub

		! ----------------------------------------------------

		call KineticUpdateASub

		! ----------------------------------------------------

		! CONTINUE WITH RK4 METHOD AT HALF TIME-STEP FORWARD FOR ALL PARTICLES:

		k1Vx(1)= Axk1(1) ! Compute k1 values
		k1Vy(1)= Ayk1(1)
		k1Vz(1)= Azk1(1)

		! Compute x, y, z half time-step forward

		xk2(1)= x(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*VxR(1)
		yk2(1)= y(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*VyR(1)
		zk2(1)= z(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*VzR(1)

		! ----------------------------------------------------

		call KineticUpdateBSub

		! ----------------------------------------------------

		! CONTINUE WITH RK4 METHOD AT FULL TIME-STEP FORWARD:

		k2Vx(1)= Axk2(1) ! Compute k2 values
		k2Vy(1)= Ayk2(1)
		k2Vz(1)= Azk2(1)
		k3Vx(1)= Axk2(1) ! Compute k3 values (Note: a_i ~= a_i(v_i) for i= x, y, z, par)
		k3Vy(1)= Ayk2(1)
		k3Vz(1)= Azk2(1)

		! Compute x, y, z full time-step forward

		xk4(1)= x(j)+ SpecieT(s)%FluxTubeT(f)%hT(1)*VxR(1)+ &
			((SpecieT(s)%FluxTubeT(f)%hT(1)**2d0)/2d0)*Axk1(1)
		yk4(1)= y(j)+ SpecieT(s)%FluxTubeT(f)%hT(1)*VyR(1)+ &
			((SpecieT(s)%FluxTubeT(f)%hT(1)**2d0)/2d0)*Ayk1(1)
		zk4(1)= z(j)+ SpecieT(s)%FluxTubeT(f)%hT(1)*VzR(1)+ &
			((SpecieT(s)%FluxTubeT(f)%hT(1)**2d0)/2d0)*Azk1(1)

		! ----------------------------------------------------

		call KineticUpdateCSub

		! ----------------------------------------------------

		! INTEGRATE Ax, Ay, Az FOR Vx, Vy, Vz AT FULL TIME-STEP FORWARD:

		k4Vx(1)= Axk4(1) ! Compute k4 values
		k4Vy(1)= Ayk4(1)
		k4Vz(1)= Azk4(1)

		VxNp(1)= VxR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1Vx(1)+ 2d0*k2Vx(1)+ &
			2d0*k3Vx(1)+ k4Vx(1)) ! Compute VxN values
		VyNp(1)= VyR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1Vy(1)+ 2d0*k2Vy(1)+ &
			2d0*k3Vy(1)+ k4Vy(1)) ! Compute VyN values
		VzNp(1)= VzR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1Vz(1)+ 2d0*k2Vz(1)+ &
			2d0*k3Vz(1)+ k4Vz(1)) ! Compute VzN values

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(VxNp(1))) .eqv. .true.) .or. &
			(size(VxNp(:)) /= 1d0)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxNp= ', VxNp(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(VyNp(1))) .eqv. .true.) .or. &
			(size(VyNp(:)) /= 1d0)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyNp= ', VyNp(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(VzNp(1))) .eqv. .true.) .or. &
			(size(VzNp(:)) /= 1d0)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzNp= ', VzNp(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! INTEGRATE Vx, Vy, Vz FOR x, y, z AT FULL TIME-STEP FORWARD:

		k1x(1)= VxR(1) ! Compute k1 values
		k1y(1)= VyR(1)
		k1z(1)= VzR(1)

		Vxk2(1)= VxR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*Axk1(1) ! Compute Vx, Vy, Vz half-step forward
		Vyk2(1)= VyR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*Ayk1(1)
		Vzk2(1)= VzR(1)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/2d0)*Azk1(1)

		k2x(1)= Vxk2(1) ! Compute k2 values
		k2y(1)= Vyk2(1)
		k2z(1)= Vzk2(1)
		k3x(1)= Vxk2(1) ! Compute k3 values (Note: v_i ~= v_i(r_i) for i= x, y, z)
		k3y(1)= Vyk2(1)
		k3z(1)= Vzk2(1)

		Vxk4(1)= VxR(1)+ SpecieT(s)%FluxTubeT(f)%hT(1)*Axk2(1)
		Vyk4(1)= VyR(1)+ SpecieT(s)%FluxTubeT(f)%hT(1)*Ayk2(1)
		Vzk4(1)= VzR(1)+ SpecieT(s)%FluxTubeT(f)%hT(1)*Azk2(1)

		k4x(1)= Vxk4(1) ! Compute k4 values
		k4y(1)= Vyk4(1)
		k4z(1)= Vzk4(1)

		! ----------------------------------------------------

		! UPDATE NEW ION POSITIONS ON CORRECT L-SHELL:

		if (ENAflag(j) .eqv. .false.) then
			xNp(1)= x(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1x(1)+ 2d0*k2x(1)+ 2d0*k3x(1)+ k4x(1))
				! Compute xNp values
			yNp(1)= y(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1y(1)+ 2d0*k2y(1)+ 2d0*k3y(1)+ k4y(1))
				! Compute yNp values
			zNp(1)= z(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1z(1)+ 2d0*k2z(1)+ 2d0*k3z(1)+ k4z(1))
				! Compute zNp values

			! Set displacement dr= ds and update new q position accordingly with reset p (and xN, yN, zN)

			call rsub(rNp(1), xNp(1), yNp(1), zNp(1))
			call thetasub(thetaNp(1), zNp(1), rNp(1))
			call qsub(qNp(1), rNp(1), thetaNp(1))

			pNp(1)= SpecieT(s)%FluxTubeT(f)%pGCT(nnind, 1) ! Reset particle on correct L-shell

			! ----------------------------------------------------

			call RK4DipolePolynomialSolverSub

			! ----------------------------------------------------

	  	xN(j)= xfinalRK4(1)
	  	yN(j)= yfinalRK4(1)
	  	zN(j)= zfinalRK4(1)

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		! UPDATE NEW ENA POSITIONS:

		if (ENAflag(j) .eqv. .true.) then
			xN(j)= x(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1x(1)+ 2d0*k2x(1)+ 2d0*k3x(1)+ k4x(1))
				! Compute xN values
			yN(j)= y(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1y(1)+ 2d0*k2y(1)+ 2d0*k3y(1)+ k4y(1))
				! Compute yN values
			zN(j)= z(j)+ (SpecieT(s)%FluxTubeT(f)%hT(1)/6d0)*(k1z(1)+ 2d0*k2z(1)+ 2d0*k3z(1)+ k4z(1))
				! Compute zN values
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(xN(j))) .eqv. .true.) .or. &
			(size(xN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xN= ', xN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(yN(j))) .eqv. .true.) .or. &
			(size(yN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yN= ', yN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(zN(j))) .eqv. .true.) .or. &
			(size(zN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zN= ', zN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE UPDATED VELOCITY COMPONENTS:

		call rsub(rN(1), xN(j), yN(j), zN(j))
		call thetasub(thetaN(1), zN(j), rN(1))

		if (ENAflag(j) .eqv. .false.) then
			! Note: Without cross L-shell convection, solar forcing or other azimuthal asymmetries, ion motion is phi-invariant.
			phiN(1)= SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1)

		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			call phisub(phiN(1), xN(j), yN(j))
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

		if (ENAflag(j) .eqv. .false.) then
			if (phiN(1) /= SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ION phiN= ', phiN(1), &
					' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCT(nnind, 1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
					' UPDATE SUBROUTINE' // achar(27) // '[0m.'
			end if
			if (pNp(1) /= SpecieT(s)%FluxTubeT(f)%pGCT(nnind, 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ION pNp= ', pNp(1), &
					' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%pGCT(nnind, 1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
					' UPDATE SUBROUTINE' // achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

		call ellsub(ellN(1), thetaN(1))

		! ----------------------------------------------------

		! RE-ALIGN VELOCITY VECTOR INTO LOCAL DIPOLE COORDINATES:

		VparNp(1)= 3d0*cos(thetaN(1))*sin(thetaN(1))*(VxNp(1)* &
			cos(phiN(1))+ VyNp(1)*sin(phiN(1)))/sqrt(ellN(1))+ &
			VzNp(1)*(3d0*(cos(thetaN(1))**2d0)- 1d0)/sqrt(ellN(1))

		! ----------------------------------------------------

		if (ENAflag(j) .eqv. .false.) then
			VpNp(1)= 0d0
			VphiNp(1)= 0d0

			! ----------------------------------------------------

		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			VpNp(1)= abs((1d0- 3d0*(cos(thetaN(1))**2d0))*(VxNp(1)* &
				cos(phiN(1))+ VyNp(1)*sin(phiN(1)))/sqrt(ellN(1))+ 3d0*VzNp(1)* &
				cos(thetaN(1))*sin(thetaN(1))/sqrt(ellN(1)))
			VphiNp(1)= -VxNp(1)*sin(phiN(1))+ VyNp(1)*cos(phiN(1))

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		! COMPUTE NEW FIELD-ALIGNED CARTESIAN VELOCITY COMPONENTS BY COMPUTATIONAL L-SHELL DRIFT LIMIT:

		VxN(j)= cos(phiN(1))*(3d0*VparNp(1)*cos(thetaN(1))*sin(thetaN(1))+ &
			VpNp(1)*(1d0- 3d0*(cos(thetaN(1))**2d0)))/sqrt(ellN(1))- &
			VphiNp(1)*sin(phiN(1))
		VyN(j)= sin(phiN(1))*(3d0*VparNp(1)*cos(thetaN(1))*sin(thetaN(1))+ &
			VpNp(1)*(1d0- 3d0*(cos(thetaN(1))**2d0)))/sqrt(ellN(1))+ &
			VphiNp(1)*cos(phiN(1))
		VzN(j)= VparNp(1)*(3d0*(cos(thetaN(1))**2d0)- 1d0)/sqrt(ellN(1))+ &
			3d0*VpNp(1)*cos(thetaN(1))*sin(thetaN(1))/sqrt(ellN(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(VxN(j))) .eqv. .true.) .or. &
			(size(VxN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxN= ', VxN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(VyN(j))) .eqv. .true.) .or. &
			(size(VyN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyN= ', VyN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((isnan(real(VzN(j))) .eqv. .true.) .or. &
			(size(VzN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzN= ', VzN(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC RK4', &
				' UPDATE SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if ((n == 1) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S81End)
				write(S81string, '(i10)')  nint(S81End)
				write(*, *) trim('%% 8.11- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S81string)) // &
					trim(' s. INITIAL KINETIC RK4 COMPLETE %%')
			end if
		end if

		if (rank == 0) then
			if ((n == 2) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S81End)
				write(S81string, '(i10)')  nint(S81End)
				write(*, *) trim('%% 8.11- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S81string)) // &
					trim(' s. SECOND KINETIC RK4 COMPLETE %%')
			end if
		end if

		if (rank == 0) then
			if ((n == SpecieT(s)%FluxTubeT(f)%ndatfacT(1)- 1) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S81End)
				write(S81string, '(i10)')  nint(S81End)
				write(*, *) trim('%% 8.11- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S81string)) // &
					trim(' s. PENULTIMATE TO SECOND STATISTICAL TIME-STEP KINETIC RK4 COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticRK4UpdateSub

end module KineticRK4Update
