module KineticUpdateC

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.3 KINETIC UPDATE C:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! UPDATE CARTESIAN KINETIC VARIABLES BY RK4 INTEGRATION:

	subroutine KineticUpdateCSub

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT FULL TIME-STEP FORWARD:

		call rsub(rk4(j), xk4(1), yk4(1), zk4(1))
		call thetasub(thetak4(j), zk4(1), rk4(j))

		if (ENAflag(j) .eqv. .false.) then
			! Note: Without cross L-shell convection, solar forcing or other azimuthal asymmetries, ion motion is phi-invariant.
			phik4(j)= SpecieT(s)%FluxTubeT(f)%phiGCGT(nnind, 1)
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			call phisub(phik4(j), xk4(1), yk4(1))
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

		if (ENAflag(j) .eqv. .false.) then
			if (phik4(j) /= SpecieT(s)%FluxTubeT(f)%phiGCGT(nnind, 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ION phik4= ', phik4(j), &
					' AND phiGCGT VALUE= ', SpecieT(s)%FluxTubeT(f)%phiGCGT(nnind, 1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
					// achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

		call qsub(qk4(j), rk4(j), thetak4(j))
		pk4(j)= SpecieT(s)%FluxTubeT(f)%pGCGT(nnind, 1)
		!call psub(pk4(j), rk4(j), thetak4(j))
		call ellsub(ellk4(j), thetak4(j))
		call Bmagsub(Bmagk4(1), rk4(j), ellk4(j))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(rk4(j))) .eqv. .true.) .or. &
			(size(rk4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rk4= ', rk4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(thetak4(j))) .eqv. .true.) .or. &
			(size(thetak4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetak4= ', thetak4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(phik4(j))) .eqv. .true.) .or. &
			(size(phik4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phik4= ', phik4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(qk4(j))) .eqv. .true.) .or. &
			(size(qk4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qk4= ', qk4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(pk4(j))) .eqv. .true.) .or. &
			(size(pk4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pk4= ', pk4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(ellk4(j))) .eqv. .true.) .or. &
			(size(ellk4(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellk4= ', ellk4(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Bmagk4(1))) .eqv. .true.) .or. &
			(size(Bmagk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Bmagk4= ', Bmagk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE ION DRIFTS OVER ENTIRE SIMULATION:

		!if (j < SpecieT(s)%FluxTubeT(f)%NsnT(1)) then
		!	do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		!		if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
		!			(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

		!			pdriftion(j)= abs(pk4(j)- SpecieT(s)%FluxTubeT(f)%p0T(1))

		!		end if
		!	end do

		!	if ((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. (ENAflag(j) .eqv. .false.)) then

		!		qdriftion(j)= abs(qk4(j)- SpecieT(s)%FluxTubeT(f)%q0T(j))
		!		phidriftion(j)= abs(phik4(j)- SpecieT(s)%FluxTubeT(f)%phi0T(j))

		!	end if
		!end if

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT LAST SPLIT TIME-STEP:

		call dBdssub(dBdsk4(1), rk4(j), thetak4(j), ellk4(j))

		! ----------------------------------------------------

		if (ENAflag(j) .eqv. .false.) then
			call Rperpsub(Rperpk4(1), pk4(j), phik4(j))
			call OmegaGsub(OmegaGk4(1), SpecieT(s)%qsT(1), SpecieT(s)%msT(1), Bmagk4(1))
			call musub(muk4(1), SpecieT(s)%msT(1), Bmagk4(1), Vperp(j))
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			Rperpk4(1)= 0d0
			OmegaGk4(1)= 0d0
			muk4(1)= 0d0
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(dBdsk4(1))) .eqv. .true.) .or. &
			(size(dBdsk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dBdsk4= ', dBdsk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Rperpk4(1))) .eqv. .true.) .or. &
			(size(Rperpk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Rperpk4= ', Rperpk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(OmegaGk4(1))) .eqv. .true.) .or. &
			(size(OmegaGk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaGk4= ', OmegaGk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(muk4(1))) .eqv. .true.) .or. &
			(size(muk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' muk4= ', muk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN ACCELERATION COMPONENTS AT FULL TIME-STEP FORWARD:

		if ((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.)) then

			call AMxsub(AMxk4p(1), SpecieT(s)%msT(1), muk4(1), &
				dBdsk4(1), thetak4(j), phik4(j), ellk4(j))
			call AMysub(AMyk4p(1), SpecieT(s)%msT(1), muk4(1), &
				dBdsk4(1), thetak4(j), phik4(j), ellk4(j))
			call AMzsub(AMzk4p(1), SpecieT(s)%msT(1), muk4(1), &
				dBdsk4(1), thetak4(j), ellk4(j))

			if (qk4(j) <= 0d0) then  ! S Mag Hemisphere
				AMpark4(1)= -abs(3d0*cos(thetak4(j))*sin(thetak4(j))*(AMxk4p(1)* &
					cos(phik4(j))+ AMyk4p(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
					AMzk4p(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j)))
			else if (qk4(j) > 0d0) then ! N Mag Hemisphere
				AMpark4(1)= abs(3d0*cos(thetak4(j))*sin(thetak4(j))*(AMxk4p(1)* &
					cos(phik4(j))+ AMyk4p(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
					AMzk4p(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j)))
			end if

			AMpk4(1)= 0d0
			AMphik4(1)= 0d0

			AMxk4(1)= cos(phik4(j))*(3d0*AMpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AMpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))- &
				AMphik4(1)*sin(phik4(j))
			AMyk4(1)= sin(phik4(j))*(3d0*AMpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AMpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))+ &
				AMphik4(1)*cos(phik4(j))
			AMzk4(1)= AMpark4(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))+ &
				3d0*AMpk4(1)*cos(thetak4(j))*sin(thetak4(j))/sqrt(ellk4(j))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk4(j) <= 0) .and. (AMpark4(1) > 0d0)) .or. &
				((qk4(j) > 0) .and. (AMpark4(1) < 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AMpark4= ', AMpark4(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE C SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AMxk4p(1)= 0d0
			AMyk4p(1)= 0d0
			AMzk4p(1)= 0d0
			AMpark4(1)= 0d0
			AMpk4(1)= 0d0
			AMphik4(1)= 0d0
			AMxk4(1)= 0d0
			AMyk4(1)= 0d0
			AMzk4(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AMxk4p(1))) .eqv. .true.) .or. &
			(size(AMxk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk4p= ', AMxk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk4p(1))) .eqv. .true.) .or. &
			(size(AMyk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk4p= ', AMyk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk4p(1))) .eqv. .true.) .or. &
			(size(AMzk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk4p= ', AMzk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMphik4(1))) .eqv. .true.) .or. &
			(size(AMphik4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMphik4= ', AMphik4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpark4(1))) .eqv. .true.) .or. &
			(size(AMpark4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpark4= ', AMpark4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpk4(1))) .eqv. .true.) .or. &
			(size(AMpk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpk4= ', AMpk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMxk4(1))) .eqv. .true.) .or. &
			(size(AMxk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk4= ', AMxk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk4(1))) .eqv. .true.) .or. &
			(size(AMyk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk4= ', AMyk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk4(1))) .eqv. .true.) .or. &
			(size(AMzk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk4= ', AMzk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN GRAVITATIONAL FORCE ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.) &
			.and. (Qindk1(j) /= 0d0) .and. (Qindk1(j) /= -1d0)) then

			if (qk4(j) <= 0d0) then ! SMH
				AGmagSk4(1)= 1d0
			else if (qk4(j) > 0d0) then ! NMH
				AGmagSk4(1)= -1d0
			end if

			call AGxsub(AGxk4p(1), AGmagSk4(1), AGmag(j), thetak4(j), phik4(j), ellk4(j))
			call AGysub(AGyk4p(1), AGmagSk4(1), AGmag(j), thetak4(j), phik4(j), ellk4(j))
			call AGzsub(AGzk4p(1), AGmagSk4(1), AGmag(j), thetak4(j), ellk4(j))

			AGpark4(1)= (3d0*cos(thetak4(j))*sin(thetak4(j))*(AGxk4p(1)* &
				cos(phik4(j))+ AGyk4p(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
				AGzk4p(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j)))

			AGpk4(1)= 0d0
			AGphik4(1)= 0d0

			AGxk4(1)= cos(phik4(j))*(3d0*AGpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AGpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))- &
				AGphik4(1)*sin(phik4(j))
			AGyk4(1)= sin(phik4(j))*(3d0*AGpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AGpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))+ &
				AGphik4(1)*cos(phik4(j))
			AGzk4(1)= AGpark4(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))+ &
				3d0*AGpk4(1)*cos(thetak4(j))*sin(thetak4(j))/sqrt(ellk4(j))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk4(j) <= 0) .and. (AGpark4(1) < 0d0)) .or. &
				((qk4(j) > 0) .and. (AGpark4(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark4= ', AGpark4(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE C SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then

			call AGENAxsub(AGxk4(1), mNeut, rk4(j), &
				thetak4(j), phik4(j))
			call AGENAysub(AGyk4(1), mNeut, rk4(j), &
				thetak4(j), phik4(j))
			call AGENAzsub(AGzk4(1), mNeut, rk4(j), &
				thetak4(j))

			AGpark4(1)= 3d0*cos(thetak4(j))*sin(thetak4(j))*(AGxk4(1)* &
				cos(phik4(j))+ AGyk4(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
				AGzk4(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk4(j) <= 0) .and. (AGpark4(1) < 0d0)) .or. &
				((qk4(j) > 0) .and. (AGpark4(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark4= ', AGpark4(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE C SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AGxk4p(1)= 0d0
			AGyk4p(1)= 0d0
			AGzk4p(1)= 0d0
			AGpark4(1)= 0d0
			AGpk4(1)= 0d0
			AGphik4(1)= 0d0
			AGxk4(1)= 0d0
			AGyk4(1)= 0d0
			AGzk4(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AGxk4p(1))) .eqv. .true.) .or. &
			(size(AGxk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk4p= ', AGxk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk4p(1))) .eqv. .true.) .or. &
			(size(AGyk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk4p= ', AGyk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk4p(1))) .eqv. .true.) .or. &
			(size(AGzk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk4p= ', AGzk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGphik4(1))) .eqv. .true.) .or. &
			(size(AGphik4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGphik4= ', AGphik4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpark4(1))) .eqv. .true.) .or. &
			(size(AGpark4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpark4= ', AGpark4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpk4(1))) .eqv. .true.) .or. &
			(size(AGpk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpk4= ', AGpk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGxk4(1))) .eqv. .true.) .or. &
			(size(AGxk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk4= ', AGxk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk4(1))) .eqv. .true.) .or. &
			(size(AGyk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk4= ', AGyk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk4(1))) .eqv. .true.) .or. &
			(size(AGzk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk4= ', AGzk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN AMBIPOLAR ELECTRIC FIELD ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (qk4(j) <= 0d0) then ! SMH
					AEAmagSk4(1)= -1d0
				else if (qk4(j) > 0d0) then ! NMH
					AEAmagSk4(1)= 1d0
				end if
			end if
			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
				AEAmagSk4(1)= 1d0
			end if

			call AEAxsub(AEAxk4p(1), AEAmagSk4(1), AEAmag(j), thetak4(j), phik4(j), ellk4(j))
			call AEAysub(AEAyk4p(1), AEAmagSk4(1), AEAmag(j), thetak4(j), phik4(j), ellk4(j))
			call AEAzsub(AEAzk4p(1), AEAmagSk4(1), AEAmag(j), thetak4(j), ellk4(j))

			AEApark4(1)= (3d0*cos(thetak4(j))*sin(thetak4(j))*(AEAxk4p(1)* &
				cos(phik4(j))+ AEAyk4p(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
				AEAzk4p(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j)))

			AEApk4(1)= 0d0
			AEAphik4(1)= 0d0

			AEAxk4(1)= cos(phik4(j))*(3d0*AEApark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AEApk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))- &
				AEAphik4(1)*sin(phik4(j))
			AEAyk4(1)= sin(phik4(j))*(3d0*AEApark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AEApk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))+ &
				AEAphik4(1)*cos(phik4(j))
			AEAzk4(1)= AEApark4(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))+ &
				3d0*AEApk4(1)*cos(thetak4(j))*sin(thetak4(j))/sqrt(ellk4(j))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER AMBIPOLAR ACCELERATION SIGN:

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (((qk4(j) <= 0) .and. (AEApark4(1) > 0d0)) .or. &
					((qk4(j) > 0) .and. (AEApark4(1) < 0d0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
						' AEApark4= ', AEApark4(1), 'qk4= ', qk4(j), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
						' IN KINETIC UPDATE C SUBROUTINE' // achar(27) // '[0m.'
				end if
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEAxk4p(1)= 0d0
			AEAyk4p(1)= 0d0
			AEAzk4p(1)= 0d0
			AEApark4(1)= 0d0
			AEApk4(1)= 0d0
			AEAphik4(1)= 0d0
			AEAxk4(1)= 0d0
			AEAyk4(1)= 0d0
			AEAzk4(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEAxk4p(1))) .eqv. .true.) .or. &
			(size(AEAxk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk4p', AEAxk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk4p(1))) .eqv. .true.) .or. &
			(size(AEAyk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk4p', AEAyk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk4p(1))) .eqv. .true.) .or. &
			(size(AEAzk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk4p', AEAzk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAphik4(1))) .eqv. .true.) .or. &
			(size(AEAphik4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAphik4= ', AEAphik4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApark4(1))) .eqv. .true.) .or. &
			(size(AEApark4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApark4= ', AEApark4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApk4(1))) .eqv. .true.) .or. &
			(size(AEApk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApk4= ', AEApk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAxk4(1))) .eqv. .true.) .or. &
			(size(AEAxk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk4', AEAxk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk4(1))) .eqv. .true.) .or. &
			(size(AEAyk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk4', AEAyk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk4(1))) .eqv. .true.) .or. &
			(size(AEAzk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk4', AEAzk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN PARALLEL ELECTRIC FIELD ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			if (qk4(j) <= 0d0) then ! SMH
				AEPmagSk4(1)= 1d0
			else if (qk4(j) > 0d0) then ! NMH
				AEPmagSk4(1)= -1d0
			end if

			if (n < EParlim) then
				AEPmag(j)= 0d0
			end if
			if (n >= EParlim) then
				AEPmag(j)= AEPmag(j)
			end if

			call AEPxsub(AEPxk4p(1), AEPmagSk4(1), AEPmag(j), thetak4(j), phik4(j), ellk4(j))
			call AEPysub(AEPyk4p(1), AEPmagSk4(1), AEPmag(j), thetak4(j), phik4(j), ellk4(j))
			call AEPzsub(AEPzk4p(1), AEPmagSk4(1), AEPmag(j), thetak4(j), ellk4(j))

			AEPpark4(1)= (3d0*cos(thetak4(j))*sin(thetak4(j))*(AEPxk4p(1)* &
				cos(phik4(j))+ AEPyk4p(1)*sin(phik4(j)))/sqrt(ellk4(j))+ &
				AEPzk4p(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j)))

			AEPpk4(1)= 0d0
			AEPphik4(1)= 0d0

			AEPxk4(1)= cos(phik4(j))*(3d0*AEPpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AEPpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))- &
				AEPphik4(1)*sin(phik4(j))
			AEPyk4(1)= sin(phik4(j))*(3d0*AEPpark4(1)*cos(thetak4(j))*sin(thetak4(j))+ &
				AEPpk4(1)*(1d0- 3d0*(cos(thetak4(j))**2d0)))/sqrt(ellk4(j))+ &
				AEPphik4(1)*cos(phik4(j))
			AEPzk4(1)= AEPpark4(1)*(3d0*(cos(thetak4(j))**2d0)- 1d0)/sqrt(ellk4(j))+ &
				3d0*AEPpk4(1)*cos(thetak4(j))*sin(thetak4(j))/sqrt(ellk4(j))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk4(j) <= 0) .and. (AEPpark4(1) < 0d0)) .or. &
				((qk4(j) > 0) .and. (AEPpark4(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AEPpark4= ', AEPpark4(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE C SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEPxk4p(1)= 0d0
			AEPyk4p(1)= 0d0
			AEPzk4p(1)= 0d0
			AEPpark4(1)= 0d0
			AEPpk4(1)= 0d0
			AEPphik4(1)= 0d0
			AEPxk4(1)= 0d0
			AEPyk4(1)= 0d0
			AEPzk4(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEPxk4p(1))) .eqv. .true.) .or. &
			(size(AEPxk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk4p', AEPxk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk4p(1))) .eqv. .true.) .or. &
			(size(AEPyk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk4p', AEPyk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk4p(1))) .eqv. .true.) .or. &
			(size(AEPzk4p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk4p', AEPzk4p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPphik4(1))) .eqv. .true.) .or. &
			(size(AEPphik4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPphik4= ', AEPphik4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpark4(1))) .eqv. .true.) .or. &
			(size(AEPpark4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpark4= ', AEPpark4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpk4(1))) .eqv. .true.) .or. &
			(size(AEPpk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpk4= ', AEPpk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPxk4(1))) .eqv. .true.) .or. &
			(size(AEPxk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk4', AEPxk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk4(1))) .eqv. .true.) .or. &
			(size(AEPyk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk4', AEPyk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk4(1))) .eqv. .true.) .or. &
			(size(AEPzk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk4', AEPzk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN NET FORCE ACCELERATION COMPONENTS:

		call Axsub(Axk4(1), AMxk4(1), AGxk4(1), AEAxk4(1), AEPxk4(1))
		call Aysub(Ayk4(1), AMyk4(1), AGyk4(1), AEAyk4(1), AEPyk4(1))
		call Azsub(Azk4(1), AMzk4(1), AGzk4(1), AEAzk4(1), AEPzk4(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(Axk4(1))) .eqv. .true.) .or. &
			(size(Axk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Axk4= ', Axk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Ayk4(1))) .eqv. .true.) .or. &
			(size(Ayk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Ayk4= ', Ayk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Azk4(1))) .eqv. .true.) .or. &
			(size(Azk4(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Azk4= ', Azk4(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE C SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if ((n == 1) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S813End)
				write(S813string, '(F10.4)')  S813End
				write(*, *) trim('%% 8.11.3- REAL-TIME= ' // adjustl(S813string)) // &
					trim(' s. INITIAL KINETIC UPDATE C COMPLETE %%')
			end if
		end if

		if (rank == 0) then
			if ((n == 2) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S813End)
				write(S813string, '(F10.4)')  S813End
				write(*, *) trim('%% 8.11.3- REAL-TIME= ' // adjustl(S813string)) // &
					trim(' s. SECOND KINETIC UPDATE C COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticUpdateCSub

end module KineticUpdateC
