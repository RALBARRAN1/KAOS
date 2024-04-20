module KineticUpdateB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.2 KINETIC UPDATE B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! UPDATE CARTESIAN KINETIC VARIABLES BY RK4 INTEGRATION:

	subroutine KineticUpdateBSub

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT HALF TIME-STEP FORWARD:

		call rsub(rk2(1), xk2(1), yk2(1), zk2(1))
		call thetasub(thetak2(1), zk2(1), rk2(1))

		if (ENAflag(j) .eqv. .false.) then
			! Note: Without cross L-shell convection, solar forcing or other azimuthal asymmetries, ion motion is phi-invariant.
			phik2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1)
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			call phisub(phik2(1), xk2(1), yk2(1))
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

		if (ENAflag(j) .eqv. .false.) then
			if (phik2(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ION phik2= ', phik2(1), &
					' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
					// achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

		call qsub(qk2(1), rk2(1), thetak2(1))
		call psub(pk2(1), rk2(1), thetak2(1))
		call ellsub(ellk2(1), thetak2(1))
		call Bmagsub(Bmagk2(1), rk2(1), ellk2(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(rk2(1))) .eqv. .true.) .or. &
			(size(rk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rk2= ', rk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(thetak2(1))) .eqv. .true.) .or. &
			(size(thetak2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetak2= ', thetak2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(phik2(1))) .eqv. .true.) .or. &
			(size(phik2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phik2= ', phik2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(qk2(1))) .eqv. .true.) .or. &
			(size(qk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qk2= ', qk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(pk2(1))) .eqv. .true.) .or. &
			(size(pk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pk2= ', pk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(ellk2(1))) .eqv. .true.) .or. &
			(size(ellk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellk2= ', ellk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Bmagk2(1))) .eqv. .true.) .or. &
			(size(Bmagk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Bmagk2= ', Bmagk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT SECOND SPLIT TIME-STEP:

		call dBdssub(dBdsk2(1), rk2(1), thetak2(1), ellk2(1))

		! ----------------------------------------------------

		if (ENAflag(j) .eqv. .false.) then
			call Rperpsub(Rperpk2(1), pk2(1), phik2(1))
			call OmegaGsub(OmegaGk2(1), SpecieT(s)%qsT(1), SpecieT(s)%msT(1), Bmagk2(1))
			call musub(muk2(1), SpecieT(s)%msT(1), Bmagk2(1), Vperp(j))
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			Rperpk2(1)= 0d0
			OmegaGk2(1)= 0d0
			muk2(1)= 0d0
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(dBdsk2(1))) .eqv. .true.) .or. &
			(size(dBdsk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dBdsk2= ', dBdsk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Rperpk2(1))) .eqv. .true.) .or. &
			(size(Rperpk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Rperpk2= ', Rperpk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(OmegaGk2(1))) .eqv. .true.) .or. &
			(size(OmegaGk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaGk2= ', OmegaGk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(muk2(1))) .eqv. .true.) .or. &
			(size(muk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' muk2= ', muk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN MIRROR FORCE ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.)) then

			call AMxsub(AMxk2p(1), SpecieT(s)%msT(1), muk2(1), &
				dBdsk2(1), thetak2(1), phik2(1), ellk2(1))
			call AMysub(AMyk2p(1), SpecieT(s)%msT(1), muk2(1), &
				dBdsk2(1), thetak2(1), phik2(1), ellk2(1))
			call AMzsub(AMzk2p(1), SpecieT(s)%msT(1), muk2(1), &
				dBdsk2(1), thetak2(1), ellk2(1))

			if (qk2(1) <= 0d0) then  ! N Mag Hemisphere
				AMpark2(1)= abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AMxk2p(1)* &
					cos(phik2(1))+ AMyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AMzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			else if (qk2(1) > 0d0) then ! S Mag Hemisphere
				AMpark2(1)= -abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AMxk2p(1)* &
					cos(phik2(1))+ AMyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AMzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			end if

			AMpk2(1)= 0d0
			AMphik2(1)= 0d0

			AMxk2(1)= cos(phik2(1))*(3d0*AMpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AMpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))- &
				AMphik2(1)*sin(phik2(1))
			AMyk2(1)= sin(phik2(1))*(3d0*AMpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AMpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))+ &
				AMphik2(1)*cos(phik2(1))
			AMzk2(1)= AMpark2(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))+ &
				3d0*AMpk2(1)*cos(thetak2(1))*sin(thetak2(1))/sqrt(ellk2(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk2(1) <= 0) .and. (AMpark2(1) < 0d0)) .or. &
				((qk2(1) > 0) .and. (AMpark2(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AMpark2= ', AMpark2(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE B SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AMxk2p(1)= 0d0
			AMyk2p(1)= 0d0
			AMzk2p(1)= 0d0
			AMpark2(1)= 0d0
			AMpk2(1)= 0d0
			AMphik2(1)= 0d0
			AMxk2(1)= 0d0
			AMyk2(1)= 0d0
			AMzk2(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AMxk2p(1))) .eqv. .true.) .or. &
			(size(AMxk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk2p= ', AMxk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk2p(1))) .eqv. .true.) .or. &
			(size(AMyk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk2p= ', AMyk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk2p(1))) .eqv. .true.) .or. &
			(size(AMzk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk2p= ', AMzk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMphik2(1))) .eqv. .true.) .or. &
			(size(AMphik2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMphik2= ', AMphik2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpark2(1))) .eqv. .true.) .or. &
			(size(AMpark2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpark2= ', AMpark2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpk2(1))) .eqv. .true.) .or. &
			(size(AMpk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpk2= ', AMpk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMxk2(1))) .eqv. .true.) .or. &
			(size(AMxk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk2= ', AMxk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk2(1))) .eqv. .true.) .or. &
			(size(AMyk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk2= ', AMyk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk2(1))) .eqv. .true.) .or. &
			(size(AMzk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk2= ', AMzk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN GRAVITATIONAL FORCE ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.)) then

			call AGxsub(AGxk2p(1), AGmag(j), thetak2(1), phik2(1), ellk2(1))
			call AGysub(AGyk2p(1), AGmag(j), thetak2(1), phik2(1), ellk2(1))
			call AGzsub(AGzk2p(1), AGmag(j), thetak2(1), ellk2(1))

			if (qk2(1) <= 0d0) then  ! N Mag Hemisphere
				AGpark2(1)= -abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AGxk2p(1)* &
					cos(phik2(1))+ AGyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AGzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			end if
			if (qk2(1) > 0d0) then  ! S Mag Hemisphere
				AGpark2(1)= abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AGxk2p(1)* &
					cos(phik2(1))+ AGyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AGzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			end if

			AGpk2(1)= 0d0
			AGphik2(1)= 0d0

			AGxk2(1)= cos(phik2(1))*(3d0*AGpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AGpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))- &
				AGphik2(1)*sin(phik2(1))
			AGyk2(1)= sin(phik2(1))*(3d0*AGpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AGpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))+ &
				AGphik2(1)*cos(phik2(1))
			AGzk2(1)= AGpark2(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))+ &
				3d0*AGpk2(1)*cos(thetak2(1))*sin(thetak2(1))/sqrt(ellk2(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk2(1) <= 0) .and. (AGpark2(1) > 0d0)) .or. &
				((qk2(1) > 0) .and. (AGpark2(1) < 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark2= ', AGpark2(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE B SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

			else if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. &
				(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
				(ENAflag(j) .eqv. .true.)) then

			call AGENAxsub(AGxk2(1), mNeut, rk2(1), &
				thetak2(1), phik2(1))
			call AGENAysub(AGyk2(1), mNeut, rk2(1), &
				thetak2(1), phik2(1))
			call AGENAzsub(AGzk2(1), mNeut, rk2(1), &
				thetak2(1))

			AGpark2(1)= 3d0*cos(thetak2(1))*sin(thetak2(1))*(AGxk2(1)* &
				cos(phik2(1))+ AGyk2(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
				AGzk2(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk2(1) <= 0) .and. (AGpark2(1) > 0d0)) .or. &
				((qk2(1) > 0) .and. (AGpark2(1) < 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark2= ', AGpark2(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE B SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AGxk2p(1)= 0d0
			AGyk2p(1)= 0d0
			AGzk2p(1)= 0d0
			AGpark2(1)= 0d0
			AGpk2(1)= 0d0
			AGphik2(1)= 0d0
			AGxk2(1)= 0d0
			AGyk2(1)= 0d0
			AGzk2(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AGxk2p(1))) .eqv. .true.) .or. &
			(size(AGxk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk2p= ', AGxk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk2p(1))) .eqv. .true.) .or. &
			(size(AGyk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk2p= ', AGyk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk2p(1))) .eqv. .true.) .or. &
			(size(AGzk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk2p= ', AGzk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGphik2(1))) .eqv. .true.) .or. &
			(size(AGphik2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGphik2= ', AGphik2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpark2(1))) .eqv. .true.) .or. &
			(size(AGpark2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpark2= ', AGpark2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpk2(1))) .eqv. .true.) .or. &
			(size(AGpk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpk2= ', AGpk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGxk2(1))) .eqv. .true.) .or. &
			(size(AGxk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk2= ', AGxk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk2(1))) .eqv. .true.) .or. &
			(size(AGyk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk2= ', AGyk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk2(1))) .eqv. .true.) .or. &
			(size(AGzk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk2= ', AGzk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN AMBIPOLAR ELECTRIC FIELD ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			call AEAxsub(AEAxk2p(1), AEAmag(j), thetak2(1), phik2(1), ellk2(1))
			call AEAysub(AEAyk2p(1), AEAmag(j), thetak2(1), phik2(1), ellk2(1))
			call AEAzsub(AEAzk2p(1), AEAmag(j), thetak2(1), ellk2(1))

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (qk2(1) <= 0d0) then  ! N Mag Hemisphere
					AEApark2(1)= abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AEAxk2p(1)* &
						cos(phik2(1))+ AEAyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
						AEAzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
				end if
				if (qk2(1) > 0d0) then  ! S Mag Hemisphere
					AEApark2(1)= -abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AEAxk2p(1)* &
						cos(phik2(1))+ AEAyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
						AEAzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
				end if
			end if

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
				AEApark2(1)= 3d0*cos(thetak2(1))*sin(thetak2(1))*(AEAxk2p(1)* &
					cos(phik2(1))+ AEAyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AEAzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))
			end if

			AEApk2(1)= 0d0
			AEAphik2(1)= 0d0

			AEAxk2(1)= cos(phik2(1))*(3d0*AEApark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AEApk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))- &
				AEAphik2(1)*sin(phik2(1))
			AEAyk2(1)= sin(phik2(1))*(3d0*AEApark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AEApk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))+ &
				AEAphik2(1)*cos(phik2(1))
			AEAzk2(1)= AEApark2(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))+ &
				3d0*AEApk2(1)*cos(thetak2(1))*sin(thetak2(1))/sqrt(ellk2(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER AMBIPOLAR ACCELERATION SIGN:

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (((qk2(1) <= 0) .and. (AEApark2(1) < 0d0)) .or. &
					((qk2(1) > 0) .and. (AEApark2(1) > 0d0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
						' AEApark2= ', AEApark2(1), 'qk2= ', qk2(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
						' IN KINETIC UPDATE B SUBROUTINE' // achar(27) // '[0m.'
				end if
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEAxk2p(1)= 0d0
			AEAyk2p(1)= 0d0
			AEAzk2p(1)= 0d0
			AEApark2(1)= 0d0
			AEApk2(1)= 0d0
			AEAphik2(1)= 0d0
			AEAxk2(1)= 0d0
			AEAyk2(1)= 0d0
			AEAzk2(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEAxk2p(1))) .eqv. .true.) .or. &
			(size(AEAxk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk2p', AEAxk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk2p(1))) .eqv. .true.) .or. &
			(size(AEAyk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk2p', AEAyk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk2p(1))) .eqv. .true.) .or. &
			(size(AEAzk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk2p', AEAzk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAphik2(1))) .eqv. .true.) .or. &
			(size(AEAphik2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAphik2= ', AEAphik2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApark2(1))) .eqv. .true.) .or. &
			(size(AEApark2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApark2= ', AEApark2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApk2(1))) .eqv. .true.) .or. &
			(size(AEApk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApk2= ', AEApk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAxk2(1))) .eqv. .true.) .or. &
			(size(AEAxk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk2', AEAxk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk2(1))) .eqv. .true.) .or. &
			(size(AEAyk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk2', AEAyk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk2(1))) .eqv. .true.) .or. &
			(size(AEAzk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk2', AEAzk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN PARALLEL ELECTRIC FIELD ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			if (n < EParlim) then
				AEPmag(j)= 0d0
			end if
			if (n >= EParlim) then
				AEPmag(j)= AEPmag(j)
			end if

			call AEPxsub(AEPxk2p(1), AEPmag(j), thetak2(1), phik2(1), ellk2(1))
			call AEPysub(AEPyk2p(1), AEPmag(j), thetak2(1), phik2(1), ellk2(1))
			call AEPzsub(AEPzk2p(1), AEPmag(j), thetak2(1), ellk2(1))

			if (qk2(1) <= 0d0) then  ! N Mag Hemisphere
				AEPpark2(1)= -abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AEPxk2p(1)* &
					cos(phik2(1))+ AEPyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AEPzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			end if
			if (qk2(1) > 0d0) then  ! S Mag Hemisphere
				AEPpark2(1)= abs(3d0*cos(thetak2(1))*sin(thetak2(1))*(AEPxk2p(1)* &
					cos(phik2(1))+ AEPyk2p(1)*sin(phik2(1)))/sqrt(ellk2(1))+ &
					AEPzk2p(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1)))
			end if

			AEPpk2(1)= 0d0
			AEPphik2(1)= 0d0

			AEPxk2(1)= cos(phik2(1))*(3d0*AEPpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AEPpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))- &
				AEPphik2(1)*sin(phik2(1))
			AEPyk2(1)= sin(phik2(1))*(3d0*AEPpark2(1)*cos(thetak2(1))*sin(thetak2(1))+ &
				AEPpk2(1)*(1d0- 3d0*(cos(thetak2(1))**2d0)))/sqrt(ellk2(1))+ &
				AEPphik2(1)*cos(phik2(1))
			AEPzk2(1)= AEPpark2(1)*(3d0*(cos(thetak2(1))**2d0)- 1d0)/sqrt(ellk2(1))+ &
				3d0*AEPpk2(1)*cos(thetak2(1))*sin(thetak2(1))/sqrt(ellk2(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk2(1) <= 0) .and. (AEPpark2(1) > 0d0)) .or. &
				((qk2(1) > 0) .and. (AEPpark2(1) < 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AEPpark2= ', AEPpark2(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE B SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEPxk2p(1)= 0d0
			AEPyk2p(1)= 0d0
			AEPzk2p(1)= 0d0
			AEPpark2(1)= 0d0
			AEPpk2(1)= 0d0
			AEPphik2(1)= 0d0
			AEPxk2(1)= 0d0
			AEPyk2(1)= 0d0
			AEPzk2(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEPxk2p(1))) .eqv. .true.) .or. &
			(size(AEPxk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk2p', AEPxk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk2p(1))) .eqv. .true.) .or. &
			(size(AEPyk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk2p', AEPyk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk2p(1))) .eqv. .true.) .or. &
			(size(AEPzk2p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk2p', AEPzk2p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPphik2(1))) .eqv. .true.) .or. &
			(size(AEPphik2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPphik2= ', AEPphik2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpark2(1))) .eqv. .true.) .or. &
			(size(AEPpark2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpark2= ', AEPpark2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpk2(1))) .eqv. .true.) .or. &
			(size(AEPpk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpk2= ', AEPpk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPxk2(1))) .eqv. .true.) .or. &
			(size(AEPxk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk2', AEPxk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk2(1))) .eqv. .true.) .or. &
			(size(AEPyk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk2', AEPyk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk2(1))) .eqv. .true.) .or. &
			(size(AEPzk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk2', AEPzk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN NET FORCE ACCELERATION COMPONENTS:

		call Axsub(Axk2(1), AMxk2(1), AGxk2(1), AEAxk2(1), AEPxk2(1))
		call Aysub(Ayk2(1), AMyk2(1), AGyk2(1), AEAyk2(1), AEPyk2(1))
		call Azsub(Azk2(1), AMzk2(1), AGzk2(1), AEAzk2(1), AEPzk2(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(Axk2(1))) .eqv. .true.) .or. &
			(size(Axk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Axk2= ', Axk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Ayk2(1))) .eqv. .true.) .or. &
			(size(Ayk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Ayk2= ', Ayk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Azk2(1))) .eqv. .true.) .or. &
			(size(Azk2(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Azk2= ', Azk2(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE B SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if ((n == 1) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S812End)
				write(S812string, '(i10)')  nint(S812End)
				write(*, *) trim('%% 8.11.2- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S812string)) // &
					trim(' s. INITIAL KINETIC UPDATE B COMPLETE %%')
			end if
		end if

		if (rank == 0) then
			if ((n == 2) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S812End)
				write(S812string, '(i10)')  nint(S812End)
				write(*, *) trim('%% 8.11.2- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S812string)) // &
					trim(' s. SECOND KINETIC UPDATE B COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticUpdateBSub

end module KineticUpdateB
