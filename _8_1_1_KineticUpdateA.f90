module KineticUpdateA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.1 KINETIC UPDATE A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use WaveHeating

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! UPDATE CARTESIAN KINETIC VARIABLES BY RK4 INTEGRATION:

	subroutine KineticUpdateASub

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(x(j))) .eqv. .true.) .or. &
			(size(x(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' x= ', x(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(y(j))) .eqv. .true.) .or. &
			(size(y(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' y= ', y(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(z(j))) .eqv. .true.) .or. &
			(size(z(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' z= ', z(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vperp1(j))) .eqv. .true.) .or. &
			(size(Vperp1(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vperp1= ', Vperp1(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vperp2(j))) .eqv. .true.) .or. &
			(size(Vperp2(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vperp2= ', Vperp2(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vperp(j))) .eqv. .true.) .or. &
			(size(Vperp(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vperp= ', Vperp(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vx(j))) .eqv. .true.) .or. &
			(size(Vx(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vx= ', Vx(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vy(j))) .eqv. .true.) .or. &
			(size(Vy(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vy= ', Vy(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vz(j))) .eqv. .true.) .or. &
			(size(Vz(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' Vz= ', Vz(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR', &
				' SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
				', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT INITIAL TIME-STEP:

		call rsub(rk1(1), x(j), y(j), z(j))
		call thetasub(thetak1(1), z(j), rk1(1))

		if (ENAflag(j) .eqv. .false.) then
			! Note: Without cross L-shell convection, solar forcing or other azimuthal asymmetries, ion motion is phi-invariant.
			phik1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1)
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			call phisub(phik1(1), x(j), y(j))
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

		if (ENAflag(j) .eqv. .false.) then
			if (phik1(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ION phik1= ', phik1(1), &
					' AND phiGCT VALUE= ', SpecieT(s)%FluxTubeT(f)%QCellT(1)%phiGCT(1), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
					// achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

		call qsub(qk1(1), rk1(1), thetak1(1))
		call psub(pk1(1), rk1(1), thetak1(1))
		call ellsub(ellk1(1), thetak1(1))
		call Bmagsub(Bmagk1(1), rk1(1), ellk1(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(Time(1))) .eqv. .true.) .or. (size(Time(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Time= ', Time(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(rk1(1))) .eqv. .true.) .or. &
			(size(rk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rk1= ', rk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ',	j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(thetak1(1))) .eqv. .true.) .or. &
			(size(thetak1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetak1= ', thetak1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(phik1(1))) .eqv. .true.) .or. &
			(size(phik1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phik1= ', phik1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(qk1(1))) .eqv. .true.) .or. &
			(size(qk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qk1= ', qk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(pk1(1))) .eqv. .true.) .or. &
			(size(pk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pk1= ', pk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(ellk1(1))) .eqv. .true.) .or. &
			(size(ellk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellk1= ', ellk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Bmagk1(1))) .eqv. .true.) .or. &
			(size(Bmagk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Bmagk1= ', Bmagk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE EACH PARTICLE CONFIGURATION-SPACE GRID CELL:

		if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then ! S. Magnetic Hemisphere
			QloopKUA1: do Qind= NqLB(1), NqUB(1), 1
				if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= qk1(1)) &
					.and. (qk1(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
					Qindk1(j)= Qind

					exit QloopKUA1

				else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) < qk1(1)) &
					.and. (qk1(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
					Qindk1(j)= Qind

					exit QloopKUA1

				end if
			end do QloopKUA1

			if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) > qk1(1)) then
				! SMH Lower boundary escape
				Qindk1(j)= 0d0
			else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) < qk1(1)) then
				! SMH Upper boundary escape
				Qindk1(j)= -1d0
			end if
		end if

		! ----------------------------------------------------

		if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then ! N. Magnetic Hemisphere
			QloopKUA2: do Qind= NqLB(1), NqUB(1), 1
				if ((Qind == NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= qk1(1)) &
					.and. (qk1(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
					Qindk1(j)= Qind

					exit QloopKUA2

				else if ((Qind /= NqLB(1)) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) > qk1(1)) &
					.and. (qk1(1) >= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then
					Qindk1(j)= Qind

					exit QloopKUA2

				end if
			end do QloopKUA2

			if (SpecieT(s)%FluxTubeT(f)%QCellT(NqLB(1))%qGLT(1) < qk1(1)) then
				! NMH Lower boundary escape
				Qindk1(j)= 0d0
			else if (SpecieT(s)%FluxTubeT(f)%QCellT(NqUB(1))%qGHT(1) > qk1(1)) then
				! NMH Upper boundary escape
				Qindk1(j)= -1d0
			end if
		end if

		! ----------------------------------------------------

		! CALL CARTESIAN FORCE SUBROUTINES AT INITIAL SPLIT TIME-STEP:

		call dBdssub(dBdsk1(1), rk1(1), thetak1(1), ellk1(1))

		! ----------------------------------------------------

		if (ENAflag(j) .eqv. .false.) then
			call Rperpsub(Rperpk1(1), pk1(1), phik1(1))
			call OmegaGsub(OmegaGk1(1), SpecieT(s)%qsT(1), SpecieT(s)%msT(1), Bmagk1(1))
			call musub(muk1(1), SpecieT(s)%msT(1), Bmagk1(1), Vperp(j))
		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then
			Rperpk1(1)= 0d0
			OmegaGk1(1)= 0d0
			muk1(1)= 0d0
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(dBdsk1(1))) .eqv. .true.) .or. &
			(size(dBdsk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dBdsk1= ', dBdsk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Rperpk1(1))) .eqv. .true.) .or. &
			(size(Rperpk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Rperpk1= ', Rperpk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(OmegaGk1(1))) .eqv. .true.) .or. &
			(size(OmegaGk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaGk1= ', OmegaGk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(muk1(1))) .eqv. .true.) .or. &
			(size(muk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' muk1= ', muk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN, SPHERICAL, AND DIPOLE MIRROR FORCE
		! ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.)) then

			call AMxsub(AMxk1p(1), SpecieT(s)%msT(1), muk1(1), &
				dBdsk1(1), thetak1(1), phik1(1), ellk1(1))
			call AMysub(AMyk1p(1), SpecieT(s)%msT(1), muk1(1), &
				dBdsk1(1), thetak1(1), phik1(1), ellk1(1))
			call AMzsub(AMzk1p(1), SpecieT(s)%msT(1), muk1(1), &
				dBdsk1(1), thetak1(1), ellk1(1))

			if (qk1(1) <= 0d0) then  ! S Mag Hemisphere
				AMpark1(1)= -abs(3d0*cos(thetak1(1))*sin(thetak1(1))*(AMxk1p(1)* &
					cos(phik1(1))+ AMyk1p(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
					AMzk1p(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1)))
			else if (qk1(1) > 0d0) then ! N Mag Hemisphere
				AMpark1(1)= abs(3d0*cos(thetak1(1))*sin(thetak1(1))*(AMxk1p(1)* &
					cos(phik1(1))+ AMyk1p(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
					AMzk1p(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1)))
			end if

			AMpk1(1)= 0d0
			AMphik1(1)= 0d0

			AMxk1(1)= cos(phik1(1))*(3d0*AMpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AMpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))- &
				AMphik1(1)*sin(phik1(1))
			AMyk1(1)= sin(phik1(1))*(3d0*AMpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AMpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))+ &
				AMphik1(1)*cos(phik1(1))
			AMzk1(1)= AMpark1(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))+ &
				3d0*AMpk1(1)*cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk1(1) <= 0) .and. (AMpark1(1) > 0d0)) .or. &
				((qk1(1) > 0) .and. (AMpark1(1) < 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AMpark1= ', AMpark1(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE A SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AMxk1p(1)= 0d0
			AMyk1p(1)= 0d0
			AMzk1p(1)= 0d0
			AMpark1(1)= 0d0
			AMpk1(1)= 0d0
			AMphik1(1)= 0d0
			AMxk1(1)= 0d0
			AMyk1(1)= 0d0
			AMzk1(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AMxk1p(1))) .eqv. .true.) .or. &
			(size(AMxk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk1p', AMxk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk1p(1))) .eqv. .true.) .or. &
			(size(AMyk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk1p', AMyk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk1p(1))) .eqv. .true.) .or. &
			(size(AMzk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk1p= ', AMzk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMphik1(1))) .eqv. .true.) .or. &
			(size(AMphik1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMphik1= ', AMphik1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpark1(1))) .eqv. .true.) .or. &
			(size(AMpark1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpark1= ', AMpark1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMpk1(1))) .eqv. .true.) .or. &
			(size(AMpk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMpk1= ', AMpk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMxk1(1))) .eqv. .true.) .or. &
			(size(AMxk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMxk1', AMxk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMyk1(1))) .eqv. .true.) .or. &
			(size(AMyk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMyk1', AMyk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AMzk1(1))) .eqv. .true.) .or. &
			(size(AMzk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AMzk1= ', AMzk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN, SPHERICAL, AND DIPOLE GRAVITATIONAL FORCE
		! ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.)) then

			if (qk1(1) <= 0d0) then  ! S Mag Hemisphere
				AGmagSk1(1)= 1d0
			else if (qk1(1) > 0d0) then ! N Mag Hemisphere
				AGmagSk1(1)= -1d0
			end if

			call AGxsub(AGxk1p(1), AGmagSk1(1), AGmag(j), thetak1(1), phik1(1), ellk1(1))
			call AGysub(AGyk1p(1), AGmagSk1(1), AGmag(j), thetak1(1), phik1(1), ellk1(1))
			call AGzsub(AGzk1p(1), AGmagSk1(1), AGmag(j), thetak1(1), ellk1(1))

			AGpark1(1)= 3d0*cos(thetak1(1))*sin(thetak1(1))*(AGxk1p(1)* &
				cos(phik1(1))+ AGyk1p(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
				AGzk1p(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))
			AGpk1(1)= 0d0
			AGphik1(1)= 0d0

			AGxk1(1)= cos(phik1(1))*(3d0*AGpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AGpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))- &
				AGphik1(1)*sin(phik1(1))
			AGyk1(1)= sin(phik1(1))*(3d0*AGpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AGpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))+ &
				AGphik1(1)*cos(phik1(1))
			AGzk1(1)= AGpark1(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))+ &
				3d0*AGpk1(1)*cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk1(1) <= 0) .and. (AGpark1(1) < 0d0)) .or. &
				((qk1(1) > 0) .and. (AGpark1(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark1= ', AGpark1(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE A SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if ((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.)) then

			call AGENAxsub(AGxk1(1), mNeut, rk1(1), &
				thetak1(1), phik1(1))
			call AGENAysub(AGyk1(1), mNeut, rk1(1), &
				thetak1(1), phik1(1))
			call AGENAzsub(AGzk1(1), mNeut, rk1(1), &
				thetak1(1))

			AGpark1(1)= 3d0*cos(thetak1(1))*sin(thetak1(1))*(AGxk1(1)* &
				cos(phik1(1))+ AGyk1(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
				AGzk1(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk1(1) <= 0) .and. (AGpark1(1) < 0d0)) .or. &
				((qk1(1) > 0) .and. (AGpark1(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AGpark1= ', AGpark1(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE A SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. (ENAflag(j) .eqv. .true.))) then

			AGxk1p(1)= 0d0
			AGyk1p(1)= 0d0
			AGzk1p(1)= 0d0
			AGpark1(1)= 0d0
			AGpk1(1)= 0d0
			AGphik1(1)= 0d0
			AGxk1(1)= 0d0
			AGyk1(1)= 0d0
			AGzk1(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AGxk1p(1))) .eqv. .true.) .or. &
			(size(AGxk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk1p= ', AGxk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk1p(1))) .eqv. .true.) .or. &
			(size(AGyk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk1p= ', AGyk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk1p(1))) .eqv. .true.) .or. &
			(size(AGzk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk1p= ', AGzk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGphik1(1))) .eqv. .true.) .or. &
			(size(AGphik1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGphik1= ', AGphik1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpark1(1))) .eqv. .true.) .or. &
			(size(AGpark1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpark1= ', AGpark1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGpk1(1))) .eqv. .true.) .or. &
			(size(AGpk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGpk1= ', AGpk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGxk1(1))) .eqv. .true.) .or. &
			(size(AGxk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGxk1= ', AGxk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGyk1(1))) .eqv. .true.) .or. &
			(size(AGyk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGyk1= ', AGyk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AGzk1(1))) .eqv. .true.) .or. &
			(size(AGzk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGzk1= ', AGzk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN, SPHERICAL, AND DIPOLE AMBIPOLAR ELECTRIC FIELD
		! ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (qk1(1) <= 0d0) then  ! S Mag Hemisphere
					AEAmagSk1(1)= -1d0
				else if (qk1(1) > 0d0) then ! N Mag Hemisphere
					AEAmagSk1(1)= 1d0
				end if
			end if
			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1) then
				AEAmagSk1(1)= 1d0
			end if

			call AEAxsub(AEAxk1p(1), AEAmagSk1(1), AEAmag(j), thetak1(1), phik1(1), ellk1(1))
			call AEAysub(AEAyk1p(1), AEAmagSk1(1), AEAmag(j), thetak1(1), phik1(1), ellk1(1))
			call AEAzsub(AEAzk1p(1), AEAmagSk1(1), AEAmag(j), thetak1(1), ellk1(1))

			AEApark1(1)= 3d0*cos(thetak1(1))*sin(thetak1(1))*(AEAxk1p(1)* &
				cos(phik1(1))+ AEAyk1p(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
				AEAzk1p(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))
			AEApk1(1)= 0d0
			AEAphik1(1)= 0d0

			AEAxk1(1)= cos(phik1(1))*(3d0*AEApark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AEApk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))- &
				AEAphik1(1)*sin(phik1(1))
			AEAyk1(1)= sin(phik1(1))*(3d0*AEApark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AEApk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))+ &
				AEAphik1(1)*cos(phik1(1))
			AEAzk1(1)= AEApark1(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))+ &
				3d0*AEApk1(1)*cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER AMBIPOLAR ACCELERATION SIGN:

			if (SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 0) then
				if (((qk1(1) <= 0) .and. (AEApark1(1) > 0d0)) .or. &
					((qk1(1) > 0) .and. (AEApark1(1) < 0d0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
						' AEApark1= ', AEApark1(1), 'qk1= ', qk1(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
						' IN KINETIC UPDATE A SUBROUTINE' // achar(27) // '[0m.'
				end if
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEAxk1p(1)= 0d0
			AEAyk1p(1)= 0d0
			AEAzk1p(1)= 0d0
			AEApark1(1)= 0d0
			AEApk1(1)= 0d0
			AEAphik1(1)= 0d0
			AEAxk1(1)= 0d0
			AEAyk1(1)= 0d0
			AEAzk1(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEAxk1p(1))) .eqv. .true.) .or. &
			(size(AEAxk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk1p', AEAxk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk1p(1))) .eqv. .true.) .or. &
			(size(AEAyk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk1p', AEAyk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk1p(1))) .eqv. .true.) .or. &
			(size(AEAzk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk1p', AEAzk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAphik1(1))) .eqv. .true.) .or. &
			(size(AEAphik1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAphik1= ', AEAphik1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApark1(1))) .eqv. .true.) .or. &
			(size(AEApark1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApark1= ', AEApark1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEApk1(1))) .eqv. .true.) .or. &
			(size(AEApk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEApk1= ', AEApk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAxk1(1))) .eqv. .true.) .or. &
			(size(AEAxk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAxk1', AEAxk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAyk1(1))) .eqv. .true.) .or. &
			(size(AEAyk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAyk1', AEAyk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEAzk1(1))) .eqv. .true.) .or. &
			(size(AEAzk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEAzk1', AEAzk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN, SPHERICAL, AND DIPOLE PARALLEL ELECTRIC FIELD
		! ACCELERATION COMPONENTS:

		if ((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) &
			.and. (ENAflag(j) .eqv. .false.)) then

			if (qk1(1) <= 0d0) then  ! S Mag Hemisphere
				AEPmagSk1(1)= 1d0
			else if (qk1(1) > 0d0) then ! N Mag Hemisphere
				AEPmagSk1(1)= -1d0
			end if

			if (n < EParlim) then
				AEPmag(j)= 0d0
			end if
			if (n >= EParlim) then
				AEPmag(j)= AEPmag(j)
			end if

			call AEPxsub(AEPxk1p(1), AEPmagSk1(1), AEPmag(j), thetak1(1), phik1(1), ellk1(1))
			call AEPysub(AEPyk1p(1), AEPmagSk1(1), AEPmag(j), thetak1(1), phik1(1), ellk1(1))
			call AEPzsub(AEPzk1p(1), AEPmagSk1(1), AEPmag(j), thetak1(1), ellk1(1))

			AEPpark1(1)= 3d0*cos(thetak1(1))*sin(thetak1(1))*(AEPxk1p(1)* &
				cos(phik1(1))+ AEPyk1p(1)*sin(phik1(1)))/sqrt(ellk1(1))+ &
				AEPzk1p(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))
			AEPpk1(1)= 0d0
			AEPphik1(1)= 0d0

			AEPxk1(1)= cos(phik1(1))*(3d0*AEPpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AEPpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))- &
				AEPphik1(1)*sin(phik1(1))
			AEPyk1(1)= sin(phik1(1))*(3d0*AEPpark1(1)*cos(thetak1(1))*sin(thetak1(1))+ &
				AEPpk1(1)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))+ &
				AEPphik1(1)*cos(phik1(1))
			AEPzk1(1)= AEPpark1(1)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))+ &
				3d0*AEPpk1(1)*cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER PARALLEL ACCELERATION SIGN:

			if (((qk1(1) <= 0) .and. (AEPpark1(1) < 0d0)) .or. &
				((qk1(1) > 0) .and. (AEPpark1(1) > 0d0))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT', &
					' AEPpark1= ', AEPpark1(1), ' WITH MAGNETIC HEMISPHERE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
					' IN KINETIC UPDATE A SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if (((SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) &
			.and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) then

			AEPxk1p(1)= 0d0
			AEPyk1p(1)= 0d0
			AEPzk1p(1)= 0d0
			AEPpark1(1)= 0d0
			AEPpk1(1)= 0d0
			AEPphik1(1)= 0d0
			AEPxk1(1)= 0d0
			AEPyk1(1)= 0d0
			AEPzk1(1)= 0d0

		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(AEPxk1p(1))) .eqv. .true.) .or. &
			(size(AEPxk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk1p', AEPxk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk1p(1))) .eqv. .true.) .or. &
			(size(AEPyk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk1p', AEPyk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk1p(1))) .eqv. .true.) .or. &
			(size(AEPzk1p(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk1p', AEPzk1p(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPphik1(1))) .eqv. .true.) .or. &
			(size(AEPphik1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPphik1= ', AEPphik1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpark1(1))) .eqv. .true.) .or. &
			(size(AEPpark1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpark1= ', AEPpark1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPpk1(1))) .eqv. .true.) .or. &
			(size(AEPpk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPpk1= ', AEPpk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPxk1(1))) .eqv. .true.) .or. &
			(size(AEPxk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPxk1', AEPxk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPyk1(1))) .eqv. .true.) .or. &
			(size(AEPyk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPyk1', AEPyk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(AEPzk1(1))) .eqv. .true.) .or. &
			(size(AEPzk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPzk1', AEPzk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! CALL ALL CARTESIAN, SPHERICAL, AND DIPOLE NET FORCE
		! ACCELERATION COMPONENTS FOR KINETIC AND FLUID FORCES:

		call Axsub(Axk1(1), AMxk1(1), AGxk1(1), AEAxk1(1), AEPxk1(1))
		call Aysub(Ayk1(1), AMyk1(1), AGyk1(1), AEAyk1(1), AEPyk1(1))
		call Azsub(Azk1(1), AMzk1(1), AGzk1(1), AEAzk1(1), AEPzk1(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(Axk1(1))) .eqv. .true.) .or. &
			(size(Axk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Axk1= ', Axk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Ayk1(1))) .eqv. .true.) .or. &
			(size(Ayk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Ayk1= ', Ayk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Azk1(1))) .eqv. .true.) .or. &
			(size(Azk1(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Azk1= ', Azk1(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE SPHERICAL VELOCITY COMPONENTS:

		Vr(1)= Vx(j)*sin(thetak1(1))*cos(phik1(1))+ &
			Vy(j)*sin(thetak1(1))*sin(phik1(1))+ &
			Vz(j)*cos(thetak1(1))
		Vtheta(1)= Vx(j)*cos(thetak1(1))*cos(phik1(1))+ &
			Vy(j)*cos(thetak1(1))*sin(phik1(1))- &
			Vz(j)*sin(thetak1(1))

		! ----------------------------------------------------

		! RE-ALIGN VELOCITY VECTOR INTO LOCAL DIPOLE COORDINATES:

		Vpar(j)= 3d0*cos(thetak1(1))*sin(thetak1(1))*(Vx(j)* &
			cos(phik1(1))+ Vy(j)*sin(phik1(1)))/sqrt(ellk1(1))+ &
			Vz(j)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))

		! ----------------------------------------------------

		EpsilonPar(1)= (6.242d18)*(SpecieT(s)%msT(1)/2d0)*Vpar(j)**2d0 ! Parallel kinetic energy [eV]

		if (ENAflag(j) .eqv. .false.) then
			Vp(j)= 0d0
			Vphi(j)= 0d0
			EpsilonPerp(1)= (6.242d18)*(SpecieT(s)%msT(1)/2d0)*Vperp(j)**2d0 ! Perp kinetic energy [eV]
			alpha(1)= atan2(Vperp(j), Vpar(j)) ! Velocity pitch angle [rads]

			if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

				VparloopKUA1: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

					if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, 1)%VparGLT(1) > Vpar(j)) &
						.or. (Vpar(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
						Vparindk1(j)= 0d0

						exit VparloopKUA1

					end if
					if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, Vparind)%VparGLT(1) <= Vpar(j)) &
						.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
						Vparindk1(j)= Vparind

						exit VparloopKUA1

					end if
					if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, Vparind)%VparGLT(1) < Vpar(j)) &
						.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
						Vparindk1(j)= Vparind

						exit VparloopKUA1

					end if
				end do VparloopKUA1

				Vperp1loopKUA: do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
					if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, 1)%Vperp1GLT(1) > Vperp1(j)) &
						.or. (Vperp1(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1))) then
						Vperp1indk1(j)= 0d0

						exit Vperp1loopKUA

					end if
					if ((Vperp1ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) <= Vperp1(j)) &
						.and. (Vperp1(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
						Vperp1indk1(j)= Vperp1ind

						exit Vperp1loopKUA

					end if
					if ((Vperp1ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) < Vperp1(j)) &
						.and. (Vperp1(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
						Vperp1indk1(j)= Vperp1ind

						exit Vperp1loopKUA

					end if
				end do Vperp1loopKUA

				Vperp2loopKUA: do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
					if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, 1, 1)%Vperp2GLT(1) > Vperp2(j)) &
						.or. (Vperp2(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1))) then
						Vperp2indk1(j)= 0d0

						exit Vperp2loopKUA

					end if
					if ((Vperp2ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) <= Vperp2(j)) &
						.and. (Vperp2(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
						Vperp2indk1(j)= Vperp2ind

						exit Vperp2loopKUA

					end if
					if ((Vperp2ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) < Vperp2(j)) &
						.and. (Vperp2(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
						Vperp2indk1(j)= Vperp2ind

						exit Vperp2loopKUA

					end if
				end do Vperp2loopKUA

			else

				VparloopKUA2: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
					if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(1, Vparind)%VparGLT(1) <= Vpar(j)) &
						.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(1, Vparind)%VparGHT(1))) then
						Vparindk1(j)= Vparind

						exit VparloopKUA2

					end if
					if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(1, Vparind)%VparGLT(1) < Vpar(j)) &
						.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(1, Vparind)%VparGHT(1))) then
						Vparindk1(j)= Vparind

						exit VparloopKUA2

					end if
				end do VparloopKUA2

				VperploopKUA: do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
					if ((Vperpind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(Vperpind, 1)%VperpGLT(1) <= Vperp(j)) &
						.and. (Vperp(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(Vperpind, 1)%VperpGHT(1))) then
						Vperpindk1(j)= Vperpind

						exit VperploopKUA

					end if
					if ((Vperpind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(Vperpind, 1)%VperpGLT(1) < Vperp(j)) &
						.and. (Vperp(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(Vperpind, 1)%VperpGHT(1))) then
						Vperpindk1(j)= Vperpind

						exit VperploopKUA

					end if
				end do VperploopKUA

			end if

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR ALL ION Vpar VALUES WITHIN VELOCITY-SPACE GRID:

			!if ((Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	VCellT(1, 1)%VparGLT(1)) .or. (Vpar(j) > SpecieT(s)%FluxTubeT(f)% &
			!	QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' Vpar= ', Vpar(j), ' ION VALUE OUT OF GRID WITH VparGL= ', &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VparGLT(1), &
			!		' AND VparGH= ', SpecieT(s)%FluxTubeT(f)% &
			!		QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1), &
			!		'FOR SPECIE= ', s, &
			!		', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
			!		', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
			!		// achar(27) // '[0m.'
				!call MPI_FINALIZE(ierr)
				!stop
			!end if

			! ----------------------------------------------------

		else if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.) .and. (Qindk1(j) /= -1) .and. (Qindk1(j) /= 0)) then
			Vp(j)= abs((1d0- 3d0*(cos(thetak1(1))**2d0))*(Vx(j)* &
				cos(phik1(1))+ Vy(j)*sin(phik1(1)))/sqrt(ellk1(1))+ 3d0*Vz(j)* &
				cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1)))
			Vphi(j)= -Vx(j)*sin(phik1(1))+ Vy(j)*cos(phik1(1))
			EpsilonPerp(1)= 0d0 ! Perp kinetic energy [eV]
			alpha(1)= 0d0 ! Velocity pitch angle [rads]

			VploopKUA: do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
				if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, 1, 1)%VpGLT(1) > Vp(j)) &
					.and. (Vp(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1, 1)%VpGHT(1))) then
					Vpindk1(j)= 0d0

					exit VploopKUA

				end if
				if ((Vpind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(Vpind, 1, 1)%VpGLT(1) <= Vp(j)) &
					.and. (Vp(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(Vpind, 1, 1)%VpGHT(1))) then
					Vpindk1(j)= Vpind

					exit VploopKUA

				end if
				if ((Vpind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(Vpind, 1, 1)%VpGLT(1) < Vp(j)) &
					.and. (Vp(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(Vpind, 1, 1)%VpGHT(1))) then
					Vpindk1(j)= Vpind

					exit VploopKUA

				end if
			end do VploopKUA

			VqloopKUA: do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
				if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, 1, 1)%VqGLT(1) > Vpar(j)) &
					.and. (Vpar(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1)%VqGHT(1))) then
					Vqindk1(j)= 0d0

					exit VqloopKUA

				end if
				if ((Vqind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, Vqind, 1)%VqGLT(1) <= Vpar(j)) &
					.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, Vqind, 1)%VqGHT(1))) then
					Vqindk1(j)= Vqind

					exit VqloopKUA

				end if
				if ((Vqind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, Vqind, 1)%VqGLT(1) < Vpar(j)) &
					.and. (Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, Vqind, 1)%VqGHT(1))) then
					Vqindk1(j)= Vqind

					exit VqloopKUA

				end if
			end do VqloopKUA

			VphiloopKUA: do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1
				if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, 1, 1)%VphiGLT(1) > Vphi(j)) &
					.and. (Vphi(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1))%VphiGHT(1))) then
					Vphiindk1(j)= 0d0

					exit VphiloopKUA

				end if
				if ((Vphiind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, 1, Vphiind)%VphiGLT(1) <= Vphi(j)) &
					.and. (Vphi(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, 1, Vphiind)%VphiGHT(1))) then
					Vphiindk1(j)= Vphiind

					exit VphiloopKUA

				end if
				if ((Vphiind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					V3CellT(1, 1, Vphiind)%VphiGLT(1) < Vphi(j)) &
					.and. (Vphi(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(1, 1, Vphiind)%VphiGHT(1))) then
					Vphiindk1(j)= Vphiind

					exit VphiloopKUA

				end if
			end do VphiloopKUA

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR ALL ENA Vp, Vq, Vphi VALUES WITHIN VELOCITY-SPACE GRID:

			!if ((Vpar(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(1, 1, 1)%VqGLT(1)) &
			!	.or. (Vpar(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VqGHT(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' Vpar= ', Vpar(j), ' ENA VALUE OUT OF GRID WITH VparGL= ', &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(1, 1, 1)%VqGLT(1), &
			!		' AND VparGH= ', SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VqGHT(1), &
			!		'FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
			!		', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
			!		// achar(27) // '[0m.'
				!call MPI_FINALIZE(ierr)
				!stop
			!end if

			!if ((Vp(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(1, 1, 1)%VpGLT(1)) &
			!	.or. (Vp(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VpGHT(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' Vp= ', Vp(j), ' ENA VALUE OUT OF GRID WITH VpGL= ', &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(1, 1, 1)%VpGLT(1), &
			!		' AND VpGH= ', SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VpGHT(1), &
			!		'FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
			!		', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
			!		// achar(27) // '[0m.'
				!call MPI_FINALIZE(ierr)
				!stop
			!end if

			!if ((Vphi(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(1, 1, 1)%VphiGLT(1)) &
			!	.or. (Vphi(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!	V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VphiGHT(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' Vphi= ', Vphi(j), ' ENA VALUE OUT OF GRID WITH VphiGL= ', &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(1, 1, 1)%VphiGLT(1), &
			!		' AND VphiGH= ', SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
			!		V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))%VphiGHT(1), &
			!		'FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
			!		', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
			!		// achar(27) // '[0m.'
				!call MPI_FINALIZE(ierr)
				!stop
			!end if

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		! COMPUTE NEW FIELD-ALIGNED CARTESIAN VELOCITY COMPONENTS BY COMPUTATIONAL L-SHELL DRIFT LIMIT:

		VxR(1)= cos(phik1(1))*(3d0*Vpar(j)*cos(thetak1(1))*sin(thetak1(1))+ &
			Vp(j)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))- &
			Vphi(j)*sin(phik1(1))
		VyR(1)= sin(phik1(1))*(3d0*Vpar(j)*cos(thetak1(1))*sin(thetak1(1))+ &
			Vp(j)*(1d0- 3d0*(cos(thetak1(1))**2d0)))/sqrt(ellk1(1))+ &
			Vphi(j)*cos(phik1(1))
		VzR(1)= Vpar(j)*(3d0*(cos(thetak1(1))**2d0)- 1d0)/sqrt(ellk1(1))+ &
			3d0*Vp(j)*cos(thetak1(1))*sin(thetak1(1))/sqrt(ellk1(1))

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

		if ((isnan(real(Vr(1))) .eqv. .true.) .or. &
			(size(Vr(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vr= ', Vr(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vtheta(1))) .eqv. .true.) .or. &
			(size(Vtheta(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vtheta= ', Vtheta(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vphi(j))) .eqv. .true.) .or. &
			(size(Vphi(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vphi= ', Vphi(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vpar(j))) .eqv. .true.) .or. &
			(size(Vpar(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vpar= ', Vpar(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(Vp(j))) .eqv. .true.) .or. &
			(size(Vp(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vp= ', Vp(j), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(EpsilonPar(1))) .eqv. .true.) .or. &
			(size(EpsilonPar(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EpsilonPar= ', EpsilonPar(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(EpsilonPerp(1))) .eqv. .true.) .or. &
			(size(EpsilonPerp(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EpsilonPerp= ', EpsilonPerp(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(alpha(1))) .eqv. .true.) .or. &
			(size(alpha(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' alpha= ', alpha(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(VxR(1))) .eqv. .true.) .or. &
			(size(VxR(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VxR= ', VxR(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(VyR(1))) .eqv. .true.) .or. &
			(size(VyR(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VyR= ', VyR(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(VzR(1))) .eqv. .true.) .or. &
			(size(VzR(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VzR= ', VzR(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN KINETIC UPDATE A SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! PERFORM ICR WAVE HEATING:

		call WaveHeatingSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if ((n == 1) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S811End)
				write(S811string, '(i10)')  nint(S811End)
				write(*, *) trim('%% 8.11.1- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S811string)) // &
					trim(' s. INITIAL KINETIC UPDATE A COMPLETE %%')
			end if
		end if

		if (rank == 0) then
			if ((n == 2) .and. (j == SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				call cpu_time(S811End)
				write(S811string, '(i10)')  nint(S811End)
				write(*, *) trim('%% 8.11.1- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S811string)) // &
					trim(' s. SECOND KINETIC UPDATE A COMPLETE %%')
			end if
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticUpdateASub

end module KineticUpdateA
