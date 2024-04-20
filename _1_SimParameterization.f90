module SimParameterization

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	1- SIMULATION PARAMETERIZATION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use ConfigGridGenerator
use VelGridGenerator

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SET ALL PHASE-SPACE GRID, INITIALIZATION, AND SIMULATION PARAMETERS PER PARTICLE SPECIES,
! FLUX TUBE, AND PHASE-SPACE GRID CELL:

	subroutine SimParameterizationSub

		! ----------------------------------------------------

		! ALLOCATE PARTICLE SPECIES DERIVED DATA TYPE:

		i= (0, 1) ! Define sqrt(-1)
		allocate(SpecieT(Stot))

		! ----------------------------------------------------

		! ALLOCATE FLUX TUBE DERIVED DATA TYPE NESTED IN PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1

			SpecieT(s)%NfT(1)= Nf
			allocate(SpecieT(s)%FluxTubeT(SpecieT(s)%NfT(1)))

		end do

		! ----------------------------------------------------

		! SET RANGE ALONG FLUX TUBES FOR INITIAL POSITIONS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				if (s == 1) then
					! Starting and ending Q cell for density initialization
					! (including selected cells and does not include lower boundary ghost cell)

					NqLB(1)= NqICA
					NqUB(1)= NqICB

					NqGpF= (NqUB(1)- NqLB(1)+ 4)

				end if

				NqIC(1)= abs(NqICB- NqICA)+ 1d0 ! Q cell range for density initialization

				SpecieT(s)%FluxTubeT(f)%NqICAT= NqICA ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%NqICBT= NqICB
				SpecieT(s)%FluxTubeT(f)%NqICT= NqIC

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((NqICA /= SpecieT(s)%FluxTubeT(f)%NqICAT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NqICAT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICAT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICAT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((NqICB /= SpecieT(s)%FluxTubeT(f)%NqICBT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NqICBT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICBT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICBT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((NqIC(1) /= SpecieT(s)%FluxTubeT(f)%NqICT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NqICT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! ALLOCATE INITIALIZATION CONFIG-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE
		! TYPES NESTED IN PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! Allocate QCellICT(Qind) derived data type nested in FluxTubeT(f)
				! nested in SpecieT(s)
				allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

			end do
		end do

		! ----------------------------------------------------

		! SET PARTICLE SPECIES INITIAL DENSITY PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				SpecieT(s)%FluxTubeT(f)%zns0T= zns0 ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%ns0T= ns0
				SpecieT(s)%FluxTubeT(f)%nsnormfacT= nsnormfac

				SpecieT(s)%FluxTubeT(f)%zns0NeutT= zns0Neut ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%ns0NeutT= ns0Neut

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((zns0 /= SpecieT(s)%FluxTubeT(f)%zns0T(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%zns0T(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%zns0T(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zns0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ns0 /= SpecieT(s)%FluxTubeT(f)%ns0T(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ns0T(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ns0T(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ns0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((nsnormfac /= SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%nsnormfacT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormfacT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((zns0Neut /= SpecieT(s)%FluxTubeT(f)%zns0NeutT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%zns0NeutT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%zns0NeutT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zns0NeutT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ns0Neut /= SpecieT(s)%FluxTubeT(f)%ns0NeutT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ns0NeutT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ns0NeutT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ns0NeutT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! ALLOCATE CONFIGURATION-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE TYPES NESTED IN
		! PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! Allocate QCellT(Qind) derived data type nested in FluxTubeT(f)
				! nested in SpecieT(s)
				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(((NqUB(1)- NqLB(1))+ 1)))
				allocate(SpecieT(s)%FluxTubeT(f)%QCell0T(((NqUB(1)- NqLB(1))+ 3)))

			end do
		end do

		! ----------------------------------------------------

		! CONSTRUCT INITIAL CONFIGURATION-SPACE GRID:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				call ConfigGridGeneratorSub

			end do
		end do

		! ----------------------------------------------------

		! CONFIGURATION-SPACE GRID DIMENSIONS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				do Qind= NqLB(1), NqUB(1)+ 2, 1
					if ((size(SpecieT(s)%FluxTubeT(f)%Te0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%Te0T(Qind))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0T HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (SpecieT(s)%FluxTubeT(f)%Te0T(Qind) > dNTeEND) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INITIAL Te0T IS GREATER THAN Te CAP FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if
				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! PARTICLE MASS, CHARGE, AND ATOMIC RADIUS PER PARTICLE SPECIES:

		do s= 1, Stot, 1

			SpecieT(s)%msT= mO
			SpecieT(s)%qsT= qO

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(SpecieT(s)%Qindns0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%Qindns0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0T HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((mO /= SpecieT(s)%msT(1)) .or. (size(SpecieT(s)%msT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%msT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' msT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((qO /= SpecieT(s)%qsT(1)) .or. (size(SpecieT(s)%qsT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%qsT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qsT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		! PARTICLE SPECIES INITIAL TEMPERATURES AND CONFIGURATION-SPACE GRID
 		! PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! COMPUTE TOTAL FIELD-LINE ARC LENGTH:

				allocate(SpecieT(s)%FluxTubeT(f)%ellqCT(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))
				SpecieT(s)%FluxTubeT(f)%ellqCT(:)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)*SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)
				SpecieT(s)%FluxTubeT(f)%SUMellqCT(1)= sum(SpecieT(s)%FluxTubeT(f)%ellqCT(:))

				! ----------------------------------------------------

				do Qind= NqLB(1), NqUB(1)+ 2, 1

					! ----------------------------------------------------

					! Total initial ion temperature [K]
					! Note: Isotropic temperature from grid data (Ti= Tperp= Tpar= Tperp1= Tperp2)
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%Ts0T(1)= (1d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))+ &
						(2d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerpT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsParT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGLT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGHT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) < 0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE ', &
							' d3xCT FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', &
							Qind, ' IN SIMULATION PARAMETERIZATION SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET SIMULATION DURATION AND TIME-STEP FOR MOMENT COMPUTATION:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------
				h(1)= (B- A)/Nt ! Time-step size [s]

				! Injection time-step [s] (must be > 2 and excludes initial time-step)
				! Mean transit time through LB ghost cell
				Q0ndatfac(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1))%hqC0T(1)* &
					SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1))%dqC0T(1)* &
					(abs(sqrt(SpecieT(s)%msT(1)/ &
					((kB/1d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1))%TsPar0T(1))))))/h(1)

				Q0NNt(1)= Nt/Q0ndatfac(1) ! Number of time-steps for injection
				!Q0NNt(1)= Nt/10d0 ! Number of time-steps for injection

				! Time-step interval for moment computation (must be > 2 and excludes initial time-step)
				NNt(1)= Q0NNt(1)
				ndatfac(1)= Nt/NNt(1)

				SpecieT(s)%FluxTubeT(f)%AT= A ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%BT= B
				SpecieT(s)%FluxTubeT(f)%NtT= Nt
				SpecieT(s)%FluxTubeT(f)%hT= h
				SpecieT(s)%FluxTubeT(f)%Q0ndatfacT= Q0ndatfac
				SpecieT(s)%FluxTubeT(f)%Q0NNtT= Q0NNt
				SpecieT(s)%FluxTubeT(f)%NNtT= NNt
				SpecieT(s)%FluxTubeT(f)%ndatfacT= ndatfac

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((A /= SpecieT(s)%FluxTubeT(f)%AT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%AT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%AT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((B /= SpecieT(s)%FluxTubeT(f)%BT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%BT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%BT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((Nt /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NtT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NtT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NtT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((h(1) /= SpecieT(s)%FluxTubeT(f)%hT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%hT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%hT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((Q0NNt(1) /= SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%Q0NNtT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Q0NNtT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((Q0ndatfac(1) /= SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Q0ndatfacT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((NNt(1) /= SpecieT(s)%FluxTubeT(f)%NNtT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NNtT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NNtT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NNtT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ndatfac(1) /= SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ndatfacT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ndatfacT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(s)%FluxTubeT(f)%ndatfacT(1) < 3d0) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ndatfacT= ', &
					SpecieT(s)%FluxTubeT(f)%ndatfacT(1), &
					' IS NOT COMPATIBLE WITH KINETIC SOLVER TIME FORMAT FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1) < 3d0) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Q0ndatfacT= ', &
					SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1), &
					' IS NOT COMPATIBLE WITH KINETIC SOLVER TIME FORMAT FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET SIMULATION FLAGS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! ION POPULATION FLAGS:

				! Create nested derived data types
 				SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT= FLUIDIONEXPORTflag
				SpecieT(s)%FluxTubeT(f)%SYMVPARflagT= SYMVPARflag
				SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT= LBCONDITIONflag
				SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT= UBCONDITIONflag
				SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT= LBREPLENISHflag
				SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT= UBREPLENISHflag
				SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT= DENSITYPROFILEflag
				SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT= STATICINJECTIONflag
				SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT= DENSITYOUTPUTflag
				SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT= DENSITYINPUTflag

				SpecieT(s)%FluxTubeT(f)%IONNOISEflagT= IONNOISEflag
 				SpecieT(s)%FluxTubeT(f)%ICRflagT= ICRflag
				SpecieT(s)%FluxTubeT(f)%MIRRORflagT= MIRRORflag
				SpecieT(s)%FluxTubeT(f)%GRAVflagT= GRAVflag
				SpecieT(s)%FluxTubeT(f)%EAMBflagT= EAMBflag
				SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT= EAMBSELFCONSISTflag
				SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT= EAMBSIGNflag
				SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT= EAINERTIALflag
				SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT= EAPRESSUREflag

				SpecieT(s)%FluxTubeT(f)%EPARflagT= EPARflag

				SpecieT(s)%FluxTubeT(s)%MOMENTFILTERflagT= MOMENTFILTERflag
				SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT= ION2VPERPflag
				SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT= PHASEIONDISTRIBflag
				SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT= PHASEDENSITYIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT= PHASEVELPERPIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT= PHASEVELPARIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT= PHASEENERGYIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT= &
					PHASEENERGYPERPIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT= &
					PHASEENERGYPARIONMOMENTflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT= &
					PHASEENERGYPERPIONMOMENTCENTERflag
				SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT= &
					PHASEENERGYPARIONMOMENTCENTERflag
				SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT= FLUIDIONREFflag

				SpecieT(s)%FluxTubeT(f)%ENANOISEflagT= ENANOISEflag
				SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT= QEXCHANGEflag

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((FLUIDIONEXPORTflag /= SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' FLUIDIONEXPORTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SYMVPARflag /= SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' SYMVPARflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((LBCONDITIONflag /= SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBCONDITIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((UBCONDITIONflag /= SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBCONDITIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((LBREPLENISHflag /= SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBREPLENISHflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((UBREPLENISHflag /= SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBREPLENISHflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((DENSITYPROFILEflag /= SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DENSITYPROFILEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((STATICINJECTIONflag /= SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' STATICINJECTIONflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((DENSITYOUTPUTflag /= SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DENSITYOUTPUTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((DENSITYINPUTflag /= SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DENSITYINPUTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((IONNOISEflag /= SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' IONNOISEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ICRflag /= SpecieT(s)%FluxTubeT(f)%ICRflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ICRflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ICRflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ICRflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((MIRRORflag /= SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%MIRRORflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' MIRRORflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((GRAVflag /= SpecieT(s)%FluxTubeT(f)%GRAVflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%GRAVflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%GRAVflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' GRAVflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBflag /= SpecieT(s)%FluxTubeT(f)%EAMBflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBSELFCONSISTflag /= SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBSELFCONSISTflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAMBSIGNflag /= SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAMBSIGNflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAINERTIALflag /= SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAINERTIALflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EAPRESSUREflag /= SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EAPRESSUREflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((EPARflag /= SpecieT(s)%FluxTubeT(f)%EPARflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%EPARflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%EPARflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EPARflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((MOMENTFILTERflag /= SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' MOMENTFILTERflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ION2VPERPflag /= SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION2VPERPflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEIONDISTRIBflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEIONDISTRIBflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' PHASEIONDISTRIBflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEDENSITYIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEDENSITYIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEDENSITYIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEVELPERPIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEVELPERPIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEVELPERPIONMOMENTflagflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEVELPARIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEVELPARIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEVELPARIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPERPIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPERPIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPARIONMOMENTflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPARIONMOMENTflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPERPIONMOMENTCENTERflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTCENTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPERPIONMOMENTCENTERflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPERPIONMOMENTCENTERflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((PHASEENERGYPARIONMOMENTCENTERflag /= SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTCENTERflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)% &
					PHASEENERGYPARIONMOMENTCENTERflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' PHASEENERGYPARIONMOMENTCENTERflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((FLUIDIONREFflag /= SpecieT(s)%FluxTubeT(f)% &
					FLUIDIONREFflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' FLUIDIONREFflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((ENANOISEflag /= SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENANOISEflagT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((QEXCHANGEflag /= SpecieT(s)%FluxTubeT(f)% &
					QEXCHANGEflagT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(:)) /= 1) &
					.or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1))) &
					.eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' QEXCHANGEflagT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
					' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR MOMENT FLAG CONSISTENCY:

				if (((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 1)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1) == 1))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND ION MOMENT FLAGS ', &
						' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
						' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT LBREPLENISHflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT UBREPLENISHflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT UBCONDITIONflagT AND DENSITYPROFILEflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEVELPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEVELPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND PHASEENERGYPERPIONMOMENTCENTERflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND PHASEENERGYPARIONMOMENTCENTERflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND FLUIDIONREFflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEIONDISTRIBflagT AND FLUIDIONEXPORTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0)) .or. &
					((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAINETIALflagT AND EAPRESSUREflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAINETIALflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAPRESSUREflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBSELFCONSISTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT EAMBSIGNflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEDENSITYIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEVELPERPIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEVELPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT MOMENTFILTERflagT AND PHASEENERGYPARIONMOMENTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 0) .and. &
					(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT PHASEENERGYIONMOMENTflagT AND QEXCHANGEflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT ENANOISEflagT AND QEXCHANGEflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT DENSITYOUTPUTflagT AND DENSITYINPUTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1) .and. &
					(SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 0)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT DENSITYOUTPUTflagT AND FLUIDIONEXPORTflagT FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET NOISE FILTERS FOR ION MOMENTS:

		if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS CONSISTENT MOMENT FILTER MOVING AVERAGE POINTS:

					if (M0MAfilterPt > ((NqUB(1)- NqLB(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION DENSITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1Perp1MAfilterPt > ((NqUB(1)- NqLB(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PERP1 VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1Perp2MAfilterPt > ((NqUB(1)- NqLB(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PERP2 VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M1ParMAfilterPt > ((NqUB(1)- NqLB(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PARALLEL VELOCITY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (M2ParMAfilterPt > ((NqUB(1)- NqLB(1))+ 1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION PARALLEL ENERGY MOMENT FILTER', &
							' MOVING AVERAGE POINT IS GREATER THAN NUMBER OF CONFIGURATION-SPACE GRID CELLS', &
							' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do
			end do
		end if

		! ----------------------------------------------------

		! SET NOISE LIMIT FOR ION POPULATION STATISTICS:

		if (SpecieT(1)%FluxTubeT(1)%IONNOISEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1 ! 0.2d0 limits N= 25 and 0.141 limits N= 50
					IonNoiseLimit(1)= 1d0/(sqrt(IonNoiseLimitNph)) ! Noise limit for ion distribution functions
				end do
			end do
		end if
		if (SpecieT(1)%FluxTubeT(1)%IONNOISEflagT(1) == 0) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					IonNoiseLimit(1)= 0d0 ! Noise limit for ion distribution functions
				end do
			end do
		end if
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT= IonNoiseLimit

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((IonNoiseLimit(1) /= SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' IonNoiseLimitT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET NOISE LIMIT FOR ENA POPULATION STATISTICS:

		if (SpecieT(1)%FluxTubeT(1)%ENANOISEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1 ! 0.2d0 limits N= 25 and 0.141 limits N= 50
						ENANoiseLimit(1)= 1d0/(sqrt(ENANoiseLimitNph)) ! Noise limit for ENA distribution functions
				end do
			end do
		end if
		if (SpecieT(1)%FluxTubeT(1)%ENANOISEflagT(1) == 0) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					ENANoiseLimit(1)= 0d0 ! Noise limit for ENA distribution functions
				end do
			end do
		end if
		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT= ENANoiseLimit

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((ENANoiseLimit(1) /= SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ENANoiseLimitT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
					' SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! IMPORT EQUILIBRIUM DENSITY PARAMETERS FROM PREVIOUS SIMULATION:

		if (SpecieT(1)%FluxTubeT(1)%DENSITYINPUTflagT(1) == 1) then

			! ----------------------------------------------------

			allocate(DensityInput(NqUB(1)- NqLB(1)+ 1), TemperatureInput(NqUB(1)- NqLB(1)+ 1), &
				EAInertialInput(NqUB(1)- NqLB(1)+ 1), EAPressureInput(NqUB(1)- NqLB(1)+ 1), &
				EAmagInput(NqUB(1)- NqLB(1)+ 1))

			! ----------------------------------------------------

			do s= 1, Stot, 1
				write(sstring, '(I5)') s
				do f= 1, SpecieT(s)%NfT(1), 1
					write(fstring, '(I5)') f

					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'nsnormCLBOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) nsnormCLBInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'nsnormCUBOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) nsnormCUBInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'DensityOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) DensityInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'TemperatureOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) TemperatureInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAInertialOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAInertialInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAPressureOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAPressureInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'EAmagOutputRTfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) EAmagInput(:)
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'LBREPLENISHflagTOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) LBREPLENISHflagInput
					close(0)
					open (unit= 0, file= Densitydatadir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // 'UBREPLENISHflagTOutputfort.bin')), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) UBREPLENISHflagInput
					close(0)

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%nsnormCLBInputT(1)= nsnormCLBInput
					SpecieT(s)%FluxTubeT(f)%nsnormCUBInputT(1)= nsnormCUBInput

					do Qind= NqLB(1), NqUB(1), 1
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%DensityInputT(1)= DensityInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TemperatureInputT(1)= TemperatureInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAInertialInputT(1)= EAInertialInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAPressureInputT(1)= EAPressureInput(Qind)
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EAmagInputT(1)= EAmagInput(Qind)
					end do

					! ----------------------------------------------------

					if (LBREPLENISHflagInput /= SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCONSISTENT LBREPLENISHflagInputT FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if
					if (UBREPLENISHflagInput /= SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCONSISTENT UBREPLENISHflagInputT FOR SPECIE= ', s, &
							' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do
			end do
		end if

		! ----------------------------------------------------

		! CONSTRUCT INITIAL VELOCITY-SPACE GRID:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				call VelGridGeneratorSub

			end do
		end do

		! ----------------------------------------------------

		! ION VELOCITY-SPACE GRID DIMENSIONS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp1GT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp2GT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				else

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperpGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

				do Qind= NqLB(1), NqUB(1), 1

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVparGT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! ENA VELOCITY-SPACE GRID DIMENSIONS:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVpGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVqGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVphiGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do
				end do
			end do
		end if

		! ----------------------------------------------------

		! COMPUTE NUMBER OF MPI RANKS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				call MPI_COMM_SIZE(MPI_COMM_WORLD, ranksize(1), ierr)
			end do
		end do

		allocate(RankT(ranksize(1)))

		! ----------------------------------------------------

		! ION EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do Qind= NqLB(1), NqUB(1), 1
					if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then

						! ----------------------------------------------------

						do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
								do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

									! ----------------------------------------------------

									if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(1) == 0d0) then
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(1)= 1d-15
									end if

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp1GT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp2GT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3vCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do
						end do

						! ----------------------------------------------------

					else

						! ----------------------------------------------------

						do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
							do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

								! ----------------------------------------------------

								if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCT(1) == 0d0) then
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCT(1)= 1d-15
								end if

								! ----------------------------------------------------

								! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVperpGT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVperpGT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperpGT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVparGT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVparGT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGLT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGLT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGLT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGHT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGHT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGHT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGCT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGCT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGCT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGLT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGLT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGHT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGHT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGCT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGCT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%d3vCT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%d3vCT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3vCT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end do
						end do

						! ----------------------------------------------------

					end if

				end do
			end do
		end do

		! ----------------------------------------------------

		! ENA EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1
					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
							do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
								do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS,
									! SIZES, AND FINITE VALUES:

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VpGLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VpGHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VpGCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VphiGLT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', &
											Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VphiGHT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE FOR', &
											' SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', Vpind= ', Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VphiGCT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', &
											Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hVpCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hVqCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' hVphiCT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
											' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', &
											Qind, ', Vpind= ', Vpind, ', Vqind= ', Vqind, &
											', AND Vphiind= ', Vphiind, ' IN SIMULATION', &
											' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVpCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVpCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
											', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
											' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVqCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVqCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' dVqCT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE FOR', &
											' SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', Vpind= ', Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' dVphiCT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE FOR', &
											' SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', Vpind= ', Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%d33vCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' d33vCT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE FOR', &
											' SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
											', Vpind= ', Vpind, ', Vqind= ', Vqind, ', AND Vphiind= ', &
											Vphiind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do
							end do
						end do

						! ----------------------------------------------------

					end do
				end do
			end do
		end if

		! ----------------------------------------------------

		! SET ALL DRIFT LIMITS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				! L-shell drift limit over each time-step [RE]
				pdriftLim(1)= SpecieT(s)%FluxTubeT(f)%QCellT(1)%dpCT(1)
				! ENA R drift limit over entire simulation [m]
				rdriftLim(1)= 5d0*RE
				! ENA THETA drift limit over entire simulation [rads]
				thetadriftLim(1)= 5d-2
				! ENA PHI drift limit over entire simulation [rads]
				phidriftLim(1)= 5d-2

				SpecieT(s)%FluxTubeT(f)%pdriftLimT= pdriftLim ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%rdriftLimT= rdriftLim
				SpecieT(s)%FluxTubeT(f)%thetadriftLimT= thetadriftLim
				SpecieT(s)%FluxTubeT(f)%phidriftLimT= phidriftLim

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((pdriftLim(1) /= SpecieT(s)%FluxTubeT(f)%pdriftLimT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%pdriftLimT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%pdriftLimT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pdriftLimT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((rdriftLim(1) /= SpecieT(s)%FluxTubeT(f)%rdriftLimT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%rdriftLimT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%rdriftLimT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rdriftLimT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((thetadriftLim(1) /= SpecieT(s)%FluxTubeT(f)%thetadriftLimT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%thetadriftLimT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%thetadriftLimT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetadriftLimT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((phidriftLim(1) /= SpecieT(s)%FluxTubeT(f)%phidriftLimT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%phidriftLimT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%phidriftLimT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phidriftLimT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION', &
						' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET REFERENCE PARALLEL POTENTIAL DROPS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

					! ----------------------------------------------------

					EPar0p(1)= PhiPar0p/dPhiPar0p ! [V/m]

					! Create nested derived data types
					SpecieT(s)%FluxTubeT(f)%PhiPar0BT= EPar0p(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((EPar0p(1) /= SpecieT(s)%FluxTubeT(f)%PhiPar0BT(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%PhiPar0BT(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%PhiPar0BT(1))) &
						.eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' PhiPar0BT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', AND FLUX TUBE= ', f, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if
			end do
		end do

		! ----------------------------------------------------

		! SET ICR HEATING PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) then
					do Qind= 1, 1, 1

						! ----------------------------------------------------

						OmegaG0p(1)= 2d0*pi*f0p ! Ref. gyrofreq. corresponding to EPerp0 [rads/s]

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%lambdaPerppT= lambdaPerpp
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EtaLHpT= EtaLHp
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp1pT= XiPerp1p
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp2pT= XiPerp2p
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%S0pT= S0p
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%OmegaG0pT= OmegaG0p
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp1pT= ChiPerp1p
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp2pT= ChiPerp2p

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((lambdaPerpp /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%lambdaPerppT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%lambdaPerppT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%lambdaPerppT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' lambdaPerppT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((EtaLHp /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EtaLHpT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EtaLHpT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%EtaLHpT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EtaLHpT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((XiPerp1p /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp1pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp1pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp1pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp1pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((XiPerp2p /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp2pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp2pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%XiPerp2pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp2pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((S0p /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%S0pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%S0pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%S0pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' S0pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((OmegaG0p(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%OmegaG0pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%OmegaG0pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%OmegaG0pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaG0pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((ChiPerp1p /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp1pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp1pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp1pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp1pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((ChiPerp2p /= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp2pT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp2pT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ChiPerp2pT(1))) &
							.eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp2pT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) '----------------------------------------------------'
					write(*, *) 'Kinetic model of Auroral ion OutflowS (KAOS)'
					write(*, *) 'Robert M. Albarran II, Ph.D. contact: albarran1@atmos.ucla.edu'
					write(*, *) 'Department of Atmospheric and Oceanic Sciences, University of California, Los Angeles'

		      call date_and_time(datechar(1), datechar(2), datechar(3), dateint)
					write(monthstring, '(i10)') dateint(2)
					write(daystring, '(i10)') dateint(3)
					write(yearstring, '(i10)') dateint(1)
					write(hourstring, '(i10)') dateint(5)
					write(minutestring, '(i10)') dateint(6)
					write(secondstring, '(i10)') dateint(7)
	 				write(*, *) trim('Month= ' // adjustl(monthstring))
					write(*, *) trim('Day= ' // adjustl(daystring))
					write(*, *) trim('Year= ' // adjustl(yearstring))
					write(*, *) trim('Hour= ' // adjustl(hourstring))
					write(*, *) trim('Minute= ' // adjustl(minutestring))
					write(*, *) trim('Second= ' // adjustl(secondstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					if (SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1) then
						write(*, *) trim('KAOS SPIN-UP')
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'SIMULATION PARAMETERS:'
					write(*, *)

					write(paramstring, '(i10)') ranksize(1)
					write(*, *) trim('Number of MPI ranks= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QCell0T(1)%qGL0T(1) <= 0d0) then
						write(*, *) 'Northern Magnetic Hemisphere'
					else if (SpecieT(s)%FluxTubeT(f)%QCell0T(1)%qGL0T(1) > 0d0) then
						write(*, *) 'Southern Magnetic Hemisphere'
					end if

					write(*, *) trim('Data Export Path= ' // adjustl(dataexportdir))

					if ((SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 1) .or. &
						(SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1)) then
						write(*, *) trim('Spin-up Data I/O Path= ' // adjustl(Densitydatadir))
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'CONFIGURATION-SPACE PARAMETERS:'
					write(*, *)

					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T(1)%pGC0T
					write(*, *) trim('L-shell [RE]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)% &
						QCell0T(2)%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Lower Grid Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T &
					(((NqUB(1)- NqLB(1))+ 3))%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Upper Grid Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)% &
						QCell0T(SpecieT(s)%FluxTubeT(f)%NqICAT(1)+ 1)%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Lower Ion Injection Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T &
					(SpecieT(s)%FluxTubeT(f)%NqICBT(1)+ 1)%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Upper Ion Initialization Altitude [km]= ' // adjustl(paramstring))

					write(paramstring, '(i10)') Stot
					write(*, *) trim('Total Number of Particle Species= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(1)%NfT(1)
					write(*, *) trim('Total Number of Flux Tubes= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqGT(1)
					write(*, *) trim('Total Number of Ion/ENA Config-Space Grid Cells= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION VELOCITY-SPACE PARAMETERS:'
					write(*, *)

					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp1GLT(1)*1d-3
						write(*, *) trim('Lower Vperp1 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp2GLT(1)*1d-3
						write(*, *) trim('Lower Vperp2 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1)*1d-3
						write(*, *) trim('Upper Vperp1 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1)*1d-3
						write(*, *) trim('Upper Vperp2 Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%VparGLT(1)*1d-3
						write(*, *) trim('Lower Vpar Limit [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1)*1d-3
						write(*, *) trim('Upper Vpar Limit [km/s]= ' // adjustl(paramstring))

						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)
						write(*, *) trim('Total Number of Ion Vperp1 Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)
						write(*, *) trim('Total Number of Ion Vperp2 Grid Cells= ' // adjustl(paramstring))
					else
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)
						write(*, *) trim('Total Number of Ion Vperp Grid Cells= ' // adjustl(paramstring))
					end if
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)
					write(*, *) trim('Total Number of Ion Vpar Grid Cells= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12sigma*1d-3
					write(*, *) trim('Vperp1 and Vperp2 Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12sigmaFac
					write(*, *) trim('Number of Linear Vperp1 and Vperp2 Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(i10)') Vperp12NlinRange
					write(*, *) trim('Number of Vperp1 and Vperp2 Linear Grid Cells Spanning to &
						Last Standard Deviation= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vparsigma*1d-3
					write(*, *) trim('Vpar Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VparsigmaFac
					write(*, *) trim('Number of Linear Vpar Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(i10)') VparNlinRange
					write(*, *) trim('Number of Vpar Linear Grid Cells Spanning to &
						Last Standard Deviation= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ENA VELOCITY-SPACE PARAMETERS:'
						write(*, *)

						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)
						write(*, *) trim('Total Number of ENA Vp Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)
						write(*, *) trim('Total Number of ENA Vq Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)
						write(*, *) trim('Total Number of ENA Vphi Grid Cells= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') Vpphisigma*1d-3
						write(*, *) trim('Vp and Vphi Standard Deviation [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphisigmaFac
						write(*, *) trim('Number of Linear Vp and Vpphi Standard Deviations &
							Spanned by Linear Grid= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphiNlinRange
						write(*, *) trim('Number of Vp and Vphi Linear Grid Cells Spanning to &
							Last Standard Deviation= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') Vqsigma*1d-3
						write(*, *) trim('Vq Standard Deviation [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VqsigmaFac
						write(*, *) trim('Number of Linear Vq Standard Deviations &
							Spanned by Linear Grid= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VqNlinRange
						write(*, *) trim('Number of Vq Linear Grid Cells Spanning to &
							Last Standard Deviation= ' // adjustl(paramstring))
					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'TIME PARAMETERS:'
					write(*, *)

					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%AT(1)
					write(*, *) trim('Beginning Simulation Time [s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%BT(1)
					write(*, *) trim('End Simulation Time [s]= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NtT(1)
					write(*, *) trim('Total Number of Base Time-Steps= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%hT(1)
					write(*, *) trim('Base Time-Step Duration [s]= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)
					write(*, *) trim('Total Number of Injection Time-Steps= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NNtT(1)
					write(*, *) trim('Total Number of Statistical Time-Steps= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)*SpecieT(s)%FluxTubeT(f)%hT(1)
					write(*, *) trim('Injection Time-Step Duration [s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ndatfacT(1)*SpecieT(s)%FluxTubeT(f)%hT(1)
					write(*, *) trim('Statistical Time-Step Duration [s]= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION SIMULATION FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%IONNOISEflagT(1)
					write(*, *) trim('Ion Statistical Noise Reduction Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1)
					write(*, *) trim('Ion Fluid Data Export Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'BOUNDARY CONDITIONS FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1)
					write(*, *) trim('Lower Boundary Condition Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1)
					write(*, *) trim('Upper Boundary Condition Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1)
					write(*, *) trim('Lower Boundary Density Replenish Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1)
					write(*, *) trim('Upper Boundary Density Replenish Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%DENSITYPROFILEflagT(1)
					write(*, *) trim('Initial Density Profile Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%STATICINJECTIONflagT(1)
					write(*, *) trim('Static Injection Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1)
					write(*, *) trim('Density Output Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1)
					write(*, *) trim('Density Input Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION FORCE FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ICRflagT(1)
					write(*, *) trim('Ion Cyclotron Wave Heating Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%MIRRORflagT(1)
					write(*, *) trim('Mirror Force Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%GRAVflagT(1)
					write(*, *) trim('Gravitational Force Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBflagT(1)
					write(*, *) trim('Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBSELFCONSISTflagT(1)
					write(*, *) trim('Self-Consistent Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAMBSIGNflagT(1)
					write(*, *) trim('Self-Consistent Sign of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAINERTIALflagT(1)
					write(*, *) trim('Inertial Term of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EAPRESSUREflagT(1)
					write(*, *) trim('Pressure Term of Ambipolar Electric Field Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%EPARflagT(1)
					write(*, *) trim('Parallel Electric Field Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION MOMENT FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1)
					write(*, *) trim('Moment Filter Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1)
					write(*, *) trim('2D Perpendicular Velocity Space Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1)
					write(*, *) trim('Ion Phase-Space Distribution Function Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Density Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Velocity Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Velocity Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Total Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Thermal Energy Moment Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTCENTERflagT(1)
					write(*, *) trim('Ion Phase-Space Perpendicular Thermal Energy Moment Centering Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTCENTERflagT(1)
					write(*, *) trim('Ion Phase-Space Parallel Thermal Energy Moment Centering Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1)
					write(*, *) trim('Ion Reference Moments Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ENA FLAGS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%ENANOISEflagT(1)
					write(*, *) trim('ENA Statistical Noise Reduction Flag= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1)
					write(*, *) trim('Charge-Exchange Flag= ' // adjustl(paramstring))

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION INITIALIZATION PARAMETERS:'
					write(*, *)

					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqICAT(1)
					write(*, *) trim('Ion Initialization LB Grid Cell= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NqICBT(1)
					write(*, *) trim('Ion Initialization UB Grid Cell= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T( &
						SpecieT(s)%FluxTubeT(f)%NqICAT(1))%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Ion Initialization Lower Altitude Boundary [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T( &
						SpecieT(s)%FluxTubeT(f)%NqICBT(1))%rGC0T(1)*1d-3- RE*1d-3
					write(*, *) trim('Ion Initialization Upper Altitude Boundary [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%zns0T(1)- RE)*1d-3
					write(*, *) trim('Ion Initialization Reference Altitude [km]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ns0T(1)
					write(*, *) trim('Ion Initialization Reference Density [m^-3]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T(1)%Ts0T(1)
					write(*, *) trim('Ion Initialization Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T(1)%TsPerp0T(1)
					write(*, *) trim('Ion Initialization Perpendicular Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCell0T(1)%TsPar0T(1)
					write(*, *) trim('Ion Initialization Parallel Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%Te0T(1)
					write(*, *) trim('Electron Initialization Temperature [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') dNTe
					write(*, *) trim('Additive Increment of Electron Temperature on Statistical Time-steps [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') dNTeEND
					write(*, *) trim('Additive Increment Cap of Electron Temperature on Statistical Time-steps [K]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') dns0
					write(*, *) trim('Multiplicative Factor of Reference Ion Density= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%Qindns0T(1)
					write(*, *) trim('Reference Density Grid Cell= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'NEUTRAL OXYGEN INITIALIZATION PARAMETERS:'
						write(*, *)

						write(paramstring, '(D10.2)') (SpecieT(s)%FluxTubeT(f)%zns0neutT(1)- RE)*1d-3
						write(*, *) trim('Neutral Oxygen Initialization Reference Altitude [km]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ns0neutT(1)
						write(*, *) trim('Neutral Oxygen Initialization Reference Density [m^-3]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') TNeut
						write(*, *) trim('Neutral Oxygen Initialization Temperature [K]= ' // adjustl(paramstring))

					end if

					if (SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ICR HEATING PARAMETERS:'
						write(*, *)

						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%lambdaPerppT(1)
						write(*, *) trim('BBELF Wavelength [m]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%EtaLHpT(1)
						write(*, *) trim('BBELF Wave Power LHP Fraction= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%XiPerp1pT(1)
						write(*, *) trim('BBELF Wave Power Fraction Along Vperp1= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%XiPerp2pT(1)
						write(*, *) trim('BBELF Wave Power Fraction Along Vperp2= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%S0pT(1)
						write(*, *) trim('Wave Spectral Energy Density [(V^2/m^2)/Hz]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%OmegaG0pT(1)/(2d0*pi)
						write(*, *) trim('Reference Ion Cyclotron Frequency [Hz]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%ChiPerp1pT(1)
						write(*, *) trim('Wave Spectral Index Along Vperp1= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%QCellT(1)%ChiPerp2pT(1)
						write(*, *) trim('Wave Spectral Index Along Vperp2= ' // adjustl(paramstring))

					end if

					if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'PARALLEL POTENTIAL PARAMETERS:'
						write(*, *)

						write(paramstring, '(i10)') Eparlim
						write(*, *) trim('Computational Time-step of Activated Potential Drop= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') PhiPar0p
						write(*, *) trim('Reference Parallel Potential Drop [V]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') dPhiPar0p
						write(*, *) trim('Reference Parallel Potential Distance [m]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%PhiPar0BT
						write(*, *) trim('Reference Parallel Electric Field [V/m]= ' // adjustl(paramstring))

					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *) 'ION LIMITS:'
					write(*, *)

					if (SpecieT(1)%FluxTubeT(1)%MOMENTFILTERflagT(1) == 1) then
						write(paramstring, '(i10)') M0MAfilterPt
						write(*, *) trim('Density Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1Perp1MAfilterPt
						write(*, *) trim('Perp1 Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1Perp2MAfilterPt
						write(*, *) trim('Perp2 Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M1ParMAfilterPt
						write(*, *) trim('Parallel Velocity Moment Moving Average Point= ' // adjustl(paramstring))
						write(paramstring, '(i10)') M2ParMAfilterPt
						write(*, *) trim('Parallel Energy Moment Moving Average Point= ' // adjustl(paramstring))
					end if
					write(paramstring, '(D10.2)') IonNoiseLimitNph
					write(*, *) trim('Minimum Phase-Space Ion Macroparticle Number= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%IonNoiseLimitT(1)
					write(*, *) trim('Ion Noise Limit= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%pdriftLimT(1)
					write(*, *) trim('Ion Iterative L-Shell Drift Limit [RE]= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

						write(*, *)
						write(*, *) '----------------------------------------------------'
						write(*, *) 'ENA LIMITS:'
						write(*, *)

						write(paramstring, '(D10.2)') ENANoiseLimitNph
						write(*, *) trim('Minimum Phase-Space ENA Macroparticle Number= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%ENANoiseLimitT(1)
						write(*, *) trim('ENA Noise Limit= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%pdriftLimT(1)
						write(*, *) trim('ENA Total R Drift Limit [m]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%thetadriftLimT(1)
						write(*, *) trim('ENA Total Theta Drift Limit [rads]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%phidriftLimT(1)
						write(*, *) trim('ENA Total Phi Drift Limit [rads]= ' // adjustl(paramstring))

					end if

					write(*, *)
					write(*, *) '----------------------------------------------------'
					write(*, *)

				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S1End)
					write(S1string, '(i10)')  nint(S1End)
					write(*, *) trim('%% 1- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S1string)) // &
						trim(' s. SIMULATION PARAMETERIZATION COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine SimParameterizationSub

end module SimParameterization
