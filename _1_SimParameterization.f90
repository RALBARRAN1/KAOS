module SimParameterization

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	1- SIMULATION PARAMETERIZATION:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

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

		allocate(SpecieT(Stot))

		! ----------------------------------------------------

		! IMPORT NUMBER OF FLUX TUBES PER PARTICLE SPECIES FROM MATLAB:

		! ----------------------------------------------------

		allocate(Nfp(Stot))

		! ----------------------------------------------------

		do s= 1, Stot, 1
			write(sstring, '(I5)') s

			! ----------------------------------------------------

			open (unit= 0, file= dataimportdir // adjustl(adjustr(sstring) // '_' // 'Nfmat.bin'), &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) Nfp(s)
			close(0)

			Nf(1)= Nfp(s)

			SpecieT(s)%NfT= Nf ! Create nested derived data types

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((Nf(1) /= SpecieT(s)%NfT(1)) .or. (size(SpecieT(s)%NfT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%NfT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NfT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		deallocate(Nfp)

		! ----------------------------------------------------

		! ALLOCATE FLUX TUBE DERIVED DATA TYPE NESTED IN PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1

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

					! For VISIONS-1 CS:
					!NqICA(1)= 1d0
					!NqICB(1)= 36d0

					! For (Barakat '83) VS01:
					!NqICA(1)= 1d0
					!NqICB(1)= 26d0

					! For (Barakat '83) VS02:
					!NqICA(1)= 1d0
					!NqICB(1)= 36d0

					! For (Barghouthi '98) VS03 VS04:
					!NqICA(1)= 1d0
					!NqICB(1)= 26d0

					! For (Barghouthi '97 case1 and Crew '90, Barghouthi '97 case2) VS1 VS2:
					!NqICA(1)= 1d0
					!NqICB(1)= 15d0

					! For (Barghouthi '94) VS3:
					!NqICA(1)= 1d0
					!NqICB(1)= 26d0

					! For (Wilson '90) VS4:
					!NqICA(1)= 1d0
					!NqICB(1)= 26d0

					! For (Wilson '92 case1) VS5:
					!NqICA(1)= 1d0
					!NqICB(1)= 16d0

					! For (Wilson '92 case2) VS6:
					!NqICA(1)= 1d0
					!NqICB(1)= 12d0

					! For (Wu '02 and Brown '95) VS8 VS9 VS10 VS11 VS12:
					!NqICA(1)= 1d0
					!NqICB(1)= 31d0

					! For (Wu '02 and Brown '95) VS8 VS9 VS10 VS11 VS12:
					!NqICA(1)= 1d0
					!NqICB(1)= 31d0

					! For Parametric Study:
					!NqICA(1)= 1d0 ! T1
					!NqICB(1)= 2d0 ! 19d0 ! T1

					NqLB(1)= NqICA
					NqUB(1)= NqICB
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

				! For (Barakat '83 case1 and case2) VS01 VS02:
				!(Ti= 3000 K, Te= 3000 K) and (Ti= 3000 K, Te= 10000 K):
				!zns0(1)= RE+ 10870d3 ! Reference O+ r value [m]
				!ns0(1)= 5d7 ! Reference O+ number density at zns0 [m^-3]

				! For (Barghouthi '98 case1 and case2) VS03 VS04:
				!(Ti= 3000 K, Te= 3000 K):
				!zns0(1)= RE+ 10870d3 ! Reference O+ r value [m]
				!ns0(1)= 1d8 ! Reference O+ number density at zns0 [m^-3]

				! For (Barghouthi '97 case1 and Crew '90, Barghouthi '97 case2) VS1 VS2:
				!(Ti= 2321 K, Te= 1000 K):
				!zns0(1)= RE+ 7656d3 ! Reference O+ r value [m]
				!ns0(1)= 5d9 ! Reference O+ number density at zns0 [m^-3]

				! For (Barghouthi '94) VS3:
				!(Ti= 3000 K, Te= 3000 K):
				!zns0(1)= RE+ 10846d3 ! Reference O+ r value [m]
				!ns0(1)= 1d8 ! Reference O+ number density at zns0 [m^-3]

				! For (Wilson '90) VS4:
				!(Ti= 3000 K, Te= 3000 K):
				!zns0(1)= RE+ 10846d3 ! Reference O+ r value [m]
				!ns0(1)= 5d7 ! Reference O+ number density at zns0 [m^-3]

				! For (Wilson '92 case1) VS5:
				!(Ti= 3000 K, Te= 3000 K):
				!zns0(1)= RE+ 1000d3 ! Reference O+ r value [m]
				!ns0(1)= 5d10 ! Reference O+ number density at zns0 [m^-3]

				! For (Wilson '92 case2) VS6:
				!(Ti= 2000 K, Te= 2000 K):
				!zns0(1)= RE+ 1000d3 ! Reference O+ r value [m]
				!ns0(1)= 1d10 ! Reference O+ number density at zns0 [m^-3]

				! For (Wu '02 and Brown '95) VS8 VS9 VS10 VS11 VS12:
				!(Ti= 8000 K, Te= 8000 K):
				!zns0(1)= RE+ 3*RE ! Reference O+ r value [m]
				!ns0(1)= 1d6 ! Reference O+ number density at zns0 [m^-3]

				if (s == 1) then
					!nsnormfac(1)= 6.2d15 ! Macro-particle density norm. const. (inversely proportional to density)
					!nsnormfac(1)= 4.05d12 ! For (Barakat '83 case1) VS01
					!nsnormfac(1)= 9.95d12 ! For (Barakat '83 case2) VS02
					!nsnormfac(1)= 8.10d12 ! For (Barghouthi '98 case1 and case2) VS03 VS04
					!nsnormfac(1)= 2.5d14 ! For (Barghouthi '97 case1 and Crew '90, Barghouthi '97 case2) VS1 VS2
					!nsnormfac(1)= 3.21d14 ! For (Barghouthi '97 case1 and Crew '90, Barghouthi '97 case2) VS1 VS2
					!nsnormfac(1)= 8.1d12 ! For (Barghouthi '94) VS3
					!nsnormfac(1)= 4.05d12 ! For (Wilson '90) VS4
					!nsnormfac(1)= 6.37d14 ! For (Wilson '92 case1) VS5
					!nsnormfac(1)= 9.85d13 ! For (Wilson '92 case2) VS6
					!nsnormfac(1)= 2.22d10 ! For (Wu '02 and Brown '95) VS8 VS9 VS10 VS11 VS12
				end if

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

		! IMPORT ALL CONFIGURATION-SPACE GRID DIMENSIONS FROM MATLAB:

		do s= 1, Stot, 1
			write(sstring, '(I5)') s

			! ----------------------------------------------------

			allocate(NqGp(Stot, SpecieT(s)%NfT(1)), HiCratioMeanp(Stot, SpecieT(s)%NfT(1)), &
				Tep(Stot, SpecieT(s)%NfT(1)))

			! ----------------------------------------------------

			do f= 1, SpecieT(s)%NfT(1), 1
				write(fstring, '(I5)') f

				! ----------------------------------------------------

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'NqGmat.bin')), status= 'old', form= 'unformatted', access= 'stream')
				read(0) NqGp(s, f)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'HiCratioMeanmat.bin')), status= 'old', form= 'unformatted', access= 'stream')
				read(0) HiCratioMeanp(s, f)
				close(0)

				open (unit= expint, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'Temat.bin')), status= 'old', form= 'unformatted', access= 'stream')
				read(expint) Tep(s, f)
				close(expint)

				NqG(1)= NqGp(s, f)
				HiCratioMean(1)= HiCratioMeanp(s, f)
				Te(1)= Tep(s, f)

				SpecieT(s)%FluxTubeT(f)%NqGT= NqG ! Create nested derived data types
				SpecieT(s)%FluxTubeT(f)%HiCratioMeanT= HiCratioMean
				SpecieT(s)%FluxTubeT(f)%TeT= Te

				! Re-index to account for lower and upper ghost cells (non-computational domain)
				SpecieT(s)%FluxTubeT(f)%NqG0T(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)
				SpecieT(s)%FluxTubeT(f)%NqGT(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)- 2d0

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((NqG(1) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((NqG(1) /= SpecieT(s)%FluxTubeT(f)%NqGT(1)+ 2d0) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((HiCratioMean(1) /= SpecieT(s)%FluxTubeT(f)%HiCratioMeanT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%HiCratioMeanT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%HiCratioMeanT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' HiCratioMeanT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((Te(1) /= SpecieT(s)%FluxTubeT(f)%TeT(1)) .or. &
					(size(SpecieT(s)%FluxTubeT(f)%TeT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%TeT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TeT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(s)%FluxTubeT(f)%TeT(1) > dNTeEND) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INITIAL TeT IS GREATER THAN Te CAP FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN SIMULATION PARAMETERIZATION', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if
				! ----------------------------------------------------

			end do

			! ----------------------------------------------------

			deallocate(NqGp, HiCratioMeanp, Tep)

			! ----------------------------------------------------

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

		! IMPORT PARTICLE MASS, CHARGE, AND ATOMIC RADIUS PER PARTICLE SPECIES FROM MATLAB:

		allocate(Qindns0p(Stot), msp(Stot), qsp(Stot))

		do s= 1, Stot, 1
			write(sstring, '(I5)') s

			! ----------------------------------------------------

			open (unit= 0, file= dataimportdir &
				// adjustl(adjustr(sstring) // '_' // 'Qindns0mat.bin'), &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) Qindns0p(s)
			close(0)

			open (unit= 0, file= dataimportdir &
				// adjustl(adjustr(sstring) // '_' // 'msmat.bin'), &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) msp(s)
			close(0)

			open (unit= 0, file= dataimportdir &
				// adjustl(adjustr(sstring) // '_' // 'qsmat.bin'), &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) qsp(s)
			close(0)

			Qindns0(1)= Qindns0p(s)
			ms(1)= msp(s)
			qs(1)= qsp(s)

			SpecieT(s)%Qindns0T= Qindns0 ! Create nested derived data types
			SpecieT(s)%msT= ms
			SpecieT(s)%qsT= qs

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((Qindns0(1) /= SpecieT(s)%Qindns0T(1)) .or. (size(SpecieT(s)%Qindns0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%Qindns0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0T HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((ms(1) /= SpecieT(s)%msT(1)) .or. (size(SpecieT(s)%msT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%msT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' msT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((qs(1) /= SpecieT(s)%qsT(1)) .or. (size(SpecieT(s)%qsT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%qsT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qsT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' IN SIMULATION PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		deallocate(Qindns0p, msp, qsp)

		! ----------------------------------------------------

		! IMPORT ALL PARTICLE SPECIES INITIAL TEMPERATURES AND CONFIGURATION-SPACE GRID
 		! PARAMETERS FROM MATLAB:

		do s= 1, Stot, 1
			write(sstring, '(I5)') s
			do f= 1, SpecieT(s)%NfT(1), 1
				write(fstring, '(I5)') f

				! ----------------------------------------------------

				allocate(qGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					hqCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					dpCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					dqCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					dphiCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					TsPerpp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					TsParp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					rGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					phiGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					thetaGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					ellGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					qGLp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					qGHp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					pGCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					d3xCp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 3)), &
					SpecieT(s)%FluxTubeT(f)%rGCTp(((NqUB(1)- NqLB(1))+ 3)), &
					SpecieT(s)%FluxTubeT(f)%hqCTp(((NqUB(1)- NqLB(1))+ 3)), &
					SpecieT(s)%FluxTubeT(f)%dqCTp(((NqUB(1)- NqLB(1))+ 3)))

				! ----------------------------------------------------

				! COMPUTE TOTAL FIELD-LINE ARC LENGTH:

				allocate(hq(Stot, SpecieT(s)%NfT(1), SpecieT(s)%FluxTubeT(f)%NqG0T(1)), &
					dq(Stot, SpecieT(s)%NfT(1), SpecieT(s)%FluxTubeT(f)%NqG0T(1)))
				allocate(SpecieT(s)%FluxTubeT(f)%ellqCT(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'hqCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) hq(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'dqCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) dq(s, f, :)
				close(0)

				SpecieT(s)%FluxTubeT(f)%ellqCT(:)= hq(s, f, :)*dq(s, f, :)
				SpecieT(s)%FluxTubeT(f)%SUMellqCT(1)= sum(SpecieT(s)%FluxTubeT(f)%ellqCT(:))

				! ----------------------------------------------------

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'qGCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) qGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'hqCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) hqCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'dpCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) dpCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'dqCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) dqCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'dphiCmat.bin')), &
					status= 'old', form= 'unformatted', access= 'stream')
				read(0) dphiCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'TsPerpmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) TsPerpp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'TsParmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) TsParp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'rGCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) rGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'phiGCmat.bin')), &
					status= 'old', form= 'unformatted', access= 'stream')
				read(0) phiGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'thetaGCmat.bin')), &
					status= 'old', form= 'unformatted', access= 'stream')
				read(0) thetaGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'ellGCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) ellGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'qGLmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) qGLp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'qGHmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) qGHp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'pGCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) pGCp(s, f, :)
				close(0)

				open (unit= 0, file= dataimportdir &
					// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
					// '_' // 'd3xCmat.bin')), status= 'old', &
					form= 'unformatted', access= 'stream')
				read(0) d3xCp(s, f, :)
				close(0)

				do Qind= NqLB(1), NqUB(1)+ 2, 1
					write(Qindstring, '(I5)') Qind

					! ----------------------------------------------------

					qGC(1)= qGCp(s, f, Qind)
					hqC(1)= hqCp(s, f, Qind)
					dpC(1)= dpCp(s, f, Qind)
					dqC(1)= dqCp(s, f, Qind)
					dphiC(1)= dphiCp(s, f, Qind)
					TsPerp(1)= TsPerpp(s, f, Qind)
					TsPar(1)= TsParp(s, f, Qind)
					rGC(1)= rGCp(s, f, Qind)
					phiGC(1)= phiGCp(s, f, Qind)
					thetaGC(1)= thetaGCp(s, f, Qind)
					ellGC(1)= ellGCp(s, f, Qind)
					qGL(1)= qGLp(s, f, Qind)
					qGH(1)= qGHp(s, f, Qind)
					pGC(1)= pGCp(s, f, Qind)
					d3xC(1)= d3xCp(s, f, Qind)

					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T= qGC ! Create nested derived data types (including ghost cells)
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T= hqC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T= dpC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T= dqC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T= dphiC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T= rGC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T= phiGC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T= thetaGC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T= ellGC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T= qGL
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T= qGH
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T= pGC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T= d3xC
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T= TsPerp
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T= TsPar

					! Total initial ion temperature [K]
					! Note: Isotropic temperature from grid data (Ti= Tperp= Tpar= Tperp1= Tperp2)
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%Ts0T= (1d0/3d0)*TsPar+ (2d0/3d0)*TsPerp

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((qGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((hqC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((dpC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((dqC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((dphiC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((TsPerp(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerpT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((TsPar(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsParT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((rGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((phiGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((thetaGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((ellGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((qGL(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGLT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((qGH(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGHT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((pGC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGCT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((d3xC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1) .or. &
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

				deallocate(qGCp, hqCp, dpCp, dqCp, dphiCp, TsPerpp, TsParp, &
					rGCp, phiGCp, thetaGCp, ellGCp, qGLp, qGHp, pGCp, d3xCp)

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SET SIMULATION DURATION AND TIME-STEP FOR MOMENT COMPUTATION:

		i= (0, 1) ! Define sqrt(-1)

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

		! IMPORT ION GRID PARAMETERS FROM MATLAB:

		open (unit= 0, file= dataimportdir // 'Vperp12NlinRangemat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) Vperp12NlinRange
		close(0)
		open (unit= 0, file= dataimportdir // 'VparNlinRangemat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) VparNlinRange
		close(0)
		open (unit= 0, file= dataimportdir // 'Vperp12Gridlinspacemat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) Vperp12Gridlinspace
		close(0)
		open (unit= 0, file= dataimportdir // 'Vperp12sigmamat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) Vperp12sigma
		close(0)
		open (unit= 0, file= dataimportdir // 'Vperp12sigmaFacmat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) Vperp12sigmaFac
		close(0)
		open (unit= 0, file= dataimportdir // 'VparGridlinspacemat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) VparGridlinspace
		close(0)
		open (unit= 0, file= dataimportdir // 'Vparsigmamat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) Vparsigma
		close(0)
		open (unit= 0, file= dataimportdir // 'VparsigmaFacmat.bin', &
			status= 'old', form= 'unformatted', access= 'stream')
		read(0) VparsigmaFac
		close(0)

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

		! IMPORT ENA MASS AND INITIAL TEMPERATURE FROM MATLAB:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			open (unit= 0, file= dataimportdir // 'mNeutmat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) mNeut
			close(0)
			open (unit= 0, file= dataimportdir // 'TNeutmat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) TNeut
			close(0)
			open (unit= 0, file= dataimportdir // 'VpphiNlinRangemat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VpphiNlinRange
			close(0)
			open (unit= 0, file= dataimportdir // 'VqNlinRangemat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VqNlinRange
			close(0)
			open (unit= 0, file= dataimportdir // 'VpphiGridlinspacemat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VpphiGridlinspace
			close(0)
			open (unit= 0, file= dataimportdir // 'Vpphisigmamat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) Vpphisigma
			close(0)
			open (unit= 0, file= dataimportdir // 'VpphisigmaFacmat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VpphisigmaFac
			close(0)
			open (unit= 0, file= dataimportdir // 'VqGridlinspacemat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VqGridlinspace
			close(0)
			open (unit= 0, file= dataimportdir // 'Vqsigmamat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) Vqsigma
			close(0)
			open (unit= 0, file= dataimportdir // 'VqsigmaFacmat.bin', &
				status= 'old', form= 'unformatted', access= 'stream')
			read(0) VqsigmaFac
			close(0)
		end if

		! ----------------------------------------------------

		! IMPORT ALL ION VELOCITY-SPACE GRID DIMENSIONS FROM MATLAB:

		do s= 1, Stot, 1
			write(sstring, '(I5)') s
			do f= 1, SpecieT(s)%NfT(1), 1
				write(fstring, '(I5)') f

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					! ----------------------------------------------------

					allocate(NVperp1Gp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)), &
						NVperp2Gp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)))

					! ----------------------------------------------------

					do Qind= 1, 1, 1
						! write(Qindstring, '(I5)') Qind ! Note uncomment for config grid cell dependendent velocity grid
						write(Qindstring, '(I5)') 1

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVperp1Gmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVperp1Gp(s, f, Qind)
						close(0)

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVperp2Gmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVperp2Gp(s, f, Qind)
						close(0)

						NVperp1G(1)= NVperp1Gp(s, f, Qind)
						NVperp2G(1)= NVperp2Gp(s, f, Qind)

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT= NVperp1G
						SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT= NVperp2G

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((NVperp1G(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp1GT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((NVperp2G(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp2GT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					deallocate(NVperp1Gp, NVperp2Gp)

					! ----------------------------------------------------

				else

					! ----------------------------------------------------

					allocate(NVperpGp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)))

					! ----------------------------------------------------

					do Qind= 1, 1, 1
						! write(Qindstring, '(I5)') Qind ! Note uncomment for config grid cell dependendent velocity grid
						write(Qindstring, '(I5)') 1

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVperpGmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVperpGp(s, f, Qind)
						close(0)

						NVperpG(1)= NVperpGp(s, f, Qind)

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT= NVperpG

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((NVperpG(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperpGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					deallocate(NVperpGp)

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

				allocate(NVparGp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)))

				! ----------------------------------------------------

				do Qind= 1, 1, 1
					! write(Qindstring, '(I5)') Qind ! Note uncomment for config grid cell dependendent velocity grid
					write(Qindstring, '(I5)') 1

					! ----------------------------------------------------

					open (unit= 0, file= dataimportdir &
						// adjustl(adjustr(sstring) // '_' // adjustl(adjustr(fstring) &
						// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVparGmat.bin'))), &
						status= 'old', form= 'unformatted', access= 'stream')
					read(0) NVparGp(s, f, Qind)
					close(0)

					NVparG(1)= NVparGp(s, f, Qind)

					! Create nested derived data types
					SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT= NVparG

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					if ((NVparG(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVparGT HAS', &
							' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
							' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do

				! ----------------------------------------------------

				deallocate(NVparGp)

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! IMPORT ALL ENA VELOCITY-SPACE GRID DIMENSIONS FROM MATLAB:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do s= 1, Stot, 1
				write(sstring, '(I5)') s
				do f= 1, SpecieT(s)%NfT(1), 1
					write(fstring, '(I5)') f

					! ----------------------------------------------------

					allocate(NVpGp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)), &
						NVqGp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)), &
						NVphiGp(Stot, SpecieT(s)%NfT(1), ((NqUB(1)- NqLB(1))+ 1)))

					! ----------------------------------------------------

					do Qind= 1, 1, 1
						! write(Qindstring, '(I5)') Qind ! Note uncomment for config grid cell dependendent velocity grid
						write(Qindstring, '(I5)') 1

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVpGmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVpGp(s, f, Qind)
						close(0)

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVqGmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVqGp(s, f, Qind)
						close(0)

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) &
							// '_' // adjustl(adjustr(Qindstring) // '_' // 'NVphiGmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) NVphiGp(s, f, Qind)
						close(0)

						NVpG(1)= NVpGp(s, f, Qind)
						NVqG(1)= NVqGp(s, f, Qind)
						NVphiG(1)= NVphiGp(s, f, Qind)

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT= NVpG
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT= NVqG
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT= NVphiG

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((NVpG(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVpGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((NVqG(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVqGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((NVphiG(1) /= SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVphiGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN SIMULATION', &
								' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					deallocate(NVpGp, NVqGp, NVphiGp)

					! ----------------------------------------------------

				end do
			end do
		end if

		! ----------------------------------------------------

		! ALLOCATE VELOCITY-SPACE DERIVED DATA TYPES NESTED IN CONFIGURATION-SPACE
		! TYPE NESTED IN FLUX TUBE TYPES NESTED IN PARTICLE SPECIES TYPE:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do Qind= NqLB(1), NqUB(1), 1

					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						! Allocate V2PerpCellT(Vperp1ind, Vperp2ind, Vparind) derived data type nested in
						! FluxTubeT(f) nested in SpecieT(s)
						allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))
					else
						! Allocate VCellT(Vperpind, Vparind) derived data type nested in
						! FluxTubeT(f) nested in SpecieT(s)
						allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))
					end if

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						! Allocate V3CellT(Vpind, Vqind, Vphiind) derived data type nested in
						! FluxTubeT(f) nested in SpecieT(s)
						allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)))
					end if

				end do
			end do
		end do

		! ----------------------------------------------------

		! COMPUTE NUMBER OF MPI RANKS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				call MPI_COMM_SIZE(MPI_COMM_WORLD, ranksize(1), ierr)
			end do
		end do

		allocate(RankT(ranksize(1)))

		! ----------------------------------------------------

		! IMPORT ALL ION EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES FROM MATLAB:

		do s= 1, Stot, 1
			write(sstring, '(I5)') s
			do f= 1, SpecieT(s)%NfT(1), 1
				write(fstring, '(I5)') f
				do Qind= 1, 1, 1
					write(Qindstring, '(I5)') 1

					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

						! ----------------------------------------------------

						allocate(Vperp1GLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp1GHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp1GCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperp1Gpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperp2Gpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGL2pp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGH2pp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGC2pp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVparG2pp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							d3vC2pp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))

						allocate(Vperp1GLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp1GHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp1GCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperp1Gp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							Vperp2GCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperp2Gp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGL2p(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGH2p(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGC2p(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVparG2p(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							d3vC2p(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp1GLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp1GLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp1GHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp1GHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp1GCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp1GCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVperp1Cmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVperp1Gpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp2GLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp2GLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp2GHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp2GHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'Vperp2GCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) Vperp2GCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVperp2Cmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVperp2Gpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGL2pp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGH2pp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGC2pp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVparCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVparG2pp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'd3vCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) d3vC2pp(:)
						close(0)

						Vperp1GLp(:, :, :)= reshape(Vperp1GLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						Vperp1GHp(:, :, :)= reshape(Vperp1GHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						Vperp1GCp(:, :, :)= reshape(Vperp1GCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						dVperp1Gp(:, :, :)= reshape(dVperp1Gpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						Vperp2GLp(:, :, :)= reshape(Vperp2GLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						Vperp2GHp(:, :, :)= reshape(Vperp2GHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						Vperp2GCp(:, :, :)= reshape(Vperp2GCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						dVperp2Gp(:, :, :)= reshape(dVperp2Gpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGL2p(:, :, :)= reshape(VparGL2pp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGH2p(:, :, :)= reshape(VparGH2pp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGC2p(:, :, :)= reshape(VparGC2pp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						dVparG2p(:, :, :)= reshape(dVparG2pp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						d3vC2p(:, :, :)= reshape(d3vC2pp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))

						! ----------------------------------------------------

						do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
							write(Vperp1indstring, '(I5)') Vperp1ind
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
								write(Vperp2indstring, '(I5)') Vperp2ind
								do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
									write(Vparindstring, '(I5)') Vparind

									! ----------------------------------------------------

									Vperp1GL(1)= Vperp1GLp(Vperp1ind, Vperp2ind, Vparind)
									Vperp1GH(1)= Vperp1GHp(Vperp1ind, Vperp2ind, Vparind)
									Vperp1GC(1)= Vperp1GCp(Vperp1ind, Vperp2ind, Vparind)
									dVperp1G(1)= dVperp1Gp(Vperp1ind, Vperp2ind, Vparind)
									Vperp2GL(1)= Vperp2GLp(Vperp1ind, Vperp2ind, Vparind)
									Vperp2GH(1)= Vperp2GHp(Vperp1ind, Vperp2ind, Vparind)
									Vperp2GC(1)= Vperp2GCp(Vperp1ind, Vperp2ind, Vparind)
									dVperp2G(1)= dVperp2Gp(Vperp1ind, Vperp2ind, Vparind)
									VparGL(1)= VparGL2p(Vperp1ind, Vperp2ind, Vparind)
									VparGH(1)= VparGH2p(Vperp1ind, Vperp2ind, Vparind)
									VparGC(1)= VparGC2p(Vperp1ind, Vperp2ind, Vparind)
									dVparG(1)= dVparG2p(Vperp1ind, Vperp2ind, Vparind)
									d3vC(1)= d3vC2p(Vperp1ind, Vperp2ind, Vparind)

									if (VparGC(1) == 0d0) then
										VparGC(1)= 1d-15
									end if

									! Create nested derived data types
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT= Vperp1GL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT= Vperp1GH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT= Vperp1GC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT= dVperp1G
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT= Vperp2GL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT= Vperp2GH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT= Vperp2GC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT= dVperp2G
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT= VparGL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT= VparGH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT= VparGC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT= dVparG
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT= d3vC

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

									if ((Vperp1GL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp1GLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((Vperp1GH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp1GHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((Vperp1GC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp1GCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((dVperp1G(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										dVperp1GT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp1GT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((Vperp2GL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp2GLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((Vperp2GH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp2GHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((Vperp2GC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										Vperp2GCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((dVperp2G(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										dVperp2GT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp2GT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((VparGL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										VparGLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((VparGH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										VparGHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((VparGC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										VparGCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((dVparG(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										dVparGT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(:)) &
										/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1))) .eqv. .true.)) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
											' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
											', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
											', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
											' SUBROUTINE' // achar(27) // '[0m.'
									end if

									if ((d3vC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)% &
										d3vCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

						deallocate(Vperp1GLp, Vperp1GHp, Vperp1GCp, dVperp1Gp, Vperp2GLp, Vperp2GHp, Vperp2GCp, dVperp2Gp, &
							VparGL2p, VparGH2p, VparGC2p, dVparG2p, d3vC2p)

						! ----------------------------------------------------

					else

						! ----------------------------------------------------

						allocate(VperpGLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperpGpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVparGpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VperpGHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VperpGCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							d3vCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))

						allocate(VperpGLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVperpGp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							dVparGp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VperpGHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VperpGCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							VparGCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)), &
							d3vCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)))

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVperpCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVperpGpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVparCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVparGpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VperpGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VperpGLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VperpGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VperpGHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VperpGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VperpGCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VparGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VparGCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'd3vCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) d3vCpp(:)
						close(0)

						dVperpGp(:, :)= reshape(dVperpGpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						dVparGp(:, :)= reshape(dVparGpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VperpGLp(:, :)= reshape(VperpGLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VperpGHp(:, :)= reshape(VperpGHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VperpGCp(:, :)= reshape(VperpGCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGLp(:, :)= reshape(VparGLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGHp(:, :)= reshape(VparGHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						VparGCp(:, :)= reshape(VparGCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))
						d3vCp(:, :)= reshape(d3vCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1)/))

						! ----------------------------------------------------

						do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
							write(Vperpindstring, '(I5)') Vperpind
							do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
								write(Vparindstring, '(I5)') Vparind

								! ----------------------------------------------------

								dVperpG(1)= dVperpGp(Vperpind, Vparind)
								dVparG(1)= dVparGp(Vperpind, Vparind)
								VperpGL(1)= VperpGLp(Vperpind, Vparind)
								VperpGH(1)= VperpGHp(Vperpind, Vparind)
								VperpGC(1)= VperpGCp(Vperpind, Vparind)
								VparGL(1)= VparGLp(Vperpind, Vparind)
								VparGH(1)= VparGHp(Vperpind, Vparind)
								VparGC(1)= VparGCp(Vperpind, Vparind)
								d3vC(1)= d3vCp(Vperpind, Vparind)

								if (VparGC(1) == 0d0) then
									VparGC(1)= 1d-15
								end if

								! Create nested derived data types
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVperpGT= dVperpG
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVparGT= dVparG
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGLT= VperpGL
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGHT= VperpGH
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGCT= VperpGC
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGLT= VparGL
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGHT= VparGH
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCT= VparGC
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%d3vCT= d3vC

								! ----------------------------------------------------

								! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

								if ((dVperpG(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									dVperpGT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVperpGT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVperpGT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperpGT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((dVparG(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									dVparGT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVparGT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%dVparGT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VperpGL(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VperpGLT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGLT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGLT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGLT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VperpGH(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VperpGHT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGHT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGHT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGHT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VperpGC(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VperpGCT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGCT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VperpGCT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGCT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VparGL(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VparGLT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGLT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGLT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VparGH(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VparGHT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGHT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGHT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((VparGC(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									VparGCT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGCT(:)) &
									/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%VparGCT(1))) .eqv. .true.)) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCT HAS', &
										' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
										', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
										', AND Vparind= ', Vparind, ' IN SIMULATION PARAMETERIZATION', &
										' SUBROUTINE' // achar(27) // '[0m.'
								end if

								if ((d3vC(1) /= &
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)% &
									d3vCT(1)) .or. &
									(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

						deallocate(dVperpGp, dVparGp, VperpGLp, VperpGHp, VperpGCp, VparGLp, VparGHp, &
							VparGCp, d3vCp)

						! ----------------------------------------------------

					end if

				end do
			end do
		end do

		! ----------------------------------------------------

		! IMPORT ALL ENA EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES FROM MATLAB:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do s= 1, Stot, 1
				write(sstring, '(I5)') s
				do f= 1, SpecieT(s)%NfT(1), 1
					write(fstring, '(I5)') f
					do Qind= 1, 1, 1
						write(Qindstring, '(I5)') 1

						! ----------------------------------------------------

						allocate(VpGLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VpGHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VpGCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGLpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGHpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVpCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVqCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVphiCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVpCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVqCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVphiCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							d33vCpp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1)* &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)))

						allocate(VpGLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VpGHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VpGCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VqGCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGLp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGHp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							VphiGCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVpCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVqCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							hVphiCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVpCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVqCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							dVphiCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)), &
							d33vCp(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)))

						! ----------------------------------------------------

						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VpGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VpGLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VpGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VpGHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VpGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VpGCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VqGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VqGLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VqGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VqGHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VqGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VqGCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VphiGLmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VphiGLpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VphiGHmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VphiGHpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'VphiGCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) VphiGCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'hVpCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) hVpCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'hVqCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) hVqCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'hVphiCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) hVphiCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVpCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVpCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVqCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVqCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'dVphiCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) dVphiCpp(:)
						close(0)
						open (unit= 0, file= dataimportdir &
							// adjustl(adjustr(sstring) // '_' // &
							adjustl(adjustr(fstring) // '_' // adjustl(adjustr(Qindstring) &
							// '_' // 'd33vCmat.bin'))), &
							status= 'old', form= 'unformatted', access= 'stream')
						read(0) d33vCpp(:)
						close(0)

						VpGLp(:, :, :)= reshape(VpGLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VpGHp(:, :, :)= reshape(VpGHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VpGCp(:, :, :)= reshape(VpGCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VqGLp(:, :, :)= reshape(VqGLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VqGHp(:, :, :)= reshape(VqGHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VqGCp(:, :, :)= reshape(VqGCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VphiGLp(:, :, :)= reshape(VphiGLpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VphiGHp(:, :, :)= reshape(VphiGHpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						VphiGCp(:, :, :)= reshape(VphiGCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						hVpCp(:, :, :)= reshape(hVpCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						hVqCp(:, :, :)= reshape(hVqCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						hVphiCp(:, :, :)= reshape(hVphiCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						dVpCp(:, :, :)= reshape(dVpCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						dVqCp(:, :, :)= reshape(dVqCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						dVphiCp(:, :, :)= reshape(dVphiCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))
						d33vCp(:, :, :)= reshape(d33vCpp(:), &
							(/SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1)/))

						! ----------------------------------------------------

						do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
							write(Vpindstring, '(I5)') Vpind
							do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
								write(Vqindstring, '(I5)') Vqind
								do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1
									write(Vphiindstring, '(I5)') Vphiind

									! ----------------------------------------------------

									VpGL(1)= VpGLp(Vpind, Vqind, Vphiind)
									VpGH(1)= VpGHp(Vpind, Vqind, Vphiind)
									VpGC(1)= VpGCp(Vpind, Vqind, Vphiind)
									VqGL(1)= VqGLp(Vpind, Vqind, Vphiind)
									VqGH(1)= VqGHp(Vpind, Vqind, Vphiind)
									VqGC(1)= VqGCp(Vpind, Vqind, Vphiind)
									VphiGL(1)= VphiGLp(Vpind, Vqind, Vphiind)
									VphiGH(1)= VphiGHp(Vpind, Vqind, Vphiind)
									VphiGC(1)= VphiGCp(Vpind, Vqind, Vphiind)
									hVpC(1)= hVpCp(Vpind, Vqind, Vphiind)
									hVqC(1)= hVqCp(Vpind, Vqind, Vphiind)
									hVphiC(1)= hVphiCp(Vpind, Vqind, Vphiind)
									dVpC(1)= dVpCp(Vpind, Vqind, Vphiind)
									dVqC(1)= dVqCp(Vpind, Vqind, Vphiind)
									dVphiC(1)= dVphiCp(Vpind, Vqind, Vphiind)
									d33vC(1)= d33vCp(Vpind, Vqind, Vphiind)

									! Create nested derived data types
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGLT= VpGL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGHT= VpGH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGCT= VpGC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGLT= VqGL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGHT= VqGH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGCT= VqGC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										VphiGLT= VphiGL
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										VphiGHT= VphiGH
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										VphiGCT= VphiGC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										hVpCT= hVpC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										hVqCT= hVqC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										hVphiCT= hVphiC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										dVpCT= dVpC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										dVqCT= dVqC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										dVphiCT= dVphiC
									SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)% &
										d33vCT= d33vC

									! ----------------------------------------------------

									! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS,
									! SIZES, AND FINITE VALUES:

									if ((VpGL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VpGH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VpGC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VpGCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VqGL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VqGH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VqGC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VqGCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VphiGL(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VphiGH(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((VphiGC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%VphiGCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((hVpC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((hVqC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((hVphiC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((dVpC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((dVqC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVqCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((dVphiC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

									if ((d33vC(1) /= &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
										V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1)) .or. &
										(size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
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

						deallocate(VpGLp, VpGHp, VpGCp, VqGLp, VqGHp, VqGCp, VphiGLp, VphiGHp, VphiGCp, &
							hVpCp, hVqCp, hVphiCp, dVpCp, dVqCp, dVphiCp, d33vCp)

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

						!For (Barghouthi '98 case1) VS03:
						!lambdaPerpp(1)= 0d0 ! VLF wavelength at local gyrofrequency [m] (set to zero for long wavelength limit)
						!f0p(1)= 6d0 ! Reference gyrofrequency [Hz]
						!S0p(1)= 1d-8 ! BBELF Wave Spectral Energy Density in [(V^2/m^2)/Hz]
						!ChiPerp1p(1)= 1.7d0 ! ELF wave spectral index
						!ChiPerp2p(1)= ChiPerp1p(1) ! ELF wave spectral index

						!For (Barghouthi '97 case1 and Crew '90, Barghouthi '97 case2, Barghouthi '94, Wu '02 case3) VS1 VS2 VS3 VS10:
						!lambdaPerpp(1)= 0d0 ! VLF wavelength at local gyrofrequency [m] (set to zero for long wavelength limit)
						!f0p(1)= 5.6d0 ! Reference gyrofrequency [Hz]
						!S0p(1)= 1.2d-6 ! BBELF Wave Spectral Energy Density in [(V^2/m^2)/Hz]
						!ChiPerp1p(1)= 1.7d0 ! ELF wave spectral index
						!ChiPerp2p(1)= ChiPerp1p(1) ! ELF wave spectral index

						!For (Wu '02 case1 and Brown '95) VS8 VS11 VS12:
						!lambdaPerpp(1)= 0d0 ! VLF wavelength at local gyrofrequency [m] (set to zero for long wavelength limit)
						!f0p(1)= 6.5d0 ! Reference gyrofrequency [Hz]
						!S0p(1)= 1d-8 ! BBELF Wave Spectral Energy Density in [(V^2/m^2)/Hz]
						!ChiPerp1p(1)= 1.7d0 ! ELF wave spectral index
						!ChiPerp2p(1)= ChiPerp1p(1) ! ELF wave spectral index

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
					write(*, *) 'THREE-DIMENSIONAL KINETIC MODEL OF THE IONOSPHERIC OUTFLOW (3D-KMIO)'
					write(*, *) 'Robert M. Albarran II, Experimental and Computational Laboratory for Atmospheric and Ionospheric Research'
					write(*, *) 'Department of Physical Sciences, Embry-Riddle Aeronautical University'

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
					write(*, *) 'SIMULATION PARAMETERS:'
					write(*, *)

					write(paramstring, '(i10)') ranksize(1)
					write(*, *) trim('Number of MPI ranks= ' // adjustl(paramstring))

					if (SpecieT(s)%FluxTubeT(f)%QCell0T(1)%qGL0T(1) <= 0d0) then
						write(*, *) 'Southern Magnetic Hemisphere'
					else if (SpecieT(s)%FluxTubeT(f)%QCell0T(1)%qGL0T(1) > 0d0) then
						write(*, *) 'Northern Magnetic Hemisphere'
					end if

					write(*, *) trim('Data Import Path= ' // adjustl(dataimportdir))
					write(*, *) trim('Data Export Path= ' // adjustl(dataexportdir))

					if ((SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 1) .or. &
						(SpecieT(s)%FluxTubeT(f)%DENSITYOUTPUTflagT(1) == 1)) then
						write(*, *) trim('Density Data I/O Path= ' // adjustl(Densitydatadir))
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

					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%HiCratioMeanT
					write(*, *) trim('Mean FA Grid Length to Ion Scale Height Ratio= ' // adjustl(paramstring))
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
					write(paramstring, '(D10.2)') Vperp12Gridlinspace*1d-3
					write(*, *) trim('Vperp1 and Vperp2 Linear Grid Spacing [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12sigma*1d-3
					write(*, *) trim('Vperp1 and Vperp2 Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12sigmaFac
					write(*, *) trim('Number of Linear Vperp1 and Vperp2 Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vperp12NlinRange
					write(*, *) trim('Number of Vperp1 and Vperp2 Linear Grid Cells Spanning to &
						Last Standard Deviation= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VparGridlinspace*1d-3
					write(*, *) trim('Vpar Linear Grid Spacing [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') Vparsigma*1d-3
					write(*, *) trim('Vpar Standard Deviation [km/s]= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VparsigmaFac
					write(*, *) trim('Number of Linear Vpar Standard Deviations &
						Spanned by Linear Grid= ' // adjustl(paramstring))
					write(paramstring, '(D10.2)') VparNlinRange
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
						write(paramstring, '(D10.2)') VpphiGridlinspace*1d-3
						write(*, *) trim('Vp and Vphi Linear Grid Spacing [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') Vpphisigma*1d-3
						write(*, *) trim('Vp and Vphi Standard Deviation [km/s]= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphisigmaFac
						write(*, *) trim('Number of Linear Vp and Vpphi Standard Deviations &
							Spanned by Linear Grid= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VpphiNlinRange
						write(*, *) trim('Number of Vp and Vphi Linear Grid Cells Spanning to &
							Last Standard Deviation= ' // adjustl(paramstring))
						write(paramstring, '(D10.2)') VqGridlinspace*1d-3
						write(*, *) trim('Vq Linear Grid Spacing [km/s]= ' // adjustl(paramstring))
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
					write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%TeT(1)
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
