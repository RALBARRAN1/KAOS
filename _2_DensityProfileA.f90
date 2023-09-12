module DensityProfileA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	2- DENSITY PROFILE A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use DensityProfileA1
use DensityProfileA2

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! SET FUNCTIONAL FORM OF DENSITY PROFILE ONTO CONFIG. SPACE CELL-CENTERED .m files:

	subroutine DensityProfileASub

		! ----------------------------------------------------

		! COMPUTE PARTICLE NUMBER DENSITY OVER ENTIRE CONFIGURATION-SPACE GRID:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then

					allocate(ICbbp(((NqUB(1)- NqLB(1))+ 3)))

					do Qind= NqLB(1), NqUB(1)+ 2, 1

						! ----------------------------------------------------

						! COMPUTE STEADY-STATE ION DENSITY PROFILE WITH NORMALIZATION:
						! NOTE: SET ninormfac SUCH THAT NO FA CELL IS EMPTY.

						ICbbp(Qind)= 2d0*GG*ME*((SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)* &
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)* &
							cos(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)))/ &
							((SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)**2d0)* &
							(sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))))

					end do

					IC0bb(1)= sum(ICbbp(1:SpecieT(s)%Qindns0T(1)))

					do Qind= NqLB(1), NqUB(1)+ 2, 1

						!Altitude dependent gravitational acceleration, integration over ds
						ICbb(1)= sum(ICbbp(1:Qind))
	          gC(1)= ICbb(1)- IC0bb(1)

						if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 1) then
		          argC(1)= (SpecieT(s)%msT(1)*gC(1))/ &
		           	(kB*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)+ SpecieT(s)%FluxTubeT(f)%TeT(1)))
						end if

						if (SpecieT(s)%FluxTubeT(f)%EAMBflagT(1) == 0) then
							argC(1)= (SpecieT(s)%msT(1)*gC(1))/ &
		           	(kB*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)+ SpecieT(s)%FluxTubeT(f)%TeT(1)))
						end if

						nsC(1)= SpecieT(s)%FluxTubeT(f)%ns0T(1)* &
							exp(argC(1)) ! Number of ions per flux-tube grid cell [m^-3]

						! Normalized number of macroparticles per flux-tube grid cell [unitless]
						nsnormC(1)= (nsC(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
							(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))

						! Create nested derived data types

						if (SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 0) then
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT= nsnormC
						end if
						if (SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 1) then
							if (Qind == NqLB(1)) then
								SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= SpecieT(s)%FluxTubeT(f)%nsnormCLBInputT(1)
							end if
							if (Qind == NqUB(1)+ 2) then
								SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= SpecieT(s)%FluxTubeT(f)%nsnormCUBInputT(1)
							end if
							if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1)+ 2)) then
								SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)= &
									(SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%DensityInputT(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
									(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))
							end if
						end if

						! ----------------------------------------------------

						! COMPUTE STEADY-STATE O NEUTRAL DENSITY PROFILE:

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

							gCNeut(1)= (GG*ME)/(kB*TNeut)

							argCNeut(1)= mNeut*gCNeut(1)*((1d0/SpecieT(s)%FluxTubeT(f)%zns0NeutT(1))- &
								(1d0/SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)))

							nsCNeut(1)= SpecieT(s)%FluxTubeT(f)%ns0NeutT(1)* &
								exp(-argCNeut(1)) ! Number of H neutrals per flux-tube grid cell [m^-3]

							! Normalized number of macroparticles per flux-tube grid cell [unitless]
							nsnormCNeut(1)= (nsCNeut(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
								(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))

							! Create nested derived data types
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T= nsnormCNeut

						end if

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

						if ((size(gC(:)) /= 1) .or. (isnan(real(gC(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gC HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(argC(:)) /= 1) .or. (isnan(real(argC(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argC HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(nsC(:)) /= 1) .or. (isnan(real(nsC(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsC HAS', &
								' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%DENSITYINPUTflagT(1) == 0) then
							if ((nsnormC(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(:)) /= 1) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCT(1))) .eqv. &
								.true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
									' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

							if ((size(gCNeut(:)) /= 1) .or. (isnan(real(gCNeut(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gCNeut HAS', &
									' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
									' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(argCNeut(:)) /= 1) .or. (isnan(real(argCNeut(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argCNeut HAS', &
									' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
									' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(nsnormCNeut(:)) /= 1) .or. (isnan(real(nsnormCNeut(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeut HAS', &
									' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
									' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((nsnormCNeut(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(1)) .or. &
								(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(:)) /= 1) .or. &
								(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(1))) .eqv. &
								.true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeutT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
									' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
							end if

						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				end if
			end do
		end do

		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 1) then
			call DensityProfileA1Sub
		end if
		if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 0) then
			call DensityProfileA2Sub
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR CONSISTENT TOTAL PARTICLE NUMBER:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (SpecieT(s)%FluxTubeT(f)%NsT(1) /= sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:))) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT TOTAL PARTICLE NUMBER= ', SpecieT(s)%FluxTubeT(f)%NsT(1), &
						' AND FLUX-TUBE SUMMATION= ', sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)), &
						' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
						' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
				end if

				if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 0) then
					if (rank == 0) then
						if (sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' INCONSISTENT CONFIG-SPACE PARTICLE NUMBER SUMMATION= ', sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)), &
								' AND TOTAL PARTICLE NUMBER= ', SpecieT(s)%FluxTubeT(f)%NsRRT(1), &
								' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if
					end if
				end if
			end do
		end do

		! ----------------------------------------------------

		! SET BOUNDARY GHOST CELL DENSITIES:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

					if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) then
						SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)
						SpecieT(s)%FluxTubeT(f)%sigmaLBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)/ &
							(SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%hqC0T(1)* &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%dqC0T(1))
					end if

					if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%d3xCLBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)
						SpecieT(s)%FluxTubeT(f)%sigmaLBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%d3xC0T(1)/ &
							(SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%hqC0T(1)* &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqLB(1)+ 1)%dqC0T(1))
					end if

					if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0) then
						SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)
						SpecieT(s)%FluxTubeT(f)%sigmaUBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)/ &
							(SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%hqC0T(1)* &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%dqC0T(1))
					end if

					if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%d3xCUBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)
						SpecieT(s)%FluxTubeT(f)%sigmaUBT(1)= &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%d3xC0T(1)/ &
							(SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%hqC0T(1)* &
							SpecieT(s)%FluxTubeT(f)%QCell0T(NqUB(1)+ 1)%dqC0T(1))
					end if

					SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1)= nsnormCLB(1)
					SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(1)= nsnormCUB(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS CONSISTENT LOWER AND UPPER BOUNDARY GHOST CELL DENSITY:

					if (SpecieT(1)%FluxTubeT(1)%LBCONDITIONflagT(1) == 1) then
						if (nsnormCLB(1) == 0d0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' ZERO LOWER BOUNDARY DENSITY= ', nsnormCLB(1), &
								' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1) == 0d0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' NOMINAL LOWER BOUNDARY DENSITY= ', SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(1), &
								' FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
								' IN DENSITY PROFILE A SUBROUTINE' // achar(27) // '[0m.'
						end if
					end if

					! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! RE-INDEX CONFIGURATION SPACE GRID FOR A NON-COMPUTATIONAL LOWER BOUNDARY GHOST CELL:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				do Qind= NqLB(1), NqUB(1), 1

					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%hqCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%hqC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dpCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dpC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dqCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dqC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dphiCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%dphiC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%rGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%rGC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%phiGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%phiGC0T(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR CONSISTENT PHI VALUE:

					if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1) /= SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT INITIAL phiGC0T VALUE= ', &
							SpecieT(s)%FluxTubeT(f)%QCell0T(1)%phiGC0T(1), ' AND phiGC0T VALUE= ', &
							SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Q CELL= ', Qind, &
							' IN DENSITY PROFILE A', ' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%thetaGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%thetaGC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ellGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%ellGC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGL0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%qGH0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%pGCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%pGC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%d3xCT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TsPerpT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPerp0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TsParT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%TsPar0T(1)
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TsT(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%Ts0T(1)

					SpecieT(s)%FluxTubeT(f)%rGCTp(Qind)= &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%rGCT(1)
					SpecieT(s)%FluxTubeT(f)%hqCTp(Qind)= &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%hqCT(1)
					SpecieT(s)%FluxTubeT(f)%dqCTp(Qind)= &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dqCT(1)

					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%nsnormCNeutT(1)= nsnormCNeut0(Qind+ 1)
					end if

				end do
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! PRINT TERMINAL EXPORT PARAMETERS:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					write(*, *) 'Cluster Initial Ion Macro-particle Number per Config-Space Grid Cell= ', &
						SpecieT(s)%FluxTubeT(f)%NsFARpT(:)
					write(*, *)
					write(paramstring, '(i10)') sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:))
					write(*, *) trim('Cluster Total Initial Ion Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') nsnormCLB(1)
					write(*, *) trim('Cluster LB Nominal Injection Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') nsnormCUB(1)
					write(*, *) trim('Cluster UB Nominal Injection Macro-particle Number= ' // adjustl(paramstring))
					write(paramstring, '(i10)') SpecieT(s)%FluxTubeT(f)%NsT(1)
					write(*, *) trim('Root Rank Total Initial Ion Macro-particle Number= ' // adjustl(paramstring))

					if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 1) then
						write(paramstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)
						write(*, *) trim('Macro-particle Normalization Constant= ' // adjustl(paramstring))
					end if
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					call cpu_time(S2End)
					write(S2string, '(i10)')  nint(S2End)
					write(*, *) trim('%% 2- RANK= ' // adjustl(rankstring)) // &
						trim(', REAL-TIME= ' // adjustl(S2string)) // &
						trim(' s. DENSITY PROFILE A COMPLETE %%')
				end if
			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DensityProfileASub

end module DensityProfileA
