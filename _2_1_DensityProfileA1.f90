module DensityProfileA1

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	2.1- DENSITY PROFILE A1:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MPIReduceSum

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine DensityProfileA1Sub

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				allocate(SpecieT(s)%FluxTubeT(f)%NsFARRpT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))
				allocate(SpecieT(s)%FluxTubeT(f)%NsFApT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))
				if (rank == 0) then

					! ----------------------------------------------------

					allocate(NsFARRp(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						! Number of particles in each FA cell
						NsFARRp(Qind)= nint(SpecieT(s)%FluxTubeT(f)%QCell0T( &
							SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ Qind- 1+ 1)%nsnormCT(1))

						NsFARR(1)= NsFARRp(Qind)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER GRIDDED PARTICLE NUMBERS:

						if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
							if (NsFARRp(Qind) /= nint((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%DensityInputT(1)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
								(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1)))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRp= ', NsFARRp(Qind), &
								  ' DOES NOT EQUAL DensityInputT= ', &
									nint((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%DensityInputT(1)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
									(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%d3xC0T(1))), &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
									' AND Qind= ', Qind, ' IN DENSITY PROFILE A1 SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
						end if

						if (NsFARRp(Qind) /= NsFARR(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRp DOES', &
								' NOT EQUAL NsFARRp FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								' AND Qind= ', Qind, ' IN DENSITY PROFILE A1 SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						if (NsFARRp(Qind) < 0d0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRp', &
								' IS LESS THAN OR EQUALS ZERO FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								' AND Qind= ', Qind, ' IN DENSITY PROFILE A1 SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						! Create nested derived data types
						SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT= NsFARR

						! ----------------------------------------------------

					end do

					NsRR(1)= sum(NsFARRp(:)) ! Total number of particles

					! Create nested derived data types
					SpecieT(s)%FluxTubeT(f)%NsFARRpT= NsFARRp ! Number of particles in cluster per q cell
					SpecieT(s)%FluxTubeT(f)%NsRRT= NsRR ! Number of particles in cluster in all q cells

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						if ((NsFARRp(Qind) /= SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(:)) /= 1) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1))) .eqv. &
							.true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ' AND Qind= ', Qind, &
								' IN DENSITY PROFILE A1 SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((NsFARRp(Qind) /= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)) .or. &
							(size(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:)) /= &
							SpecieT(s)%FluxTubeT(f)%NqICT(1)) .or. &
							(isnan(real(SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind))) .eqv. &
							.true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRpT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ' AND Qind= ', Qind, &
								' IN DENSITY PROFILE A1 SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					if ((NsRR(1) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) .or. &
						(size(SpecieT(s)%FluxTubeT(f)%NsRRT(:)) /= 1) .or. &
						(isnan(real(SpecieT(s)%FluxTubeT(f)%NsRRT(1))) .eqv. .true.)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsRRT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', &
						' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (SpecieT(s)%FluxTubeT(f)%NsRRT(1) <= 7) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsRRT VALUE', &
						' IS NOT COMPATIBLE WITH TEST EXPORT PARTICLE NUMBER', &
						' FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', &
						' SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:)) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRpT DOES', &
							' NOT SUM TO NsRRT FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
							' IN DENSITY PROFILE A1 SUBROUTINE' // achar(27) // '[0m.'
					end if

					if (sum(NsFARRp(:)) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NsFARRp DOES', &
							' NOT SUM TO NsRRT FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, &
							' IN DENSITY PROFILE A1 SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					deallocate(NsFARRp)

					! ----------------------------------------------------

				end if
			end do
		end do

		if (rank == 0) then

			! ----------------------------------------------------

			! COMPUTE TOTAL NUMBER OF PARTICLES OVER ALL CONFIGURATION-SPACE:

			! ----------------------------------------------------

			allocate(NsRRTotpp(Stot))

			! ----------------------------------------------------

			do s= 1, Stot, 1

				! ----------------------------------------------------

				allocate(SpecieT(s)%NsRRTotpT(SpecieT(s)%NfT(1)))

				! ----------------------------------------------------

				do f= 1, SpecieT(s)%NfT(1), 1
					SpecieT(s)%NsRRTotpT(f)= SpecieT(s)%FluxTubeT(f)%NsRRT(1)
				end do
				SpecieT(s)%NsRRTotpST(1)= sum(SpecieT(s)%NsRRTotpT(:))
				NsRRTotpp(s)= SpecieT(s)%NsRRTotpST(1)

			end do

			! Total number of simulated particles over all particle species and flux tubes
			NsRRTot(1)= sum(NsRRTotpp(:))

			! ----------------------------------------------------

			deallocate(NsRRTotpp)
			do s= 1, Stot, 1
				deallocate(SpecieT(s)%NsRRTotpT)
			end do

			! ----------------------------------------------------

			! Compute number of initial particles per mpi rank
			do s= 1, Stot, 1
				do f= 1, SpecieT(s)%NfT(1), 1

					! Number of particles per rank in all q cells
					SpecieT(s)%FluxTubeT(f)%NsT(1)= SpecieT(s)%FluxTubeT(f)%NsRRT(1)/ranksize(1)
					! Number of particles in cluster in all q cells
					SpecieT(s)%FluxTubeT(f)%NsRRT(1)= SpecieT(s)%FluxTubeT(f)%NsT(1)*ranksize(1)

					NsFARRpTSum(1)= sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:))

					if (NsFARRpTSum(1) > SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						Qloop: do Qind= SpecieT(s)%FluxTubeT(f)%NqUBT(1), SpecieT(s)%FluxTubeT(f)%NqLBT(1), -1
							if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) <= &
								(NsFARRpTSum(1)- SpecieT(s)%FluxTubeT(f)%NsRRT(1))) then
								SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) !0d0
								NsFARRpTSum(1)= sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:))
							end if
							if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) > &
								(NsFARRpTSum(1)- SpecieT(s)%FluxTubeT(f)%NsRRT(1))) then
								SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)= &
									SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)- &
									(NsFARRpTSum(1)- SpecieT(s)%FluxTubeT(f)%NsRRT(1))
								NsFARRpTSum(1)= sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:))
							end if
							if (NsFARRpTSum(1) == SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
								exit Qloop
							end if
						end do Qloop
					end if

					! ----------------------------------------------------

					allocate(NsFARRmsk(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) /= 0d0) then
							NsFARRmsk(Qind)= .true.
						else
							NsFARRmsk(Qind)= .false.
						end if
					end do

					allocate(NsFARRpTMP(SpecieT(s)%FluxTubeT(f)%NqICT(1)))
					NsFARRpTMP(:)= SpecieT(s)%FluxTubeT(f)%NsFARRpT(:)
					deallocate(SpecieT(s)%FluxTubeT(f)%NsFARRpT)

					SpecieT(s)%FluxTubeT(f)%NqICT(1)= SpecieT(s)%FluxTubeT(f)%NqICT(1)- count(NsFARRmsk .eqv. .false.)
					SpecieT(s)%FluxTubeT(f)%NqUBT(1)= SpecieT(s)%FluxTubeT(f)%NqUBT(1)- count(NsFARRmsk .eqv. .false.)

					! ----------------------------------------------------

					! DIAGNOSTIC FOR CORRECT INITIALIZATION GRID CELLS:

					if (SpecieT(s)%FluxTubeT(f)%NqICT(1) /= &
						(abs(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)))+ 1d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NqICT=  ', &
							SpecieT(s)%FluxTubeT(f)%NqICT(1), ', AND (NqICB- NqICA)+ 1= ', &
							(abs(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)))+ 1d0, &
							' VALUES FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

					allocate(SpecieT(s)%FluxTubeT(f)%NsFARRpT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)= NsFARRpTMP(Qind)

						! ----------------------------------------------------

						! DIAGNOSTIC FOR CORRECT CLUSTER PARTICLE NUMBER:

						if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) <= 0d0) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE OR ZERO NsFARRpT= ', &
								SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
								', AND Q CELL= ', Qind, ' IN DENSITY PROFILE A1 SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					! DIAGNOSTIC FOR CORRECT CLUSTER PARTICLE NUMBER:

					if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:)) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NsFARRpT SUMMATION=  ', &
							sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:)), ', AND NsRRT= ', SpecieT(s)%FluxTubeT(f)%NsRRT(1), &
							' VALUES FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end do
			end do

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		! BROADCAST NUMBER OF INITIAL PARTICLES PER RANK:

		! ----------------------------------------------------

		do s= 1, Stot, 1
		  do f= 1, SpecieT(s)%NfT(1), 1

		    ! ----------------------------------------------------

		    ! Broadcast total particle number to all ranks
		    call mpi_barrier(MPI_COMM_WORLD, ierr)
		    call mpi_bcast(SpecieT(s)%FluxTubeT(f)%NsT(1), int(1), &
		      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

		    call mpi_barrier(MPI_COMM_WORLD, ierr)
		    call mpi_bcast(SpecieT(s)%FluxTubeT(f)%NsFARRpT(:), int(SpecieT(s)%FluxTubeT(f)%NqICT(1)), &
		      MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

		    ! ----------------------------------------------------

				SpecieT(s)%FluxTubeT(f)%nsnormCLBT(1)= SpecieT(s)%FluxTubeT(f)%LBNominalDensityGT(1)
				SpecieT(s)%FluxTubeT(f)%nsnormCUBT(1)= SpecieT(s)%FluxTubeT(f)%UBNominalDensityGT(1)

		    ! ----------------------------------------------------

		    ! Broadcast neutral densities to all ranks
		    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

		      allocate(nsnormCNeut0(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 3)))

		      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1
		        nsnormCNeut0(Qind)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%nsnormCNeut0T(1)
		      end do

		      call mpi_barrier(MPI_COMM_WORLD, ierr)
		      call mpi_bcast(nsnormCNeut0(:), int((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 3), &
		        MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

		    end if

		  end do
		end do

		! ----------------------------------------------------

		! SUM TOTAL NUMBER OF PARTICLES OVER ALL RANKS:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				Ns(1)= SpecieT(s)%FluxTubeT(f)%NsT(1)
				MPIRedSumIn(1)= Ns(1)
				call MPIReduceSumSub
				NsR(1)= MPIRedSumOut(1)

				!call mpi_reduce(Ns(1), NsR(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				if (rank == 0) then
					SpecieT(s)%FluxTubeT(f)%NsRT(1)= NsR(1)
				end if

			end do
		end do

		! ----------------------------------------------------

		! DIAGNOSTIC FOR CORRECT PARTICLE NUMBER PER MPI RANK:

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1
				if (rank == 0) then
					if (SpecieT(s)%FluxTubeT(f)%NsRT(1) /= SpecieT(s)%FluxTubeT(f)%NsRRT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NsRT=  ', &
							SpecieT(s)%FluxTubeT(f)%NsRT(1), ', AND NsRRT= ', SpecieT(s)%FluxTubeT(f)%NsRRT(1), &
							' VALUES FOR SPECIE= ', s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', &
							' SUBROUTINE' // achar(27) // '[0m.'
					end if
				end if
			end do
		end do

		! ----------------------------------------------------

		! COMPUTE NUMBER OF PARTICLES IN EACH GRID CELL FOR EACH RANK:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				do rr= 0, ranksize(1)- 1, 1
					if ((rank == 0) .and. (rr == rank)) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(1:Qind)) < SpecieT(s)%FluxTubeT(f)%NsT(1)) then
								SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)
							end if
							if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(1:Qind)) >= SpecieT(s)%FluxTubeT(f)%NsT(1)) then
								if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(1:Qind- 1)) <= SpecieT(s)%FluxTubeT(f)%NsT(1)) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= &
										(SpecieT(s)%FluxTubeT(f)%NsT(1))- sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(1:Qind- 1))

									RankT(rr+ 1)%NsTqindT(1)= Qind ! Grid cell index for last particle
									RankT(rr+ 1)%NsTFAindT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind) ! FA index for last particle
									NsTqind(1)= RankT(rr+ 1)%NsTqindT(1)
									NsTFAind(1)= RankT(rr+ 1)%NsTFAindT(1)

								end if
								if (sum(SpecieT(s)%FluxTubeT(f)%NsFARRpT(1:Qind- 1)) > SpecieT(s)%FluxTubeT(f)%NsT(1)) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= 0d0
								end if
							end if
						end do

						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)

							! ----------------------------------------------------

							! DIAGNOSTIC FOR POSITIVE PARTICLE NUMBER PER CONFIG-SPACE GRID CELL:

							if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1) < 0d0) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE NsFAT=  ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), ' FOR SPECIE= ', &
									s, ' AND FLUX TUBE= ', f, ', Q CELL= ', Qind, &
									' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

						end do

						! ----------------------------------------------------

						! DIAGNOSTIC FOR CORRECT PARTICLE NUMBER PER CONFIG-SPACE GRID CELL:

						if (sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NsT=  ', &
								SpecieT(s)%FluxTubeT(f)%NsT(1), ', AND NsFAT SUMMATION= ', &
								sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:)), ' VALUES FOR SPECIE= ', &
								s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						! Send last particle position to next rank
						if (ranksize(1) /= 1) then
							call mpi_send(NsTqind, 1, MPI_INT, rank+ 1, 0, MPI_COMM_WORLD, ierr)
							call mpi_send(NsTFAind, 1, MPI_INT, rank+ 1, 0, MPI_COMM_WORLD, ierr)
						end if

						! ----------------------------------------------------

					end if
				end do

				! ----------------------------------------------------

				do rr= 0, ranksize(1)- 1, 1
					if ((rank /= 0) .and. (rr == rank)) then

						! ----------------------------------------------------

						! Receive last particle position from previous rank
						call mpi_recv(NsTqind, 1, MPI_INT, rank- 1, 0, MPI_COMM_WORLD, status, ierr)
						call mpi_recv(NsTFAind, 1, MPI_INT, rank- 1, 0, MPI_COMM_WORLD, status, ierr)
						RankT(rr)%NsTqindT(1)= NsTqind(1)
						RankT(rr)%NsTFAindT(1)= NsTFAind(1)

						! ----------------------------------------------------

						Qindloop: do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							if (Qind < NsTqind(1)) then
								SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= 0d0
							end if
							if (Qind == NsTqind(1)) then
								if ((SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)- NsTFAind(1)) &
									< (SpecieT(s)%FluxTubeT(f)%NsT(1)- &
									sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= &
										(SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)- NsTFAind(1))
								end if
								if ((SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)- NsTFAind(1)) &
									>= (SpecieT(s)%FluxTubeT(f)%NsT(1)- &
									sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= &
										(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
										sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))

									RankT(rr+ 1)%NsTqindT(1)= Qind ! Grid cell index for last particle
									RankT(rr+ 1)%NsTFAindT(1)= NsTFAind(1)+ SpecieT(s)%FluxTubeT(f)%NsFApT(Qind) ! FA index for last particle
									NsTqind(1)= RankT(rr+ 1)%NsTqindT(1)
									NsTFAind(1)= RankT(rr+ 1)%NsTFAindT(1)

									exit Qindloop

								end if
							end if

							if (Qind > NsTqind(1)) then
								if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) < &
									(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
									sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)
								end if
								if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) >= &
									(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
									sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
									if ((SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind- 1)- NsTFAind(1)) <= &
										(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
										sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
										SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= &
											(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
											sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))

										RankT(rr+ 1)%NsTqindT(1)= Qind ! Grid cell index for last particle
									  RankT(rr+ 1)%NsTFAindT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind) ! FA index for last particle
										NsTqind(1)= RankT(rr+ 1)%NsTqindT(1)
										NsTFAind(1)= RankT(rr+ 1)%NsTFAindT(1)

										exit Qindloop

									end if
									if ((SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind- 1)- NsTFAind(1)) > &
										(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
										sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
										SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= SpecieT(s)%FluxTubeT(f)%NsT(1)- &
										sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1))

										RankT(rr+ 1)%NsTqindT(1)= Qind ! Grid cell index for last particle
									  RankT(rr+ 1)%NsTFAindT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind) ! FA index for last particle
										NsTqind(1)= RankT(rr+ 1)%NsTqindT(1)
										NsTFAind(1)= RankT(rr+ 1)%NsTFAindT(1)

										exit Qindloop

									end if
								end if
							end if
							if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
								if (SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind) == &
									(SpecieT(s)%FluxTubeT(f)%NsT(1)- &
									sum(SpecieT(s)%FluxTubeT(f)%NsFApT(1:Qind- 1)))) then
									SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)

									RankT(rr+ 1)%NsTqindT(1)= Qind ! Grid cell index for last particle
									RankT(rr+ 1)%NsTFAindT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind) ! FA index for last particle
									NsTqind(1)= RankT(rr+ 1)%NsTqindT(1)
									NsTFAind(1)= RankT(rr+ 1)%NsTFAindT(1)

									exit Qindloop

								end if
							end if
						end do Qindloop
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							if (Qind > NsTqind(1)) then
								SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= 0d0
							end if

							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)

							! ----------------------------------------------------

							! DIAGNOSTIC FOR POSITIVE PARTICLE NUMBER PER CONFIG-SPACE GRID CELL:

							if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1) < 0d0) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE NsFAT=  ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1), ' FOR SPECIE= ', &
									s, ' AND FLUX TUBE= ', f, ', Q CELL= ', Qind, &
									' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

						end do

						! ----------------------------------------------------

						! DIAGNOSTIC FOR CORRECT PARTICLE NUMBER PER CONFIG-SPACE GRID CELL:

						if (sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NsT=  ', &
								SpecieT(s)%FluxTubeT(f)%NsT(1), ', AND NsFAT SUMMATION= ', &
								sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:)), ' VALUES FOR SPECIE= ', &
								s, ' AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

						! Send last particle position to next rank
						if (rr < ranksize(1)- 1) then
							call mpi_send(NsTqind, 1, MPI_INT, rank+ 1, 0, MPI_COMM_WORLD, ierr)
							call mpi_send(NsTFAind, 1, MPI_INT, rank+ 1, 0, MPI_COMM_WORLD, ierr)
						end if

						! ----------------------------------------------------

					end if
				end do

				! ----------------------------------------------------

			end do
		end do

		! ----------------------------------------------------

		! SEND ALL REMAININNG PARTICLES TO ROOT RANK:

    ! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				if (rank == 0) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
						if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1) /= SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind)) then

							SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)= &
								SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)+ &
								(SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1)- SpecieT(s)%FluxTubeT(f)%NsFARRpT(Qind))

							SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFAT(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)

						end if
					end do

					SpecieT(s)%FluxTubeT(f)%NsT(1)= sum(SpecieT(s)%FluxTubeT(f)%NsFApT(:))

				end if

			end do
		end do

		! ----------------------------------------------------

		! SUM PARTICLES PER GRID CELL OVER ALL RANKS:

		! ----------------------------------------------------

		do s= 1, Stot, 1
			do f= 1, SpecieT(s)%NfT(1), 1

				! ----------------------------------------------------

				if (rank == 0) then
					allocate(SpecieT(s)%FluxTubeT(f)%NsFARpT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))
				end if

				! ----------------------------------------------------

				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					NsFAp(1)= SpecieT(s)%FluxTubeT(f)%NsFApT(Qind)
					MPIRedSumIn(1)= NsFAp(1)
					call MPIReduceSumSub
					NsFARp(1)= MPIRedSumOut(1)

					!call mpi_reduce(NsFAp(1), NsFARp(1), 1, &
					!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

					if (rank == 0) then
						SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind)= NsFARp(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FOR CORRECT TOTAL CLUSTER PARTICLE NUMBER:

						if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1) &
							/= SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ROOT RANK NFARRT=  ', &
								SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT, &
								', AND NsFARpT= ', SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind), ' VALUES FOR SPECIE= ', &
								s, ', FLUX TUBE= ', f, ' AND Q CELL= ', Qind, &
								' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
							if (SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1) &
								/= nint(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%nsnormCT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ROOT RANK NsFARRT=  ', &
									SpecieT(s)%FluxTubeT(f)%QCellICT(Qind)%NsFARRT(1), &
									', AND nsnormCT= ', nint(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%nsnormCT(1)), ' VALUES FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ' AND Q CELL= ', Qind, &
									' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
							end if
							if (SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind) &
								/= nint(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%nsnormCT(1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ROOT RANK NsFARpT=  ', &
									SpecieT(s)%FluxTubeT(f)%NsFARpT(Qind), &
									', AND nsnormCT= ', nint(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind+ 1)%nsnormCT(1)), ' VALUES FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ' AND Q CELL= ', Qind, &
									' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

					end if
				end do

				if (rank == 0) then

					! ----------------------------------------------------

					! DIAGNOSTIC FOR CORRECT CLUSTER PARTICLE NUMBER PER CONFIG-SPACE GRID CELL:

					if (NsRR(1) /= sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT ROOT RANK NsRR=  ', &
							NsRR(1), ', AND NsFARpT SUM= ', sum(SpecieT(s)%FluxTubeT(f)%NsFARpT(:)), ' VALUES FOR SPECIE= ', &
							s, ', AND FLUX TUBE= ', f, ' IN DENSITY PROFILE A1', ' SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

			end do
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DensityProfileA1Sub

end module DensityProfileA1
