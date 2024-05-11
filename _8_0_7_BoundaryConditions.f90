module BoundaryConditions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.7 BOUNDARY CONDITIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use BCReIndex
use DataTypeAllocA
use DataTypeReAllocA
use MPIReduceSum

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! RESET TOTAL PARTICLE NUMBER PER PROCESSOR ACCORDING TO LOWER AND UPPER BOUNDARY CONDITIONS:

	subroutine BoundaryConditionsSub

    ! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((n == 1) .and. (nn == 1)) then
				NsTK(1)= SpecieT(s)%FluxTubeT(f)%NsT(1)
				SpecieT(s)%FluxTubeT(f)%NsnT(nn)= NsTK(1)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAG FOR CONSISTENT PARTICLE NUMBER:

				if (SpecieT(s)%FluxTubeT(f)%NsT(1) /= SpecieT(s)%FluxTubeT(f)%NsnT(nn)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT NsT= ', SpecieT(s)%FluxTubeT(f)%NsT(1), ' AND NsnT= ', &
						SpecieT(s)%FluxTubeT(f)%NsnT(nn), ' FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
						nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

				allocate(NsTKp(1), NsTKRp(1))

				! ----------------------------------------------------

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				NsTKp(1)= NsTK(1)
				MPIRedSumIn(1)= NsTKp(1)
				call MPIReduceSumSub
				NsTKRp(1)= MPIRedSumOut(1)

				!call mpi_reduce(NsTKp(1), NsTKRp(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				if (rank == 0) then
					SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)= NsTKRp(1)
				end if

				! ----------------------------------------------------

				deallocate(NsTKp, NsTKRp)

				! ----------------------------------------------------

				call DataTypeAllocASub

				! ----------------------------------------------------

			end if
			if ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1)))) then

        ! ----------------------------------------------------

        ! DIAGNOSTIC FOR CORRECT CONFIGURATION-SPACE COUNTING:

        do j= 1, NsTK(1), 1
          if ((Qindk1(j) /= 0d0) .and. (Qindk1(j) /= -1d0) &
            .and. ((Qindk1(j) < 1d0) .and. (Qindk1(j) > &
						((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))) then
            write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
              ' INCONSISTENT Qindk1(j)= ', Qindk1(j), ' FOR PARTICLE= ', j, &
              ', SPECIE= ', s, ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
              nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
          end if
        end do

				! ----------------------------------------------------

        ! DISCARD EXISTING BOUNDARY GRID CELL IONS:

        ! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
					jcount= 0d0
					SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)= jcount
					jloopLBIonreplenish: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if (ENAflag(j) .eqv. .false.) then
							if (SpecieT(s)%FluxTubeT(f)%NqLBT(1) == Qindk1(j)) then

	              LBreplenishIonMsk(j)= .true.

								jcount= jcount+ 1d0
								SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)= jcount

								cycle jloopLBIonreplenish

							end if
	            if (SpecieT(s)%FluxTubeT(f)%NqLBT(1) /= Qindk1(j)) then

	              LBreplenishIonMsk(j)= .false.

	            end if
						end if
					end do jloopLBIonreplenish
				end if

				if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
					jcount= 0d0
					SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)= jcount
					jloopUBIonreplenish: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if (ENAflag(j) .eqv. .false.) then
							if (SpecieT(s)%FluxTubeT(f)%NqUBT(1) == Qindk1(j)) then

	              UBreplenishIonMsk(j)= .true.

								jcount= jcount+ 1d0
								SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)= jcount

								cycle jloopUBIonreplenish

							end if
	            if (SpecieT(s)%FluxTubeT(f)%NqUBT(1) /= Qindk1(j)) then

	              UBreplenishIonMsk(j)= .false.

	            end if
						end if
					end do jloopUBIonreplenish
				end if

        ! ----------------------------------------------------

        ! COMPUTE LOWER BOUNDARY ION ESCAPE FLUXES:

        ! ----------------------------------------------------

				jcount= 0d0
				SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)= jcount
				jloopLBIonoutflux: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
					if (ENAflag(j) .eqv. .false.) then
						if (0d0 == Qindk1(j)) then

              LBoutfluxIonMsk(j)= .true.

							jcount= jcount+ 1d0
							SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)= jcount

							cycle jloopLBIonoutflux

						end if
            if (0d0 /= Qindk1(j)) then

              LBoutfluxIonMsk(j)= .false.

            end if
					end if
				end do jloopLBIonoutflux

        ! ----------------------------------------------------

				!if (rank == 0) then
				!	write(*, *) 'BC6.1, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				NqLBoutfluxIon(1)= SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)
				MPIRedSumIn(1)= NqLBoutfluxIon(1)
				call MPIReduceSumSub
				NqLBoutfluxIonR(1)= MPIRedSumOut(1)

				!call mpi_reduce(NqLBoutfluxIon(1), NqLBoutfluxIonR(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				!if (rank == 0) then
				!	write(*, *) 'BC6.2, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				if (rank == 0) then

					SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(nn)= NqLBoutfluxIonR(1)* &
					 	SpecieT(s)%FluxTubeT(f)%nsnormfacT(1) ! Number of LB escaped ions (unitless)
					SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn)= & ! LB ion escape flux [m^-2 s-1]
						SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(nn)/ &
						(SpecieT(s)%FluxTubeT(f)%sigmaLBT(nn- 1)* &
						SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1)*SpecieT(s)%FluxTubeT(f)%hT(1))

					if (isnan(real(SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn))) .eqv. .true.) then
						SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn)= 0d0
					end if

					! ----------------------------------------------------

			    ! DIAGNOSTIC FLAG FOR PROPER ARRAY SIZE AND FINITE VALUE:

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn))) .eqv. .true.) &
						.or. (size(SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(:)) /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)) then
			      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBoutfluxIonRT HAS', &
			        ' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
			        f, ', AND INJECTION TIME-STEP= ', nn, &
			        ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
			    end if

					! ----------------------------------------------------

				end if

        ! ----------------------------------------------------

        ! COMPUTE UPPER BOUNDARY ION ESCAPE FLUXES:

        ! ----------------------------------------------------

				jcount= 0d0
				SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)= jcount
				jloopUBIonoutflux: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
					if (ENAflag(j) .eqv. .false.) then
						if (-1d0 == Qindk1(j)) then

              UBoutfluxIonMsk(j)= .true.

							jcount= jcount+ 1d0
							SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)= jcount

							cycle jloopUBIonoutflux

						end if
            if (-1d0 /= Qindk1(j)) then

              UBoutfluxIonMsk(j)= .false.

            end if
					end if
				end do jloopUBIonoutflux

        ! ----------------------------------------------------

				!if (rank == 0) then
				!	write(*, *) 'BC8.1, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				NqUBoutfluxIon(1)= SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)
				MPIRedSumIn(1)= NqUBoutfluxIon(1)
				call MPIReduceSumSub
				NqUBoutfluxIonR(1)= MPIRedSumOut(1)

				!call mpi_reduce(NqUBoutfluxIon(1), NqUBoutfluxIonR(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				!if (rank == 0) then
				!	write(*, *) 'BC8.2, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				if (rank == 0) then
					SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(nn)= &
						NqUBoutfluxIonR(1)*SpecieT(s)%FluxTubeT(f)%nsnormfacT(1) ! Number of UB escaped ions (unitless)
					SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn)= & ! UB ion escape flux [m^-2 s-1]
						SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(nn)/ &
						(SpecieT(s)%FluxTubeT(f)%sigmaUBT(nn- 1)* &
						SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1)*SpecieT(s)%FluxTubeT(f)%hT(1))

					if (isnan(real(SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn))) .eqv. .true.) then
						SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn)= 0d0
					end if

					! ----------------------------------------------------

			    ! DIAGNOSTIC FLAG FOR PROPER ARRAY SIZE AND FINITE VALUE:

					if ((isnan(real(SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn))) .eqv. .true.) &
						.or. (size(SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(:)) /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)) then
			      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBoutfluxIonRT HAS', &
			        ' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
			        f, ', AND INJECTION TIME-STEP= ', nn, &
			        ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
			    end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

        ! COMPUTE LOWER BOUNDARY ENA ESCAPE FLUXES:

        ! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					jcount= 0d0
					SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)= jcount
					jloopLBENAoutflux: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if (ENAflag(j) .eqv. .true.) then
							if (0d0 == Qindk1(j)) then

	              LBoutfluxENAMsk(j)= .true.

								jcount= jcount+ 1d0
								SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)= jcount

								cycle jloopLBENAoutflux

							end if
	            if (0d0 /= Qindk1(j)) then

	              LBoutfluxENAMsk(j)= .false.

	            end if
						end if
					end do jloopLBENAoutflux
				end if

        ! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

					!if (rank == 0) then
					!	write(*, *) 'BC10.1, rank, n, nn= ', rank, n, nn, NsTK(1)
					!end if

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					NqLBoutfluxENA(1)= SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)
					MPIRedSumIn(1)= NqLBoutfluxENA(1)
					call MPIReduceSumSub
					NqLBoutfluxENAR(1)= MPIRedSumOut(1)

					!call mpi_reduce(NqLBoutfluxENA(1), NqLBoutfluxENAR(1), 1, &
					!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

					!if (rank == 0) then
					!	write(*, *) 'BC10.2, rank, n, nn= ', rank, n, nn, NsTK(1)
					!end if

					if (rank == 0) then
						SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(nn)= NqLBoutfluxENAR(1)* &
						 	SpecieT(s)%FluxTubeT(f)%nsnormfacT(1) ! Number of LB escaped ENAs (unitless)
						SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn)= & ! LB ENA escape flux [m^-2 s-1]
							SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(nn)/ &
							(SpecieT(s)%FluxTubeT(f)%sigmaLBT(nn- 1)* &
							SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1)*SpecieT(s)%FluxTubeT(f)%hT(1))

						if (isnan(real(SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn))) .eqv. .true.) then
							SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn)= 0d0
						end if

						! ----------------------------------------------------

				    ! DIAGNOSTIC FLAG FOR PROPER ARRAY SIZE AND FINITE VALUE:

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn))) .eqv. .true.) &
							.or. (size(SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(:)) /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)) then
				      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' LBoutfluxENART HAS', &
				        ' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
				        f, ', AND INJECTION TIME-STEP= ', nn, &
				        ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
				    end if

						! ----------------------------------------------------

					end if
				end if

        ! ----------------------------------------------------

        ! COMPUTE UPPER BOUNDARY ENA ESCAPE FLUXES:

        ! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					jcount= 0d0
					SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)= jcount
					jloopUBENAoutflux: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if (ENAflag(j) .eqv. .true.) then
							if (-1d0 == Qindk1(j)) then

	              UBoutfluxENAMsk(j)= .true.

								jcount= jcount+ 1d0
								SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)= jcount

								cycle jloopUBENAoutflux

							end if
	            if (-1d0 /= Qindk1(j)) then

	              UBoutfluxENAMsk(j)= .false.

	            end if
						end if
					end do jloopUBENAoutflux
				end if

        ! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

					!if (rank == 0) then
					!	write(*, *) 'BC12.1, rank, n, nn= ', rank, n, nn, NsTK(1)
					!end if

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					NqUBoutfluxENA(1)= SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)
					MPIRedSumIn(1)= NqUBoutfluxENA(1)
					call MPIReduceSumSub
					NqUBoutfluxENAR(1)= MPIRedSumOut(1)

					!call mpi_reduce(NqUBoutfluxENA(1), NqUBoutfluxENAR(1), 1, &
					!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

					!if (rank == 0) then
					!	write(*, *) 'BC12.2, rank, n, nn= ', rank, n, nn, NsTK(1)
					!end if

					if (rank == 0) then
						SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(nn)= &
							NqUBoutfluxENAR(1)*SpecieT(s)%FluxTubeT(f)%nsnormfacT(1) ! Number of UB escaped ENAs (unitless)
						SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn)= & ! UB ENA escape flux [m^-2 s-1]
							SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(nn)/ &
							(SpecieT(s)%FluxTubeT(f)%sigmaUBT(nn- 1)* &
							SpecieT(s)%FluxTubeT(f)%ndatfacT(nn- 1)*SpecieT(s)%FluxTubeT(f)%hT(1))

						if (isnan(real(SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn))) .eqv. .true.) then
							SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn)= 0d0
						end if

						! ----------------------------------------------------

				    ! DIAGNOSTIC FLAG FOR PROPER ARRAY SIZE AND FINITE VALUE:

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn))) .eqv. .true.) &
							.or. (size(SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(:)) /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)) then
				      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' UBoutfluxENART HAS', &
				        ' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
				        f, ', AND INJECTION TIME-STEP= ', nn, &
				        ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
				    end if

						! ----------------------------------------------------

					end if
				end if

				! ----------------------------------------------------

				! BROADCAST LB AND UB ESCAPE FLUXES TO ALL RANKS:

				! ----------------------------------------------------

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_bcast(SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn), 1, &
					MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_bcast(SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn), 1, &
					MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					call mpi_bcast(SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn), 1, &
						MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					call mpi_bcast(SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn), 1, &
						MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				end if

        ! ----------------------------------------------------

        ! RE-INDEX TERMS TO ACCOUNT FOR PARTICLE LOSS:

        ! ----------------------------------------------------

				call BCReIndexSub

        ! ----------------------------------------------------

				! COMPUTE TOTAL NUMBER OF DISCARDED IONS:

				! ----------------------------------------------------

				if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
					.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
					dNsTK1Ion(1)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)
				end if
				if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
					.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
					dNsTK1Ion(1)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)
				end if
				if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
					.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
					dNsTK1Ion(1)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)
				end if
				if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
					.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
					dNsTK1Ion(1)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
						SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)
				end if

				! ----------------------------------------------------

				! COMPUTE TOTAL NUMBER OF ESCAPED ENAs:

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					dNsTK1ENA(1)= &
						SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)
				end if

				! ----------------------------------------------------

				! COMPUTE TOTAL NUMBER OF DISCARDED PARTICLES (IONS AND ENAs):

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then
					dNsTK1(1)= dNsTK1Ion(1)
				end if
				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					dNsTK1(1)= dNsTK1Ion(1)+ dNsTK1ENA(1)
				end if

        ! ----------------------------------------------------

        ! DIAGNOSTICS FOR CONSISTENT PARTICLE LOSS FLUX:

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)), ' FOR SPECIE= ', s, &
		            ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)), ' FOR SPECIE= ', s, &
		            ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)), ' FOR SPECIE= ', s, &
		            ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 0) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)), ' FOR SPECIE= ', s, &
		            ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION AND ENA LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 0)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION AND ENA LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 0) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION AND ENA LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					if ((SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1)) then
		        if (dNsTK1(1) /= nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)+ &
							SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
		          SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn))) then
		          write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		            ' INCONSISTENT ION AND ENA LOSS FLUX= ', dNsTK1(1), &
								nint(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(nn)+ &
								SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(nn)+ &
	              SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(nn)), &
								' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
		            nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
		        end if
					end if
				end if

        ! ----------------------------------------------------

        ! UPDATE TOTAL NUMBER OF PARTICLES PER PROCESSOR AND RE-ALLOCATE VARIABLES:

        ! ----------------------------------------------------

        dNsTK(1)= -dNsTK1(1)
        NsTK(1)= NsTK(1)+ dNsTK(1)

        ! ----------------------------------------------------

        call DataTypeReAllocASub

				! ----------------------------------------------------

        ! COMPUTE LOWER BOUNDARY INJECTION FLUX:

        ! ----------------------------------------------------

        if (SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 0) then
					dNsTK2(1)= 0d0
				end if

				if (SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 1) then

					! ----------------------------------------------------

					if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 1) then
						SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)= nint(SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(nn- 1)/ranksize(1))
					end if
					if (SpecieT(1)%FluxTubeT(1)%DENSITYPROFILEflagT(1) == 0) then
						SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)= nint(SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(nn- 1)/ranksize(1))
					end if

          dNsTK2(1)= SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)

					! ----------------------------------------------------

					! DIAGNOSTICS FOR CONSISTENT LB ION INJECTION DENSITY:

					if (dNsTK2(1) < 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE LB ION INJECTION DENSITY= ', &
							dNsTK2(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND INJECTION TIME-STEP= ', &
							nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
					end if

					if ((SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn) > SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(nn- 1)) &
						.and. (SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn) == 0d0)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NET LB ION INJECTION DENSITY= ', &
							SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn), ', NORMALIZED DENSITY= ', &
							SpecieT(s)%FluxTubeT(f)%LBNominalDensityT(nn- 1), &
							' AND LB OUTFLUX= ', SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND INJECTION TIME-STEP= ', &
							nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

        ! ----------------------------------------------------

        ! UPDATE TOTAl NUMBER OF PARTICLES PER PROCESSOR AND RE-ALLOCATE VARIABLES:

        ! ----------------------------------------------------

				dNsTK(1)= dNsTK2(1)
        NsTK(1)= NsTK(1)+ dNsTK(1)

				! ----------------------------------------------------

				call DataTypeReAllocASub

				! ----------------------------------------------------

				! COMPUTE UPPER BOUNDARY INJECTION FLUX:

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1) == 0) then
					dNsTK3(1)= 0d0
				end if

				if (SpecieT(s)%FluxTubeT(f)%UBCONDITIONflagT(1) == 1) then

					! ----------------------------------------------------

					SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn)= nint(SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(nn- 1)/ranksize(1))

					dNsTK3(1)= SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn)

					! ----------------------------------------------------

					! DIAGNOSTICS FOR CONSISTENT UB ION INJECTION DENSITY:

					if ((SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn) > SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(nn- 1)) &
						.and. (SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn) == 0d0)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NET UB ION INJECTION DENSITY= ', &
							SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn), ', NORMALIZED DENSITY= ', &
							SpecieT(s)%FluxTubeT(f)%UBNominalDensityT(nn- 1), &
							' AND UB OUTFLUX= ', SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND INJECTION TIME-STEP= ', &
							nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

				! UPDATE TOTAl NUMBER OF PARTICLES PER PROCESSOR AND RE-ALLOCATE VARIABLES:

				! ----------------------------------------------------

				dNsTK(1)= dNsTK3(1)
				NsTK(1)= NsTK(1)+ dNsTK(1)
				SpecieT(s)%FluxTubeT(f)%NsT(1)= NsTK(1)
				SpecieT(s)%FluxTubeT(f)%NsnT(nn)= NsTK(1)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAG FOR CONSISTENT PARTICLE NUMBER:

				if (SpecieT(s)%FluxTubeT(f)%NsT(1) /= SpecieT(s)%FluxTubeT(f)%NsnT(nn)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' INCONSISTENT NsT= ', SpecieT(s)%FluxTubeT(f)%NsT(1), ' AND NsnT= ', &
						SpecieT(s)%FluxTubeT(f)%NsnT(nn), ' FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ' AND INJECTION TIME-STEP= ', &
						nn, ' IN BOUNDARY CONDITIONS SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

				allocate(NsTKp(1), NsTKRp(1))

				! ----------------------------------------------------

				!if (rank == 0) then
				!	write(*, *) 'BC24.1, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				NsTKp(1)= NsTK(1)
				MPIRedSumIn(1)= NsTKp(1)
				call MPIReduceSumSub
				NsTKRp(1)= MPIRedSumOut(1)

				!call mpi_reduce(NsTKp(1), NsTKRp(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				!if (rank == 0) then
				!	write(*, *) 'BC24.2, rank, n, nn= ', rank, n, nn, NsTK(1)
				!end if

				if (rank == 0) then
					SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)= NsTKRp(1)
				end if

				! ----------------------------------------------------

				deallocate(NsTKp, NsTKRp)

				! ----------------------------------------------------

				call DataTypeReAllocASub

				! ----------------------------------------------------

			end if

		end do

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn /= 1d0) .and. (((nn == 2d0) .and. ((n > sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 2))+ 1d0) .and. &
				(n < sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))) .or. &
				((nn > 2d0) .and. ((n >= sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 2))+ 1d0) .and. &
				(n < sum(SpecieT(s)%FluxTubeT(f)%ndatfacT(1:nn- 1))))))) then

				NsTK(1)= NsTK(1)
				SpecieT(s)%FluxTubeT(f)%NsT(1)= NsTK(1)

			end if
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine BoundaryConditionsSub

end module BoundaryConditions
