module PotentialStructure

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.8 PARALLEL ELECTRIC FIELD:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE FIELD-ALIGNED PARALLEL ELECTRIC FIELD [V/m] FOR ALL PARTICLES:

	subroutine PotentialStructureSub

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

				! ----------------------------------------------------

				! COMPUTE PARALLEL FIELD MAGNITUDES:

				! ----------------------------------------------------

				if (rank == 0) then

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) then
							SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind)= 0d0
						end if
						if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then
							SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind)= &
								SpecieT(s)%FluxTubeT(f)%PhiPar0BT(1)* &
								sum(SpecieT(s)%FluxTubeT(f)%dqCT(nn, 1:Qind)* &
								SpecieT(s)%FluxTubeT(f)%hqCT(nn, 1:Qind))
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
							EParR(1)= (2d0/(SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind+ 1)+ &
								SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind)))* &
								(abs(SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind+ 1)- &
								SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind))/ &
								abs(SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind+ 1)- &
								SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)))
						end if
						if (Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
							if (Qind < SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 5d0) then
								EParR(1)= (2d0/(SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind- 1)+ &
									SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind)))* &
									(abs(SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind- 1)- &
									SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind))/ &
									abs(SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)- &
									SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)))
							end if
							if (Qind >= SpecieT(s)%FluxTubeT(f)%NqUBT(1)- 5d0) then
								EParR(1)= 3d0*((2d0/(SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind- 1)+ &
									SpecieT(s)%FluxTubeT(f)%hqCT(nn, Qind)))* &
									(abs(SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind- 1)- &
									SpecieT(s)%FluxTubeT(f)%PhiParRT(nn, Qind))/ &
									abs(SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)- &
									SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))))
							end if
						end if

						! ----------------------------------------------------

						SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)= EParR(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))) &
							.eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%EPmagRT(:, :)) &
							/= (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' EPmagRT= ', SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
								' TIME-STEP= ', nn, ' IN PARALLEL FIELD SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then
							if (SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind) == 0d0) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' ZERO EPmagRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARALLEL FIELD SUBROUTINE' &
									// achar(27) // '[0m.'
							end if
						end if

						! ----------------------------------------------------

					end do

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! BROADCAST PARALLEL FIELD MAGNITUDES ACROSS ALL MPI RANKS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
				if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					call mpi_bcast(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, :), &
						((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
						MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				end if
			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE PARALLEL FORCE ACCELERATION MAGNITUDES ON STATISTICAL TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Parallel E field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) <= 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) < &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) < &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) < &
											SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) < &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) < &
											qk4(j)) .and. (qk4(j) < &
											SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! Compute Parallel E field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) > 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) > &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) > &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) > &
											SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) > &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1) > &
											qk4(j)) .and. (qk4(j) > &
											SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind- 1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EPmagInterp(1)= abs(yLinInterp(1))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
												AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn, Qind) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then

											EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind))

											if (EPmagInterp(1) == 0d0) then
												AEPmagN(j)= 0d0
											else
			                  AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
													', EPmagInterp= ', EPmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

					do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if ((isnan(real(AEPmagN(j))) .eqv. .true.) .or. &
							(size(AEPmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPmagN= ', &
								AEPmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) then

					AEPmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE PARALLEL FORCE ACCELERATION MAGNITUDES ON NON-STATISTICAL TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn /= 1d0) .and. (((nn == 2d0) .and. ((n > (nn- 2)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)+ 1d0) .and. &
				(n < (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) .or. &
				((nn > 2d0) .and. ((n >= (nn- 2)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)+ 1d0) .and. &
				(n < (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Parallel E field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLT(nn- 1, 1) <= 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn- 1, Qind) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1) < &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EPmagInterp(1)= abs(yLinInterp(1))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1) < &
										qk4(j)) .and. (qk4(j) < &
										SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EPmagInterp(1)= abs(yLinInterp(1))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%qGHT(nn- 1, Qind))) then

										EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! Compute Parallel E field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%qGLT(nn- 1, 1) > 0) then
						do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGLT(nn- 1, Qind) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1) > &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EPmagInterp(1)= abs(yLinInterp(1))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1) > &
										qk4(j)) .and. (qk4(j) > &
										SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind- 1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EPmagInterp(1)= abs(yLinInterp(1))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
											AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%qGCT(nn- 1, Qind) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%qGHT(nn- 1, Qind))) then

										EPmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EPmagRT(nn- 1, Qind))

										if (EPmagInterp(1) == 0d0) then
											AEPmagN(j)= 0d0
										else
		                  AEPmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EPmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AEPmagN(j) <= 0) .and. (EPmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AEPmagN= ', AEPmagN(j), &
												', EPmagInterp= ', EPmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

					do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
						if ((isnan(real(AEPmagN(j))) .eqv. .true.) .or. &
							(size(AEPmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AEPmagN= ', &
								AEPmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN PARALLEL FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 0) then

					AEPmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		 !if (SpecieT(s)%FluxTubeT(f)%EPARflagT(1) == 1) then
 			!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
 			!	if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
 			!		(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
 			!		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
			!			if (rank == 0) then
 			!				if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0) then

			!					!if (SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind) /= 0) then
			!						write(*, *) 'EPmagRT s, f, nn, Qind= ', s, f, nn, Qind, &
			!							SpecieT(s)%FluxTubeT(f)%EPmagRT(nn, Qind)
			!					!end if

 			!				end if
 			!			end if
 			!		end do
 			!	end if
 			!end do
 			 !do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			!	if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
 				!	(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
  			!		write(*, *) 'AEPmagN s, f= ', s, f, AEPmagN(:)
 				!end if
  			!end do
 		!end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine PotentialStructureSub

end module PotentialStructure
