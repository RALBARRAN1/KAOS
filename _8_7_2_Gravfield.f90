module Gravfield

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.7.2 GRAVITATIONAL FIELD:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE FIELD-ALIGNED GRAVITATIONAL FIELD [V/m] FOR ALL PARTICLES:

	subroutine GravfieldSub

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((n == 1) .and. (nn == 1)) then

				! ----------------------------------------------------

				! COMPUTE GRAVITATIONAL FIELD MAGNITUDES:

				! ----------------------------------------------------

				if (rank == 0) then

					! ----------------------------------------------------

					do Qind= NqLB(1), NqUB(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) then
							EGravR(1)= 0d0
						end if

						! ----------------------------------------------------

						! COMPUTE GRAVITATIONAL FIELD:

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then
							EGravR(1)= (-2d0*GG*ME*SpecieT(s)%msT(1)*cos(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%thetaGCT(1))/ &
	              (SpecieT(s)%qsT(1)*(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%rGCT(1)**2d0)* &
	              sqrt(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%ellGCT(1))))
						end if

						! ----------------------------------------------------

						SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)= EGravR(1)

						! ----------------------------------------------------

						! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

						if ((isnan(real(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))) &
							.eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%EGmagRT(:)) &
							/= ((NqUB(1)- NqLB(1))+ 1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
								' EGmagRT= ', SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind), &
								' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, &
								' IN GRAVITATIONAL FIELD SUBROUTINE' &
								// achar(27) // '[0m.'
						end if

						if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then
							if (SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind) == 0d0) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' ZERO EGmagRT VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, &
									' IN GRAVITATIONAL FIELD SUBROUTINE' &
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

		! BROADCAST GRAVITATIONAL FIELD MAGNITUDES ACROSS ALL MPI RANKS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((n == 1) .and. (nn == 1)) then
				if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then

					call mpi_barrier(MPI_COMM_WORLD, ierr)
					call mpi_bcast(SpecieT(s)%FluxTubeT(f)%EGmagRT(:), &
						((NqUB(1)- NqLB(1))+ 1), &
						MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

				end if
			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE GRAVITATIONAL FORCE ACCELERATION MAGNITUDES ON STATISTICAL TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Gravitational field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then
						do Qind= NqLB(1), NqUB(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == NqLB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == NqUB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) < &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) <= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == NqLB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == NqUB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
											qk4(j)) .and. (qk4(j) < &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) <= &
		                  qk4(j)) .and. (qk4(j) <= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
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

					! Compute Gravitational field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then
						do Qind= NqLB(1), NqUB(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (nn == 1) then
									if (Qind == NqLB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn == 1) then
									if (Qind == NqUB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
											SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
											(SpecieT(s)%FluxTubeT(f)%q0T(j) > &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= SpecieT(s)%FluxTubeT(f)%q0T(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) >= &
		                  SpecieT(s)%FluxTubeT(f)%q0T(j)) .and. &
		                  (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								! ----------------------------------------------------

								if (nn /= 1) then
									if (Qind == NqLB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if
									end if
								end if

								if (nn /= 1) then
									if (Qind == NqUB(1)) then
										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
											qk4(j)) .and. (qk4(j) > &
											SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

											xLinInterp(1)= qk4(j)
											xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
											xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
											yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
											yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

											call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
												yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

											EGmagInterp(1)= abs(yLinInterp(1))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
												AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
											end if

											! ----------------------------------------------------

										end if

										if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) >= &
		                  qk4(j)) .and. (qk4(j) >= &
		                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

											EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

											if (EGmagInterp(1) == 0d0) then
												AGmagN(j)= 0d0
											else
			                  AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
											end if

											! ----------------------------------------------------

											! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

											if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
												write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
													' STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
													', EGmagInterp= ', EGmagInterp(1), &
													' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
													f, ', Q GRID CELL= ', Qind, ', STATISTICAL TIME-STEP= ', nn, &
													', AND PARTICLE= ', j, &
													' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
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
						if ((isnan(real(AGmagN(j))) .eqv. .true.) .or. &
							(size(AGmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGmagN= ', &
								AGmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) then

					AGmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ALL PARTICLE GRAVITATIONAL FORCE ACCELERATION MAGNITUDES ON NON-STATISTICAL TIME-STEPS:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn /= 1d0) .and. (((nn == 2d0) .and. ((n > (nn- 2)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)+ 1d0) .and. &
				(n < (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) .or. &
				((nn > 2d0) .and. ((n >= (nn- 2)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)+ 1d0) .and. &
				(n < (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then

					! ----------------------------------------------------

					! Compute Gravitational field for Southern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) <= 0) then
						do Qind= NqLB(1), NqUB(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == NqLB(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EGmagInterp(1)= abs(yLinInterp(1))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == NqUB(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) < &
										qk4(j)) .and. (qk4(j) < &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EGmagInterp(1)= abs(yLinInterp(1))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) <= &
	                  qk4(j)) .and. (qk4(j) <= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

										EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								! ----------------------------------------------------

							end do
						end do
					end if

					! ----------------------------------------------------

					! Compute Gravitational field for Northern Magnetic Hemisphere
					if (SpecieT(s)%FluxTubeT(f)%QCellT(1)%qGLT(1) > 0) then
						do Qind= NqLB(1), NqUB(1), 1
							do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

								! ----------------------------------------------------

								if (Qind == NqLB(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGLT(1) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if ((Qind /= NqLB(1)) .and. (Qind /= NqUB(1))) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EGmagInterp(1)= abs(yLinInterp(1))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if
								end if

								if (Qind == NqUB(1)) then
									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1) > &
										qk4(j)) .and. (qk4(j) > &
										SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1))) then

										xLinInterp(1)= qk4(j)
										xLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%qGCT(1)
										xLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1)
										yLinInterp1(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind- 1)
										yLinInterp2(1)= SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind)

										call LinInterpsub(yLinInterp(1), xLinInterp(1), yLinInterp1(1), &
											yLinInterp2(1), xLinInterp1(1), xLinInterp2(1))

										EGmagInterp(1)= abs(yLinInterp(1))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
											AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
										end if

										! ----------------------------------------------------

									end if

									if ((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGCT(1) >= &
	                  qk4(j)) .and. (qk4(j) >= &
	                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%qGHT(1))) then

										EGmagInterp(1)= abs(SpecieT(s)%FluxTubeT(f)%EGmagRT(Qind))

										if (EGmagInterp(1) == 0d0) then
											AGmagN(j)= 0d0
										else
		                  AGmagN(j)= abs((SpecieT(s)%qsT(1)/SpecieT(s)%msT(1))*EGmagInterp(1))
										end if

										! ----------------------------------------------------

										! DIAGNOSTIC FLAGS FOR NON-ZERO VALUES:

										if ((AGmagN(j) <= 0) .and. (EGmagInterp(1) /= 0d0)) then
											write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
												' NON-STATISTICAL TIME-STEP HAS NEGATIVE OR ZERO AGmagN= ', AGmagN(j), &
												', EGmagInterp= ', EGmagInterp(1), &
												' VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', &
												f, ', Q GRID CELL= ', Qind, ', TIME-STEP= ', n, &
												', AND PARTICLE= ', j, &
												' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
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
						if ((isnan(real(AGmagN(j))) .eqv. .true.) .or. &
							(size(AGmagN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' AGmagN= ', &
								AGmagN(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
								' IN GRAVITATIONAL FIELD SUBROUTINE' // achar(27) // '[0m.'
						end if
					end do

					! ----------------------------------------------------

				else if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 0) then

					AGmagN(:)= 0d0

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

    !if (SpecieT(s)%FluxTubeT(f)%GRAVflagT(1) == 1) then
     !do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
     	!if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
     	!	(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
     	!	do Qind= NqLB(1), NqUB(1), 1
     	!		if (rank == 0) then
     	!			if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0) then
   						!write(*, *) 'EGmagRT s, f, nn, Qind= ', s, f, nn, Qind, &
   						!	SpecieT(s)%FluxTubeT(f)%EGmagRT(nn, Qind)
     	!			end if
     	!		end if
     	!	end do
     	!end if
     !end do
      !do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
      !  if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
      !   	(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
       		!write(*, *) 'AGmagN s, f= ', s, f, AGmagN(:)
      !  end if
      !end do
    !end if

		! ----------------------------------------------------

	end subroutine GravfieldSub

end module Gravfield
