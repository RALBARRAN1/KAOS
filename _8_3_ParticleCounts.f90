module ParticleCounts

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.3 PARTICLE COUNTS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MPIReduceSum
use IonParticleCounts
use ENAParticleCounts

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE STATISTICAL GRID COUNTS FOR ALL GRIDS:

	subroutine ParticleCountsSub

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (SpecieT(1)%FluxTubeT(1)%PHASEIONDISTRIBflagT(1) == 1) then
			call IonParticleCountsSub
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			call ENAParticleCountsSub
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE MEAN VELOCITIES AT FINAL TIME:

		! ----------------------------------------------------
		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0) then

								! ----------------------------------------------------

								allocate(Vperp1REF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									Vperp2REF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									VparREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))))

								! ----------------------------------------------------

								jjloop1: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)), 1
									do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
										if (ENAflag(j) .eqv. .false.) then

	                    if (nn == 1) then
	                      ! North Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          <= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          Vperp1REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)
														Vperp2REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop1

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          < SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          Vperp1REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)
														Vperp2REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop1

	                        end if
	                      end if
	                      ! South Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          >= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          Vperp1REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)
														Vperp2REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop1

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          > SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          Vperp1REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)
														Vperp2REF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop1

	                        end if
	                      end if
	                    end if

	                    if (nn /= 1) then
	                      ! North Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												<= qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												Vperp1REF(jj)= Vperp1(j)
														Vperp2REF(jj)= Vperp2(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop1

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												< qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												Vperp1REF(jj)= Vperp1(j)
														Vperp2REF(jj)= Vperp2(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop1

	  											end if
	  										end if
	  										! South Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												>= qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												Vperp1REF(jj)= Vperp1(j)
														Vperp2REF(jj)= Vperp2(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop1

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												> qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												Vperp1REF(jj)= Vperp1(j)
														Vperp2REF(jj)= Vperp2(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop1

	  											end if
	  										end if
	                    end if

										end if
									end do

									! ----------------------------------------------------

									! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

									if ((isnan(Vperp1REF(jj)) .eqv. .true.) .or. (size(Vperp1REF(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' Vperp1REF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(Vperp2REF(jj)) .eqv. .true.) .or. (size(Vperp2REF(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' Vperp2REF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(VparREF(jj)) .eqv. .true.) .or. (size(VparREF(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VparREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do jjloop1

								! ----------------------------------------------------

								SpecieT(s)%FluxTubeT(f)%Vperp1REFpT(nn, Qind)= sum(Vperp1REF(:))
								SpecieT(s)%FluxTubeT(f)%Vperp2REFpT(nn, Qind)= sum(Vperp2REF(:))
								SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)= sum(VparREF(:))

								! ----------------------------------------------------

								deallocate(Vperp1REF, Vperp2REF, VparREF)

								! ----------------------------------------------------

							else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0) then

								SpecieT(s)%FluxTubeT(f)%Vperp1REFpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%Vperp2REFpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)= 0d0

							end if

						end if

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

							if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0) then

								! ----------------------------------------------------

								allocate(VperpREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									VparREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))))

								! ----------------------------------------------------

								jjloop2: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)), 1
									do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
										if (ENAflag(j) .eqv. .false.) then

	                    if (nn == 1) then
	                      ! North Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          <= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          VperpREF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop2

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          < SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          VperpREF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop2

	                        end if
	                      end if
	                      ! South Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          >= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          VperpREF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop2

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          > SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	                          VperpREF(jj)= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)
	                          VparREF(jj)= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)

														cycle jjloop2

	                        end if
	                      end if
	                    end if

	                    if (nn /= 1) then
	                      ! North Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												<= qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												VperpREF(jj)= Vperp(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop2

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												< qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												VperpREF(jj)= Vperp(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop2

	  											end if
	  										end if
	  										! South Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												>= qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												VperpREF(jj)= Vperp(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop2

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												> qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

	  												VperpREF(jj)= Vperp(j)
	  												VparREF(jj)= Vpar(j)

														cycle jjloop2

	  											end if
	  										end if
	                    end if

										end if
									end do

									! ----------------------------------------------------

									! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

									if ((isnan(VperpREF(jj)) .eqv. .true.) .or. (size(VperpREF(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VperpREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(VparREF(jj)) .eqv. .true.) .or. (size(VparREF(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VparREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do jjloop2

								! ----------------------------------------------------

								SpecieT(s)%FluxTubeT(f)%VperpREFpT(nn, Qind)= sum(VperpREF(:))
								SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)= sum(VparREF(:))

								! ----------------------------------------------------

								deallocate(VperpREF, VparREF)

								! ----------------------------------------------------

							else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0) then

								SpecieT(s)%FluxTubeT(f)%VperpREFpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)= 0d0

							end if

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) /= 0) then

							! ----------------------------------------------------

							allocate(VpENAREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))), &
								VqENAREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))), &
								VphiENAREF(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))))

							! ----------------------------------------------------

							jjloop3: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)), 1
								do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
									if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
										(ENAflag(j) .eqv. .true.)) then

  										! North Magnetic Hemisphere
  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												<= qk4(j)) .and. (qk4(j) <= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

  												VpENAREF(jj)= Vp(j)
  												VqENAREF(jj)= Vpar(j)
  												VphiENAREF(jj)= Vphi(j)

													cycle jjloop3

  											end if
  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												< qk4(j)) .and. (qk4(j) <= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

  												VpENAREF(jj)= Vp(j)
  												VqENAREF(jj)= Vpar(j)
  												VphiENAREF(jj)= Vphi(j)

													cycle jjloop3

  											end if
  										end if
  										! South Magnetic Hemisphere
  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												>= qk4(j)) .and. (qk4(j) >= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

  												VpENAREF(jj)= Vp(j)
  												VqENAREF(jj)= Vpar(j)
  												VphiENAREF(jj)= Vphi(j)

													cycle jjloop3

  											end if
  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												> qk4(j)) .and. (qk4(j) >= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

  												VpENAREF(jj)= Vp(j)
  												VqENAREF(jj)= Vpar(j)
  												VphiENAREF(jj)= Vphi(j)

													cycle jjloop3

  											end if
  										end if

									end if
								end do

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(VpENAREF(jj)) .eqv. .true.) .or. (size(VpENAREF(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VpENAREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(VqENAREF(jj)) .eqv. .true.) .or. (size(VqENAREF(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VqENAREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(VphiENAREF(jj)) .eqv. .true.) .or. (size(VphiENAREF(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VphiENAREF HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end do jjloop3

							! ----------------------------------------------------

							SpecieT(s)%FluxTubeT(f)%VpENAREFpT(nn, Qind)= sum(VpENAREF(:))
							SpecieT(s)%FluxTubeT(f)%VqENAREFpT(nn, Qind)= sum(VqENAREF(:))
							SpecieT(s)%FluxTubeT(f)%VphiENAREFpT(nn, Qind)= sum(VphiENAREF(:))

							! ----------------------------------------------------

							deallocate(VpENAREF, VqENAREF, VphiENAREF)

							! ----------------------------------------------------

						else if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) == 0) then

							SpecieT(s)%FluxTubeT(f)%VpENAREFpT(nn, Qind)= 0d0
							SpecieT(s)%FluxTubeT(f)%VqENAREFpT(nn, Qind)= 0d0
							SpecieT(s)%FluxTubeT(f)%VphiENAREFpT(nn, Qind)= 0d0

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							! ----------------------------------------------------

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							Vperp1REFp(1)= SpecieT(s)%FluxTubeT(f)%Vperp1REFpT(nn, Qind)
							MPIRedSumIn(1)= Vperp1REFp(1)
							call MPIReduceSumSub
							Vperp1REFR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							Vperp2REFp(1)= SpecieT(s)%FluxTubeT(f)%Vperp2REFpT(nn, Qind)
							MPIRedSumIn(1)= Vperp2REFp(1)
							call MPIReduceSumSub
							Vperp2REFR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VparREFp(1)= SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)
							MPIRedSumIn(1)= VparREFp(1)
							call MPIReduceSumSub
							VparREFR(1)= MPIRedSumOut(1)

							!call mpi_reduce(Vperp1REFp(1), Vperp1REFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(Vperp2REFp(1), Vperp2REFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(VparREFp(1), VparREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							if (rank == 0) then

								! Mean velocities
								SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind)= &
									Vperp1REFR(1)/(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))
								SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind)= &
									Vperp2REFR(1)/(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))
								SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind)= &
									VparREFR(1)/(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp1REFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp2REFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VparREFRT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VparREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp1REFRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp2REFRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if

						end if

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

							! ----------------------------------------------------

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VperpREFp(1)= SpecieT(s)%FluxTubeT(f)%VperpREFpT(nn, Qind)
							MPIRedSumIn(1)= VperpREFp(1)
							call MPIReduceSumSub
							VperpREFR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VparREFp(1)= SpecieT(s)%FluxTubeT(f)%VparREFpT(nn, Qind)
							MPIRedSumIn(1)= VparREFp(1)
							call MPIReduceSumSub
							VparREFR(1)= MPIRedSumOut(1)

							!call mpi_reduce(VperpREFp(1), VperpREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(VparREFp(1), VparREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							if (rank == 0) then

								! Mean velocities
								SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind)= &
									VperpREFR(1)/(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))
								SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind)= &
									VparREFR(1)/(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VperpREFRT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VperpREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VparREFRT(:, :)) /= &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VparREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VperpREFRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						VpENAREFp(1)= SpecieT(s)%FluxTubeT(f)%VpENAREFpT(nn, Qind)
						MPIRedSumIn(1)= VpENAREFp(1)
						call MPIReduceSumSub
						VpENAREFR(1)= MPIRedSumOut(1)

						VqENAREFp(1)= SpecieT(s)%FluxTubeT(f)%VqENAREFpT(nn, Qind)
						MPIRedSumIn(1)= VqENAREFp(1)
						call MPIReduceSumSub
						VqENAREFR(1)= MPIRedSumOut(1)

						VphiENAREFp(1)= SpecieT(s)%FluxTubeT(f)%VphiENAREFpT(nn, Qind)
						MPIRedSumIn(1)= VphiENAREFp(1)
						call MPIReduceSumSub
						VphiENAREFR(1)= MPIRedSumOut(1)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VpENAREFp(1), VpENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VqENAREFp(1), VqENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VphiENAREFp(1), VphiENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						if (rank == 0) then

							! Mean velocities
							SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind)= &
								VpENAREFR(1)/ &
								(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))
							SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind)= &
								VqENAREFR(1)/ &
								(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))
							SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind)= &
								VphiENAREFR(1)/ &
								(ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))

							! ----------------------------------------------------

							! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VpENAREFRT(:, :)) /= &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VpENAREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VqENAREFRT(:, :)) /= &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VqENAREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(:, :)) /= &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VphiENAREFRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		! BROADCAST MEAN VELOCITIES:

		! ----------------------------------------------------
		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) .and. (rank == 0)) then
			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
				if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
					(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then

					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
					end if
					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then
						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
					end if
					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

						call mpi_barrier(MPI_COMM_WORLD, ierr)
						call mpi_bcast(SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, :), &
							((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1), &
							MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
					end if

				end if
			end do
		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

		! COMPUTE VELOCITY STANDARD DEVIATIONS:

		! ----------------------------------------------------
		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0) then

								! ----------------------------------------------------

								allocate(Vperp1REFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									Vperp2REFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									VparREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))))

								! ----------------------------------------------------

								jjloop4: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)), 1
									do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
										if (ENAflag(j) .eqv. .false.) then

	                    if (nn == 1) then
	                      ! North Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          <= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp10T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp20T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          < SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp10T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp20T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	                        end if
	                      end if
	                      ! South Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          >= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp10T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp20T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          > SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp10T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp20T(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	                        end if
	                      end if
	                    end if

	                    if (nn /= 1) then
	                      ! North Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												<= qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (Vperp1(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (Vperp2(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												< qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (Vperp1(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (Vperp2(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	  											end if
	  										end if
	  										! South Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												>= qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (Vperp1(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (Vperp2(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												> qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														Vperp1REFsig(jj)= (Vperp1(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, Qind))**2d0
														Vperp2REFsig(jj)= (Vperp2(j)- &
															SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop4

	  											end if
	  										end if
	                    end if

										end if
									end do

									! ----------------------------------------------------

									! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

									if ((isnan(Vperp1REFsig(jj)) .eqv. .true.) .or. (size(Vperp1REFsig(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' Vperp1REFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(Vperp2REFsig(jj)) .eqv. .true.) .or. (size(Vperp2REFsig(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' Vperp2REFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(VparREFsig(jj)) .eqv. .true.) .or. (size(VparREFsig(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VparREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do jjloop4

								! ----------------------------------------------------

								SpecieT(s)%FluxTubeT(f)%Vperp1REFsigpT(nn, Qind)= sum(Vperp1REFsig(:))
								SpecieT(s)%FluxTubeT(f)%Vperp2REFsigpT(nn, Qind)= sum(Vperp2REFsig(:))
								SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)= sum(VparREFsig(:))

								! ----------------------------------------------------

								deallocate(Vperp1REFsig, Vperp2REFsig, VparREFsig)

								! ----------------------------------------------------

							else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0) then

								SpecieT(s)%FluxTubeT(f)%Vperp1REFsigpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%Vperp2REFsigpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)= 0d0

							end if

						end if

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

							if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0) then

								! ----------------------------------------------------

								allocate(VperpREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))), &
									VparREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))))

								! ----------------------------------------------------

								jjloop5: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)), 1
									do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
										if (ENAflag(j) .eqv. .false.) then

	                    if (nn == 1) then
	                      ! North Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          <= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          < SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	                        end if
	                      end if
	                      ! South Magnetic Hemisphere
	                      if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	                        if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          >= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	                        end if
	                        if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	                          ((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	                          > SpecieT(s)%FluxTubeT(f)%q0T(j)) &
	                          .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= &
	                          SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vperp0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (SpecieT(s)%FluxTubeT(f)%Vpar0T(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	                        end if
	                      end if
	                    end if

	                    if (nn /= 1) then
	                      ! North Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												<= qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (Vperp(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												< qk4(j)) .and. (qk4(j) <= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (Vperp(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	  											end if
	  										end if
	  										! South Magnetic Hemisphere
	  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
	  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												>= qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (Vperp(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	  											end if
	  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
	  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
	  												> qk4(j)) .and. (qk4(j) >= &
	  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

														VperpREFsig(jj)= (Vperp(j)- &
															SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, Qind))**2d0
														VparREFsig(jj)= (Vpar(j)- &
															SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, Qind))**2d0

														cycle jjloop5

	  											end if
	  										end if
	                    end if

										end if
									end do

									! ----------------------------------------------------

									! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

									if ((isnan(VperpREFsig(jj)) .eqv. .true.) .or. (size(VperpREFsig(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VperpREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									if ((isnan(VparREFsig(jj)) .eqv. .true.) .or. (size(VparREFsig(:)) /= &
										SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind))) then
										write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
											' VparREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
											s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
											' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
											// achar(27) // '[0m.'
									end if

									! ----------------------------------------------------

								end do jjloop5

								! ----------------------------------------------------

								SpecieT(s)%FluxTubeT(f)%VperpREFsigpT(nn, Qind)= sum(VperpREFsig(:))
								SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)= sum(VparREFsig(:))

								! ----------------------------------------------------

								deallocate(VperpREFsig, VparREFsig)

								! ----------------------------------------------------

							else if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0) then

								SpecieT(s)%FluxTubeT(f)%VperpREFsigpT(nn, Qind)= 0d0
								SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)= 0d0

							end if

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) /= 0) then

							! ----------------------------------------------------

							allocate(VpENAREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))), &
								VqENAREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))), &
								VphiENAREFsig(nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))))

							! ----------------------------------------------------

							jjloop6: do jj= 1, nint(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind)), 1
								do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
									if ((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
										(ENAflag(j) .eqv. .true.)) then

  										! North Magnetic Hemisphere
  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) <= 0) then
  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												<= qk4(j)) .and. (qk4(j) <= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

													VpENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind))**2d0
													VqENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind))**2d0
													VphiENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind))**2d0

													cycle jjloop6

  											end if
  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												< qk4(j)) .and. (qk4(j) <= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

													VpENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind))**2d0
													VqENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind))**2d0
													VphiENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind))**2d0

													cycle jjloop6

  											end if
  										end if
  										! South Magnetic Hemisphere
  										if (SpecieT(s)%FluxTubeT(f)%qGLGT(nn, 1) > 0) then
  											if ((Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												>= qk4(j)) .and. (qk4(j) >= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

													VpENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind))**2d0
													VqENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind))**2d0
													VphiENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind))**2d0

													cycle jjloop6

  											end if
  											if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. &
  												((SpecieT(s)%FluxTubeT(f)%qGLGT(nn, Qind) &
  												> qk4(j)) .and. (qk4(j) >= &
  												SpecieT(s)%FluxTubeT(f)%qGHGT(nn, Qind)))) then

													VpENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, Qind))**2d0
													VqENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, Qind))**2d0
													VphiENAREFsig(jj)= (Vp(j)- &
														SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, Qind))**2d0

													cycle jjloop6

  											end if
  										end if

									end if
								end do

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(VpENAREFsig(jj)) .eqv. .true.) .or. (size(VpENAREFsig(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VpENAREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(VqENAREFsig(jj)) .eqv. .true.) .or. (size(VqENAREFsig(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VqENAREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(VphiENAREFsig(jj)) .eqv. .true.) .or. (size(VphiENAREFsig(:)) /= &
									SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VphiENAREFsig HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end do jjloop6

							! ----------------------------------------------------

							SpecieT(s)%FluxTubeT(f)%VpENAREFsigpT(nn, Qind)= sum(VpENAREFsig(:))
							SpecieT(s)%FluxTubeT(f)%VqENAREFsigpT(nn, Qind)= sum(VqENAREFsig(:))
							SpecieT(s)%FluxTubeT(f)%VphiENAREFsigpT(nn, Qind)= sum(VphiENAREFsig(:))

							! ----------------------------------------------------

							deallocate(VpENAREFsig, VqENAREFsig, VphiENAREFsig)

							! ----------------------------------------------------

						else if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qind) == 0) then

							SpecieT(s)%FluxTubeT(f)%VpENAREFsigpT(nn, Qind)= 0d0
							SpecieT(s)%FluxTubeT(f)%VqENAREFsigpT(nn, Qind)= 0d0
							SpecieT(s)%FluxTubeT(f)%VphiENAREFsigpT(nn, Qind)= 0d0

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

							! ----------------------------------------------------

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							Vperp1REFsigp(1)= SpecieT(s)%FluxTubeT(f)%Vperp1REFsigpT(nn, Qind)
							MPIRedSumIn(1)= Vperp1REFsigp(1)
							call MPIReduceSumSub
							Vperp1REFsigR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							Vperp2REFsigp(1)= SpecieT(s)%FluxTubeT(f)%Vperp2REFsigpT(nn, Qind)
							MPIRedSumIn(1)= Vperp2REFsigp(1)
							call MPIReduceSumSub
							Vperp2REFsigR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VparREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)
							MPIRedSumIn(1)= VparREFsigp(1)
							call MPIReduceSumSub
							VparREFsigR(1)= MPIRedSumOut(1)

							!call mpi_reduce(Vperp1REFp(1), Vperp1REFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(Vperp2REFp(1), Vperp2REFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(VparREFp(1), VparREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							if (rank == 0) then

								! Velocity standard deviations
								SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(nn, Qind)= &
									sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*Vperp1REFsigR(1))
								SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(nn, Qind)= &
									sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*Vperp2REFsigR(1))
								SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, Qind)= &
									sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VparREFsigR(1))

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(:, :)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp1REFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(:, :)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp2REFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(:, :)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VparREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp1REFsigRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' Vperp2REFsigRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if

						end if

						if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

							! ----------------------------------------------------

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VperpREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VperpREFsigpT(nn, Qind)
							MPIRedSumIn(1)= VperpREFsigp(1)
							call MPIReduceSumSub
							VperpREFsigR(1)= MPIRedSumOut(1)

							call mpi_barrier(MPI_COMM_WORLD, ierr)
							VparREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VparREFsigpT(nn, Qind)
							MPIRedSumIn(1)= VparREFsigp(1)
							call MPIReduceSumSub
							VparREFsigR(1)= MPIRedSumOut(1)

							!call mpi_reduce(VperpREFp(1), VperpREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							!call mpi_reduce(VparREFp(1), VparREFR(1), 1, &
							!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

							if (rank == 0) then

								! Velocity standard deviations
								SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(nn, Qind)= &
									sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VperpREFsigR(1))
								SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, Qind)= &
									sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqRT(nn)* &
									(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
									SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VparREFsigR(1))

								! ----------------------------------------------------

								! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(:, :)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VperpREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if ((isnan(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, Qind)) .eqv. &
									.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(:, :)) /= &
									(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
									((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VparREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								if (SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(nn, Qind) < 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' VperpREFsigRT HAS NEGATIVE VALUE FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
										' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
										// achar(27) // '[0m.'
								end if

								! ----------------------------------------------------

							end if

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

				if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
					do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

						! ----------------------------------------------------

						VpENAREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VpENAREFsigpT(nn, Qind)
						MPIRedSumIn(1)= VpENAREFsigp(1)
						call MPIReduceSumSub
						VpENAREFsigR(1)= MPIRedSumOut(1)

						VqENAREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VqENAREFsigpT(nn, Qind)
						MPIRedSumIn(1)= VqENAREFsigp(1)
						call MPIReduceSumSub
						VqENAREFsigR(1)= MPIRedSumOut(1)

						VphiENAREFsigp(1)= SpecieT(s)%FluxTubeT(f)%VphiENAREFsigpT(nn, Qind)
						MPIRedSumIn(1)= VphiENAREFsigp(1)
						call MPIReduceSumSub
						VphiENAREFsigR(1)= MPIRedSumOut(1)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VpENAREFp(1), VpENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VqENAREFp(1), VqENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						!call mpi_barrier(MPI_COMM_WORLD, ierr)
						!call mpi_reduce(VphiENAREFp(1), VphiENAREFR(1), 1, &
						!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

						if (rank == 0) then

							! Velocity standard deviations
							SpecieT(s)%FluxTubeT(f)%VpENAREFsigRT(nn, Qind)= &
								sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VpENAREFsigR(1))

							SpecieT(s)%FluxTubeT(f)%VqENAREFsigRT(nn, Qind)= &
								sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VqENAREFsigR(1))

							SpecieT(s)%FluxTubeT(f)%VphiENAREFsigRT(nn, Qind)= &
								sqrt((1d0/((ranksize(1)*nint(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NqENART(nn)* &
								(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/ &
								SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))))- 1d0))*VphiENAREFsigR(1))

							! ----------------------------------------------------

							! DIAGNOSTIC FLAG FOR NaN VALUES AND CONSISTENCY:

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VpENAREFsigRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VpENAREFsigRT(:, :)) /= &
								(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VpENAREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VqENAREFsigRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VqENAREFsigRT(:, :)) /= &
								(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VqENAREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							if ((isnan(SpecieT(s)%FluxTubeT(f)%VphiENAREFsigRT(nn, Qind)) .eqv. &
								.true.) .or. (size(SpecieT(s)%FluxTubeT(f)%VphiENAREFsigRT(:, :)) /= &
								(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
								((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1))) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VphiENAREFsigRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
									s, ', FLUX TUBE= ', f, ', Qind= ', Qind, ', AND STATISTICAL', &
									' TIME-STEP= ', nn, ' IN PARTICLE COUNTS SUBROUTINE' &
									// achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

						end if

						! ----------------------------------------------------

					end do
				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ION L-SHELL DRIFT EXTREMA OVER ENTIRE SIMULATION:

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
				(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then

				pdriftionMax(1)= maxval(pdriftion(:))

				! ----------------------------------------------------

				! REDUCE ION L-SHELL DRIFTS STATISTICS TO MPI ROOT RANK (0):

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(pdriftionMax(1), pdriftionMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				pdriftionMean(1)= sum(pdriftion(:))/SpecieT(s)%FluxTubeT(f)%NsT(1)
				MPIRedSumIn(1)= pdriftionMean(1)
				call MPIReduceSumSub
				pdriftionMeanR(1)= MPIRedSumOut(1)

				!call mpi_reduce(pdriftionMean(1), pdriftionMeanR(1), 1, &
				!	MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

				if (rank == 0) then

					SpecieT(s)%FluxTubeT(f)%pdriftionMaxRT(1)= pdriftionMaxR(1)
					SpecieT(s)%FluxTubeT(f)%pdriftionMeanRT(1)= pdriftionMeanR(1)/ranksize(1)

				end if

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! COMPUTE ION DRIFT EXTREMA OVER ENTIRE SIMULATION:

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

			if (n == SpecieT(s)%FluxTubeT(f)%NtT(1)) then

				qdriftionMax(1)= maxval(qdriftion(:))
				phidriftionMax(1)= maxval(phidriftion(:))
				Vperp1Min(1)= minval(Vperp1N(:))
				Vperp1Max(1)= maxval(Vperp1N(:))
				Vperp2Min(1)= minval(Vperp2N(:))
				Vperp2Max(1)= maxval(Vperp2N(:))
				VparMin(1)= minval(Vpar(:))
				VparMax(1)= maxval(Vpar(:))

				! ----------------------------------------------------

				! REDUCE ALL STATISTICAL ION DRIFTS TO MPI ROOT RANK (0):

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(qdriftionMax(1), qdriftionMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(phidriftionMax(1), phidriftionMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(Vperp1Min(1), Vperp1MinR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(Vperp1Max(1), Vperp1MaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(Vperp2Min(1), Vperp2MinR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(Vperp2Max(1), Vperp2MaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VparMin(1), VparMinR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VparMax(1), VparMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				if (rank == 0) then

					SpecieT(s)%FluxTubeT(f)%qdriftionMaxRT(1)= qdriftionMaxR(1)
					SpecieT(s)%FluxTubeT(f)%phidriftionMaxRT(1)= phidriftionMaxR(1)
					SpecieT(s)%FluxTubeT(f)%Vperp1MinRT(1)= Vperp1MinR(1)
					SpecieT(s)%FluxTubeT(f)%Vperp1MaxRT(1)= Vperp1MaxR(1)
					SpecieT(s)%FluxTubeT(f)%Vperp2MinRT(1)= Vperp2MinR(1)
					SpecieT(s)%FluxTubeT(f)%Vperp2MaxRT(1)= Vperp2MaxR(1)
					SpecieT(s)%FluxTubeT(f)%VparMinRT(1)= VparMinR(1)
					SpecieT(s)%FluxTubeT(f)%VparMaxRT(1)= VparMaxR(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR ALL MAXIMUM VELOCITY VALUES WITHIN VELOCITY-SPACE GRID:

					!if (SpecieT(s)%FluxTubeT(f)%Vperp1MaxRT(1) > SpecieT(s)%FluxTubeT(f)% &
					!	QCellT(1)%V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), &
					!	1, 1)%Vperp1GHT(1)) then
					!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					!		' Vperp1MinRT= ', SpecieT(s)%FluxTubeT(f)%Vperp1MinRT(1), &
					!		' VALUE BELOW LOWER GRID LIMIT AND/OR ', &
					!		' Vperp1MaxRT= ', SpecieT(s)%FluxTubeT(f)%Vperp1MaxRT(1), &
					!		' VALUE ABOVE UPPER GRID LIMIT FOR ', &
					!		' SPECIE= ', s, ', FLUX TUBE= ', f, ' OVER ENTIRE', &
					!		' SIMULATION IN PARTICLE COUNTS SUBROUTINE' &
					!		// achar(27) // '[0m.'
					!end if

					!if (SpecieT(s)%FluxTubeT(f)%Vperp2MaxRT(1) > SpecieT(s)%FluxTubeT(f)% &
					!	QCellT(1)%V2PerpCellT(1, &
					!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1)) then
					!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					!		' Vperp2MinRT= ', SpecieT(s)%FluxTubeT(f)%Vperp2MinRT(1), &
					!		' VALUE BELOW LOWER GRID LIMIT AND/OR ', &
					!		' Vperp2MaxRT= ', SpecieT(s)%FluxTubeT(f)%Vperp2MaxRT(1), &
					!		' VALUE ABOVE UPPER GRID LIMIT FOR ', &
					!		' SPECIE= ', s, ', FLUX TUBE= ', f, ' OVER ENTIRE', &
					!		' SIMULATION IN PARTICLE COUNTS SUBROUTINE' &
					!		// achar(27) // '[0m.'
					!end if

					!if ((SpecieT(s)%FluxTubeT(f)%VparMinRT(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					!	V2PerpCellT(1, 1, 1)%VparGLT(1)) .or. &
					!	(SpecieT(s)%FluxTubeT(f)%VparMaxRT(1) > SpecieT(s)%FluxTubeT(f)% &
					!	QCellT(1)%V2PerpCellT(1, 1, &
					!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
					!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					!		' VparMinRT= ', SpecieT(s)%FluxTubeT(f)%VparMinRT(1), &
					!		' VALUE BELOW LOWER GRID LIMIT AND/OR ', &
					!		' VparMaxRT= ', SpecieT(s)%FluxTubeT(f)%VparMaxRT(1), &
					!		' VALUE ABOVE UPPER GRID LIMIT FOR ', &
					!		' SPECIE= ', s, ', FLUX TUBE= ', f, ' OVER ENTIRE', &
					!		' SIMULATION IN PARTICLE COUNTS SUBROUTINE' &
					!		// achar(27) // '[0m.'
					!end if

					! ----------------------------------------------------

				end if

			end if

		end if

		if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 0) then

			if (n == SpecieT(s)%FluxTubeT(f)%NtT(1)) then

				qdriftionMax(1)= maxval(qdriftion(:))
				phidriftionMax(1)= maxval(phidriftion(:))
				VperpMin(1)= minval(VperpN(:))
				VperpMax(1)= maxval(VperpN(:))
				VparMin(1)= minval(Vpar(:))
				VparMax(1)= maxval(Vpar(:))

				! ----------------------------------------------------

				! REDUCE ALL STATISTICAL ION DRIFTS TO MPI ROOT RANK (0):

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(qdriftionMax(1), qdriftionMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(phidriftionMax(1), phidriftionMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VperpMin(1), VperpMinR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VperpMax(1), VperpMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VparMin(1), VparMinR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierr)

				call mpi_barrier(MPI_COMM_WORLD, ierr)
				call mpi_reduce(VparMax(1), VparMaxR(1), 1, &
					MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

				if (rank == 0) then

					SpecieT(s)%FluxTubeT(f)%qdriftionMaxRT(1)= qdriftionMaxR(1)
					SpecieT(s)%FluxTubeT(f)%phidriftionMaxRT(1)= phidriftionMaxR(1)
					SpecieT(s)%FluxTubeT(f)%VperpMinRT(1)= VperpMinR(1)
					SpecieT(s)%FluxTubeT(f)%VperpMaxRT(1)= VperpMaxR(1)
					SpecieT(s)%FluxTubeT(f)%VparMinRT(1)= VparMinR(1)
					SpecieT(s)%FluxTubeT(f)%VparMaxRT(1)= VparMaxR(1)

					! ----------------------------------------------------

					! DIAGNOSTIC FLAG FOR ALL MAXIMUM VELOCITY VALUES WITHIN VELOCITY-SPACE GRID:

					if (SpecieT(s)%FluxTubeT(f)%VperpMaxRT(1) > SpecieT(s)%FluxTubeT(f)% &
						QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1)%VparGHT(1)) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VperpMinRT= ', SpecieT(s)%FluxTubeT(f)%VperpMinRT(1), &
							' VALUE BELOW LOWER GRID LIMIT AND/OR ', &
							' VperpMaxRT= ', SpecieT(s)%FluxTubeT(f)%VperpMaxRT(1), &
							' VALUE ABOVE UPPER GRID LIMIT FOR ', &
							' SPECIE= ', s, ', FLUX TUBE= ', f, ' OVER ENTIRE', &
							' SIMULATION IN PARTICLE COUNTS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if ((SpecieT(s)%FluxTubeT(f)%VparMinRT(1) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						VCellT(1, 1)%VparGLT(1)) .or. &
						(SpecieT(s)%FluxTubeT(f)%VparMaxRT(1) > SpecieT(s)%FluxTubeT(f)% &
						QCellT(1)%VCellT(1, &
						SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' VparMinRT= ', SpecieT(s)%FluxTubeT(f)%VparMinRT(1), &
							' VALUE BELOW LOWER GRID LIMIT AND/OR ', &
							' VparMaxRT= ', SpecieT(s)%FluxTubeT(f)%VparMaxRT(1), &
							' VALUE ABOVE UPPER GRID LIMIT FOR ', &
							' SPECIE= ', s, ', FLUX TUBE= ', f, ' OVER ENTIRE', &
							' SIMULATION IN PARTICLE COUNTS SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

			end if

		end if

		! ----------------------------------------------------

		! COMPUTE ENA DRIFTS OVER ENTIRE SIMULATION:

		!do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1
		!	if (((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)) .and. &
		!		((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. (ENAflag(j) .eqv. .true.))) then

		!		rdriftENA(j)= abs(rk4(j)- SpecieT(s)%FluxTubeT(f)%r0T(j))
		!		thetadriftENA(j)= abs(thetak4(j)- SpecieT(s)%FluxTubeT(f)%theta0T(j))
		!		phidriftENA(j)= abs(phik4(j)- SpecieT(s)%FluxTubeT(f)%phi0T(j))

		!	end if
		!end do

		! ----------------------------------------------------

		! COMPUTE ENA DRIFT EXTREMA OVER ENTIRE SIMULATION:

		!if ((n == SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
		!	(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)) then
		!	rdriftENAMax(1)= maxval(rdriftENA(:))
		!	thetadriftENAMax(1)= maxval(thetadriftENA(:))
		!	phidriftENAMax(1)= maxval(phidriftENA(:))

			! ----------------------------------------------------

			! REDUCE ALL STATISTICAL ENA DRIFTS STATISTICS TO MPI ROOT RANK (0):

		!	call mpi_barrier(MPI_COMM_WORLD, ierr)
		!	call mpi_reduce(rdriftENAMax(1), rdriftENAMaxR(1), 1, &
		!		MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

		!	call mpi_barrier(MPI_COMM_WORLD, ierr)
		!	call mpi_reduce(thetadriftENAMax(1), thetadriftENAMaxR(1), 1, &
		!		MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

		! call mpi_barrier(MPI_COMM_WORLD, ierr)
		!	call mpi_reduce(phidriftENAMax(1), phidriftENAMaxR(1), 1, &
		!		MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

		!	if (rank == 0) then

		!		SpecieT(s)%FluxTubeT(f)%rdriftENAMaxRT(1)= rdriftENAMaxR(1)
		!		SpecieT(s)%FluxTubeT(f)%thetadriftENAMaxRT(1)= thetadriftENAMaxR(1)
		!		SpecieT(s)%FluxTubeT(f)%phidriftENAMaxRT(1)= phidriftENAMaxR(1)

				! ----------------------------------------------------

				! DIAGNOSTIC FLAG FOR ALL ENA POSITION DRIFTS WITHIN LIMITS:

		!		if (SpecieT(s)%FluxTubeT(f)%rdriftENAMaxRT(1) > &
		!			SpecieT(s)%FluxTubeT(f)%rdriftLimT(1)) then
		!			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		!				' MAXIMUM R DRIFT= ', SpecieT(s)%FluxTubeT(f)%rdriftENAMaxRT(1), &
		!				' EXCEEDS ENA R DRIFT LIMIT= ', SpecieT(s)%FluxTubeT(f)%rdriftLimT(1), &
		!				' OVER ENTIRE SIMULATION FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
		!				' IN PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
		!		end if

		!		if (SpecieT(s)%FluxTubeT(f)%thetadriftENAMaxRT(1) > &
		!			SpecieT(s)%FluxTubeT(f)%thetadriftLimT(1)) then
		!			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		!				' MAXIMUM THETA DRIFT= ', SpecieT(s)%FluxTubeT(f)%thetadriftENAMaxRT(1), &
		!				' EXCEEDS ENA THETA DRIFT LIMIT= ', SpecieT(s)%FluxTubeT(f)%thetadriftLimT(1), &
		!				' OVER ENTIRE SIMULATION FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
		!				' IN PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
		!		end if

		!		if (SpecieT(s)%FluxTubeT(f)%phidriftENAMaxRT(1) > &
		!			SpecieT(s)%FluxTubeT(f)%phidriftLimT(1)) then
		!			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
		!				' MAXIMUM PHI DRIFT= ', SpecieT(s)%FluxTubeT(f)%phidriftENAMaxRT(1), &
		!				' EXCEEDS ENA PHI DRIFT LIMIT= ', SpecieT(s)%FluxTubeT(f)%phidriftLimT(1), &
		!				' OVER ENTIRE SIMULATION FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
		!				' IN PARTICLE COUNTS SUBROUTINE' // achar(27) // '[0m.'
		!		end if

				! ----------------------------------------------------

		!	end if

		!end if

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then
				if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						write(*, *) 'Mean Vperp1 mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vperp1 sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vperp2 mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vperp2 sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vpar mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vpar sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
					else
						write(*, *) 'Mean Vperp mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VperpREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vperp sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vpar mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VparREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vpar sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VparREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
					end if
					if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
						write(*, *) 'Mean Vp mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VpENAREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vp sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VpENAREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vq mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VqENAREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vq sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VqENAREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vphi mean over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
						write(*, *) 'Mean Vphi sigma over all altitude [km/s]= ', &
							(sum(SpecieT(s)%FluxTubeT(f)%VphiENAREFsigRT(nn, :))/(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
								SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1))*1d-3
					end if
				end if
			end if
		end do

		! ----------------------------------------------------

		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
				(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then
				if (rank == 0) then
					if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then
						write(*, *) 'LOWER VPERP1 LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp1GLT(1)*1d-3)
						write(*, *) 'LOWER VPERP2 LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%Vperp2GLT(1)*1d-3)
						write(*, *) 'UPPER VPERP1 LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1)*1d-3)
						write(*, *) 'UPPER VPERP2 LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT( &
							1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1)*1d-3)
						write(*, *) 'LOWER VPAR LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, 1)%VparGLT(1)*1d-3)
						write(*, *) 'UPPER VPAR LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(1, 1, &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1)*1d-3)
						write(*, *) 'FINAL MPI ION MAX L-SHELL DRIFT [RE]= ', &
							SpecieT(1)%FluxTubeT(1)%pdriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MAX Q DRIFT= ', SpecieT(s)%FluxTubeT(f)%qdriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MAX PHI DRIFT [rads]= ', &
							SpecieT(1)%FluxTubeT(1)%phidriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MIN Vperp1 [m/s]= ', SpecieT(s)%FluxTubeT(f)%Vperp1MinRT(1)
						write(*, *) 'FINAL MPI ION MAX Vperp1 [m/s]= ', SpecieT(s)%FluxTubeT(f)%Vperp1MaxRT(1)
						write(*, *) 'FINAL MPI ION MIN Vperp2 [m/s]= ', SpecieT(s)%FluxTubeT(f)%Vperp2MinRT(1)
						write(*, *) 'FINAL MPI ION MAX Vperp2 [m/s]= ', SpecieT(s)%FluxTubeT(f)%Vperp2MaxRT(1)
						write(*, *) 'FINAL MPI ION MIN Vpar [m/s]= ', SpecieT(s)%FluxTubeT(f)%VparMinRT(1)
						write(*, *) 'FINAL MPI ION MAX Vpar [m/s]= ', SpecieT(s)%FluxTubeT(f)%VparMaxRT(1)
					else
						write(*, *) 'LOWER VPERP LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VperpGLT(1)*1d-3)
						write(*, *) 'UPPER VPERP LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT( &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1)%VperpGHT(1)*1d-3)
						write(*, *) 'LOWER VPAR LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VparGLT(1)*1d-3)
						write(*, *) 'UPPER VPAR LIMIT [km/s]= ', &
							nint(SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, &
							SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1)*1d-3)
						write(*, *) 'FINAL MPI ION MAX L-SHELL DRIFT [RE]= ', &
							SpecieT(1)%FluxTubeT(1)%pdriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MAX Q DRIFT= ', SpecieT(s)%FluxTubeT(f)%qdriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MAX PHI DRIFT [rads]= ', &
							SpecieT(1)%FluxTubeT(1)%phidriftionMaxRT(1)
						write(*, *) 'FINAL MPI ION MIN Vperp [m/s]= ', SpecieT(s)%FluxTubeT(f)%VperpMinRT(1)
						write(*, *) 'FINAL MPI ION MAX Vperp [m/s]= ', SpecieT(s)%FluxTubeT(f)%VperpMaxRT(1)
						write(*, *) 'FINAL MPI ION MIN Vpar [m/s]= ', SpecieT(s)%FluxTubeT(f)%VparMinRT(1)
						write(*, *) 'FINAL MPI ION MAX Vpar [m/s]= ', SpecieT(s)%FluxTubeT(f)%VparMaxRT(1)
					end if
				end if
			end if
		end do

		! ----------------------------------------------------

		!do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		!	if ((nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) .and. &
		!		(n == SpecieT(s)%FluxTubeT(f)%NtT(1))) then
		!		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

		!			if (rank == 0) then

		!				write(*, *) 'FINAL MPI ENA MAX R DRIFT [m]= ', &
		!					SpecieT(s)%FluxTubeT(f)%rdriftENAMaxRT(1)
		!				write(*, *) 'FINAL MPI ENA MAX THETA DRIFT [rads]= ', &
		!					SpecieT(s)%FluxTubeT(f)%thetadriftENAMaxRT(1)
		!				write(*, *) 'FINAL MPI ENA MAX PHI DRIFT [rads]= ', &
		!					SpecieT(s)%FluxTubeT(f)%phidriftENAMaxRT(1)

		!			end if

		!		end if
		!	end if
		!end do

		! ----------------------------------------------------
		! ----------------------------------------------------

		if (rank == 0) then
			if (n == 1) then
				call cpu_time(S83End)
				write(S83string, '(i10)')  nint(S83End)
				write(*, *) trim('%% 8.1- RANK= ' // adjustl(rankstring)) // &
					trim(', REAL-TIME= ' // adjustl(S83string)) // &
					trim(' s. INITIAL PARTICLE COUNTS COMPLETE %%')
			end if
		end if

		!if (rank == 0) then
		!	do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		!		if (((n /= 1) .and. (nn /= 1) .and. &
		!			(n == sum(SpecieT(s)%FluxTubeT(f)%ndatfacGT(1:nn- 1))))) then
		!			call cpu_time(S83End)
		!			write(S83string, '(i10)')  nint(S83End)
		!			write(nnstring, '(I5)') nn
		!			write(Timestring, '(I5)') nint(Time(1))
		!			write(pdriftmaxstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%pdriftionMaxRT(1)
		!			write(pdriftmeanstring, '(D10.2)') SpecieT(s)%FluxTubeT(f)%pdriftionMeanRT(1)
		!			write(*, *) trim('** RANK= ' // adjustl(rankstring)) // &
		!				trim(', MAX ITERATIVE L-SHELL DRIFT [RE]= ' // adjustl(pdriftmaxstring)) // &
		!				trim(', MEAN ITERATIVE L-SHELL DRIFT [RE]= ' // adjustl(pdriftmeanstring)) // &
		!				trim(', REAL-TIME= ' // adjustl(S83string)) // &
		!				trim(' s., SIM-TIME= ' // adjustl(Timestring)) // &
		!				trim(' s., STATISTICAL TIME-STEP= ' // adjustl(nnstring))
		!		end if
		!	end do
		!end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine ParticleCountsSub

end module ParticleCounts
