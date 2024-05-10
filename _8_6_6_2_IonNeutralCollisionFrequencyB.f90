module IonNeutralCollisionFrequencyB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.6.6.2 ION-NEUTRAL COLLISION FREQUENCY B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE ION-NEUTRAL COLLISION FREQUENCIES:

	subroutine IonNeutralCollisionFrequencyBSub

	! ----------------------------------------------------

	! COMPUTE POISSON DISTRIBUTED NUMBER OF ION-NEUTRAL COLLISIONS PER RANK:

	! ----------------------------------------------------

	rrloop: do rr= 0, ranksize(1)- 1, 1

    ! ----------------------------------------------------

		if ((rank == 0) .and. (rr == rank)) then

			! ----------------------------------------------------

      do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
        if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
          (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
          do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

            ! ----------------------------------------------------

        		if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > 0d0) then

    					if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0d0) then
    						if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) <= SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
									!SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)

									write(*, *) 'test= ', SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)

									SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)

									!SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0

								end if
    						if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)
    						end if

    						SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)

    					end if
    					if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0d0) then
    						SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0
    						SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= 0d0
    					end if

    				else

              SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0
              SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= 0d0

            end if

						! ----------------------------------------------------

  					! DIAGNOSTIC FLAGS FOR CORRECT TOTAL NUMBER OF ION-NEUTRAL COLLISIONS PER RANK:

            if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > 0d0) then
							if (SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
    						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
    							' NUMBER OF RANK ION-NEUTRAL COLLISIONS= ', &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind), ' IS LARGER THAN ', &
    							' NUMBER OF AVAILABLE PARTICLES= ', &
    							SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind), ' FOR SPECIE= ', &
    							s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
    							', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY B', &
    							' SUBROUTINE' // achar(27) // '[0m.'
    					end if
            end if

						! ----------------------------------------------------

          end do
        end if
      end do

      ! ----------------------------------------------------

      ! SEND TOTAL NUMBER OF ION-NEUTRAL COLLISIONS TO PROCEEDING RANK:

      ! ----------------------------------------------------

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		    if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
		      (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

		      call mpi_send(SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, :), &
						(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
						SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)), MPI_DOUBLE_PRECISION, rank+ 1, 0, MPI_COMM_WORLD, ierr)

				end if
			end do

			cycle rrloop

      ! ----------------------------------------------------

    end if

		! ----------------------------------------------------

		if ((rank > 0) .and. (rr == rank)) then

      ! ----------------------------------------------------

      ! RECEIVE TOTAL NUMBER OF ION-NEUTRAL COLLISIONS FROM PREVIOUS RANK:

      ! ----------------------------------------------------

			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
		    if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
		      (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

		      call mpi_recv(SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, :), &
						(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
						SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)), MPI_DOUBLE_PRECISION, rank- 1, 0, MPI_COMM_WORLD, status, ierr)

				end if
			end do

      ! ----------------------------------------------------

      do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
        if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
          (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then
          do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

            ! ----------------------------------------------------

            if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > 0d0) then

    					if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) /= 0d0) then
    						if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind) == &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)) then
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0

    							exit rrloop

    						end if
    						if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind) < &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)) then
    							if ((SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)- &
    								SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)) > &
    								SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
    								SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= &
    									SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)
    							end if
    							if ((SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)- &
    								SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)) <= &
    								SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
    								SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= &
    									SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)- &
    									SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)
    							end if
    						end if

    						SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)+ &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)

    					end if
    					if (SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind) == 0d0) then
    						SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0
    						SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)
    					end if

            else

              SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind)= 0d0
              SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind)= 0d0

            end if

            ! ----------------------------------------------------

  					! DIAGNOSTIC FLAGS FOR CORRECT TOTAL NUMBER OF ION-NEUTRAL COLLISIONS PER RANK:

            if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > 0d0) then
    					if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)) then
    						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
    							' NUMBER OF ION-NEUTRAL COLLISIONS OVER ALL RANKS= ', &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, Qind), ' IS LARGER THAN ', &
    							' TOTAL NUMBER OF ION-NEUTRAL COLLISIONS= ', &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind), ' FOR SPECIE= ', &
    							s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
    							', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY B', &
    							' SUBROUTINE' // achar(27) // '[0m.'
    					end if
							if (SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind)) then
    						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
    							' NUMBER OF RANK ION-NEUTRAL COLLISIONS= ', &
    							SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(nn, Qind), ' IS LARGER THAN ', &
    							' NUMBER OF AVAILABLE PARTICLES= ', &
    							SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qind), ' FOR SPECIE= ', &
    							s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
    							', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY B', &
    							' SUBROUTINE' // achar(27) // '[0m.'
    					end if
            end if

            ! ----------------------------------------------------

          end do
        end if
      end do

			! ----------------------------------------------------

			! SEND TOTAL NUMBER OF ION-NEUTRAL COLLISIONS TO PROCEEDING RANK:

			! ----------------------------------------------------

			if (rr < ranksize(1)- 1) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
			    if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
			      (n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(nn)))) then

						call mpi_send(SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(nn, :), &
							(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)), &
							MPI_DOUBLE_PRECISION, rank+ 1, 0, MPI_COMM_WORLD, ierr)

					end if
				end do

				cycle rrloop

			end if
			if (rr >= ranksize(1)- 1) then
				exit rrloop
			end if

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

	end do rrloop

  ! ----------------------------------------------------

	end subroutine IonNeutralCollisionFrequencyBSub

end module IonNeutralCollisionFrequencyB
