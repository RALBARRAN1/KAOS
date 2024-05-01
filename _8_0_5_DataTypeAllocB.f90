module DataTypeAllocB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.5- DATA TYPE ALLOCATION B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! ALLOCATE KINETIC SOLVER DERIVED DATA TYPES:

	subroutine DataTypeAllocBSub

    ! ----------------------------------------------------

    ! ALLOCATE PARTICLE-INDEPENDENT DERIVED DATA TYPES:

    ! ----------------------------------------------------

		allocate(SpecieT(s)%FluxTubeT(f)%TimeT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
		do Qind= NqLB(1), NqUB(1), 1
			allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%TeNT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
		end do

		! ----------------------------------------------------

		! ALLOCATE BOUNDARY CONDITION VARIABLES AS DERIVED DATA TYPES:

		! ----------------------------------------------------

		allocate(SpecieT(s)%FluxTubeT(f)%NsnT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))
		allocate(SpecieT(s)%FluxTubeT(f)%NsnRRT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))

		allocate(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxIonT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxIonT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))

		if (SpecieT(s)%FluxTubeT(f)%LBREPLENISHflagT(1) == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%NqReNormLBreplenishIonT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))
		end if
		if (SpecieT(s)%FluxTubeT(f)%UBREPLENISHflagT(1) == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%NqReNormUBreplenishIonT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))
		end if

		if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
			allocate(SpecieT(s)%FluxTubeT(f)%NqReNormLBoutfluxENAT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
				SpecieT(s)%FluxTubeT(f)%NqReNormUBoutfluxENAT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
				SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
				SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
				SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
				SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))
		end if

		allocate(SpecieT(s)%FluxTubeT(f)%LBNetDensityT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1), &
			SpecieT(s)%FluxTubeT(f)%UBNetDensityT(SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))

		! ----------------------------------------------------

    ! ALLOCATE KINETIC VARIABLES AS DERIVED DATA TYPES:

		! ----------------------------------------------------

		! ALLOCATE ION PARTICLE COUNTS VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then

      allocate(SpecieT(s)%FluxTubeT(f)%NqReNormT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
        ((NqUB(1)- NqLB(1))+ 1)))

      allocate(SpecieT(s)%FluxTubeT(f)% &
        NqTp(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)))

			if (rank == 0) then
	      allocate(SpecieT(s)%FluxTubeT(f)% &
	        NqRTp(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)))
			end if

      do Qind= NqLB(1), NqUB(1), 1

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
          NqT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

        if (rank == 0) then
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
            NqRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
        end if

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
						N2PerpphReNormT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

					if (rank == 0) then
		        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		          N2PerpphReNormRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

		        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%N2PerpphRTp( &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), &
		          SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

		        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
			          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

			              allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT( &
											Vperp1ind, Vperp2ind, Vparind)%N2PerpphRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

			          end do
							end do
		        end do

					end if

				else

					allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
						NphReNormT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), &
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

					if (rank == 0) then
		        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		          NphReNormRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

		        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NphRTp( &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), &
		          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), &
		          SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

		        do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
		          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1
		              allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		                VCellT(Vperpind, Vparind)%NphRT( &
		                SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
		          end do
		        end do
					end if

				end if

      end do
    end if

    if (SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) then
			if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

	      allocate(SpecieT(s)%FluxTubeT(f)%Vperp1REFpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%Vperp2REFpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)), &
	        SpecieT(s)%FluxTubeT(f)%VparREFpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)))

				allocate(SpecieT(s)%FluxTubeT(f)%Vperp1REFsigpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%Vperp2REFsigpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)), &
	        SpecieT(s)%FluxTubeT(f)%VparREFsigpT(1, &
	        ((NqUB(1)- NqLB(1))+ 1)))

				allocate(SpecieT(s)%FluxTubeT(f)%Vperp1REFRT(1, &
					((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%Vperp2REFRT(1, &
					((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%VparREFRT(1, &
					((NqUB(1)- NqLB(1))+ 1)))

	      if (rank == 0) then

					allocate(SpecieT(s)%FluxTubeT(f)%Vperp1REFsigRT(1, &
	          ((NqUB(1)- NqLB(1))+ 1)), &
						SpecieT(s)%FluxTubeT(f)%Vperp2REFsigRT(1, &
	          ((NqUB(1)- NqLB(1))+ 1)), &
	          SpecieT(s)%FluxTubeT(f)%VparREFsigRT(1, &
	          ((NqUB(1)- NqLB(1))+ 1)))

	      end if

			else

				allocate(SpecieT(s)%FluxTubeT(f)%VperpREFpT(1, &
					((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%VparREFpT(1, &
					((NqUB(1)- NqLB(1))+ 1)))

				allocate(SpecieT(s)%FluxTubeT(f)%VperpREFsigpT(1, &
					((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%VparREFsigpT(1, &
					((NqUB(1)- NqLB(1))+ 1)))

				allocate(SpecieT(s)%FluxTubeT(f)%VperpREFRT(1, &
					((NqUB(1)- NqLB(1))+ 1)), &
					SpecieT(s)%FluxTubeT(f)%VparREFRT(1, &
					((NqUB(1)- NqLB(1))+ 1)))

				if (rank == 0) then

					allocate(SpecieT(s)%FluxTubeT(f)%VperpREFsigRT(1, &
						((NqUB(1)- NqLB(1))+ 1)), &
						SpecieT(s)%FluxTubeT(f)%VparREFsigRT(1, &
						((NqUB(1)- NqLB(1))+ 1)))

				end if

			end if

    end if

    ! ----------------------------------------------------

    ! ALLOCATE ENA PARTICLE COUNTS VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

      allocate(SpecieT(s)%FluxTubeT(f)%NqReNormENAT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
        ((NqUB(1)- NqLB(1))+ 1)))
      allocate(SpecieT(s)%FluxTubeT(f)% &
        NqENATp(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)))
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)% &
          NqENARTp(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)))
      end if

      do Qind= NqLB(1), NqUB(1), 1

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
          NqENAT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
        if (rank == 0) then
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
            NqENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
        end if

				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
          NphReNormENAT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)))

				if (rank == 0) then
					allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
	          NphReNormENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), &
	          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), &
	          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)))

          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NphENARTp( &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), &
            SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
	        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
	          do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
	            do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1
	                allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
	                  V3CellT(Vpind, Vqind, Vphiind)%NphENART( &
	                  SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
	            end do
	          end do
	        end do
				end if
      end do
    end if

    if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONREFflagT(1) == 1) .and. &
			(SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1)) then
      allocate(SpecieT(s)%FluxTubeT(f)%VpENAREFpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VqENAREFpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VphiENAREFpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)))

			allocate(SpecieT(s)%FluxTubeT(f)%VpENAREFsigpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VqENAREFsigpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VphiENAREFsigpT(1, &
        ((NqUB(1)- NqLB(1))+ 1)))

			allocate(SpecieT(s)%FluxTubeT(f)%VpENAREFRT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VqENAREFRT(1, &
        ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%VphiENAREFRT(1, &
        ((NqUB(1)- NqLB(1))+ 1)))

      if (rank == 0) then

				allocate(SpecieT(s)%FluxTubeT(f)%VpENAREFsigRT(1, &
          ((NqUB(1)- NqLB(1))+ 1)), &
          SpecieT(s)%FluxTubeT(f)%VqENAREFsigRT(1, &
          ((NqUB(1)- NqLB(1))+ 1)), &
          SpecieT(s)%FluxTubeT(f)%VphiENAREFsigRT(1, &
          ((NqUB(1)- NqLB(1))+ 1)))

      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE ION DISTRIBUTION FUNCTION VARIABLES IN DERIVED DATA TYPES:

    do Qind= NqLB(1), NqUB(1), 1
      if (SpecieT(s)%FluxTubeT(f)%PHASEIONDISTRIBflagT(1) == 1) then

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

	        if (rank == 0) then
	          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
	            F2PerpphRTp(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
	            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), &
	            SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

		        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
			          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1
		              allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		                V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%F2PerpphRT( &
		                SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			          end do
							end do
		        end do

					end if

				else

					if (rank == 0) then
						allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							FphRTp(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), &
							SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
						do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
							do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1
								allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
									VCellT(Vperpind, Vparind)%FphRT( &
									SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
							end do
						end do
					end if

				end if

      end if
    end do

    ! ----------------------------------------------------

    ! ALLOCATE ENA DISTRIBUTION FUNCTION VARIABLES IN DERIVED DATA TYPES:

    do Qind= NqLB(1), NqUB(1), 1
      if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
        if (rank == 0) then
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
            FphENARTp(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), &
            SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
        end if
        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
          do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
            do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1
              if (rank == 0) then
                allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
                  V3CellT(Vpind, Vqind, Vphiind)%FphENART( &
                  SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
              end if
            end do
          end do
        end do
      end if
    end do

    ! ----------------------------------------------------

    ! ALLOCATE ZEROTH ION MOMENT VARIABLES IN DERIVED DATA TYPES:

		if (SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M0phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    if (SpecieT(s)%FluxTubeT(f)%PHASEDENSITYIONMOMENTflagT(1) == 1) then
      do Qind= NqLB(1), NqUB(1), 1

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

					if (rank == 0) then
		        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
							do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
			          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1
		              allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		                V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%g02PerpphRT( &
		                SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
			          end do
			        end do
						end do
					end if

				else

					if (rank == 0) then
						do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
		          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

		              allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
		                VCellT(Vperpind, Vparind)%g0phRT( &
		                SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))

		          end do
		        end do
					end if

				end if

      end do
    end if

    ! ----------------------------------------------------

    ! ALLOCATE FIRST PERPENDICULAR ION MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEVELPERPIONMOMENTflagT(1) == 1) then
      if (rank == 0) then

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

		      allocate(SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
		        ((NqUB(1)- NqLB(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((NqUB(1)- NqLB(1))+ 1)))

				else

					allocate(SpecieT(s)%FluxTubeT(f)%M1PerpphRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						((NqUB(1)- NqLB(1))+ 1)))

				end if

      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE FIRST PARALLEL ION MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEVELPARIONMOMENTflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M1ParphRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND ION MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYIONMOMENTflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
      if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
				if (rank == 0) then
					allocate(SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						((NqUB(1)- NqLB(1))+ 1)))
				end if
        allocate(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
				allocate(SpecieT(s)%FluxTubeT(f)%nuIonNeutPoiT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
				allocate(SpecieT(s)%FluxTubeT(f)%nuIonNeutRankSumRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND PERPENDICULAR ION MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPERPIONMOMENTflagT(1) == 1) then
      if (rank == 0) then

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

	        allocate(SpecieT(s)%FluxTubeT(f)%M2Perp1phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((NqUB(1)- NqLB(1))+ 1)))
					allocate(SpecieT(s)%FluxTubeT(f)%M2Perp2phRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	          ((NqUB(1)- NqLB(1))+ 1)))

				else

					allocate(SpecieT(s)%FluxTubeT(f)%M2PerpphRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
						((NqUB(1)- NqLB(1))+ 1)))

				end if

      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND PARALLEL ION MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%PHASEENERGYPARIONMOMENTflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2ParphRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))

				if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
					allocate(SpecieT(s)%FluxTubeT(f)%DensityOutputRT((NqUB(1)- NqLB(1))+ 1), &
						SpecieT(s)%FluxTubeT(f)%TemperatureOutputRT((NqUB(1)- NqLB(1))+ 1))
					allocate(SpecieT(s)%FluxTubeT(f)%EAInertialOutputRT((NqUB(1)- NqLB(1))+ 1), &
						SpecieT(s)%FluxTubeT(f)%EAPressureOutputRT((NqUB(1)- NqLB(1))+ 1), &
						SpecieT(s)%FluxTubeT(f)%EAmagOutputRT((NqUB(1)- NqLB(1))+ 1))
				end if

      end if
    end if

		! ----------------------------------------------------

    ! ALLOCATE MOMENT FILTER VARIABLES IN DERIVED DATA TYPES:

		if ((SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1) .and. (rank == 0)) then
			allocate(SpecieT(s)%FluxTubeT(f)%MomentFiltInT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)), &
			SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
			allocate(SpecieT(s)%FluxTubeT(f)%M0FiltAvrgRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
			allocate(SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
			allocate(SpecieT(s)%FluxTubeT(f)%M1Perp2FiltAvrgRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
			allocate(SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
			allocate(SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
				((NqUB(1)- NqLB(1))+ 1)))
		end if

    ! ----------------------------------------------------

    ! ALLOCATE AMBIPOLAR ELECTRIC FIELD VARIABLES IN DERIVED DATA TYPES:

    if (rank == 0) then
      allocate(SpecieT(s)%FluxTubeT(f)%LambdaDRT( &
				SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)), &
				SpecieT(s)%FluxTubeT(f)%EAInertialRT( &
        SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)), &
        SpecieT(s)%FluxTubeT(f)%EAPressureRT( &
        SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, ((NqUB(1)- NqLB(1))+ 1)))
    end if
    allocate(SpecieT(s)%FluxTubeT(f)%EAmagRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
      ((NqUB(1)- NqLB(1))+ 1)))
		allocate(SpecieT(s)%FluxTubeT(f)%EGmagRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
      ((NqUB(1)- NqLB(1))+ 1)))

    ! ----------------------------------------------------

    ! ALLOCATE PARALLEL ELECTRIC FIELD VARIABLES IN DERIVED DATA TYPES:

    if (rank == 0) then
      allocate(SpecieT(s)%FluxTubeT(f)%PhiParRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
	      ((NqUB(1)- NqLB(1))+ 1)))
    end if
    allocate(SpecieT(s)%FluxTubeT(f)%EPmagRT(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
      ((NqUB(1)- NqLB(1))+ 1)))

    ! ----------------------------------------------------

    ! ALLOCATE ZEROTH ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      do Qind= NqLB(1), NqUB(1), 1
        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
          do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
            do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1
              if (rank == 0) then
                allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
                  V3CellT(Vpind, Vqind, Vphiind)%g0phENART( &
                  SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))
              end if
            end do
          end do
        end do
      end do
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M0phENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE FIRST P-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M1PphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE FIRST Q-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M1QphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE FIRST PHI-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M1PHIphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2phENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND P-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2PphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND Q-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2QphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

    ! ----------------------------------------------------

    ! ALLOCATE SECOND PHI-COMPONENT ENA MOMENT VARIABLES IN DERIVED DATA TYPES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (rank == 0) then
        allocate(SpecieT(s)%FluxTubeT(f)%M2PHIphENART(SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, &
          ((NqUB(1)- NqLB(1))+ 1)))
      end if
    end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DataTypeAllocBSub

end module DataTypeAllocB
