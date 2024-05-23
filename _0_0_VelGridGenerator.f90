module VelGridGenerator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	0.0 - VELOCITY-SPACE GRID GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! CONSTRUCT PHASE-SPACE GRID:

	subroutine VelGridGeneratorSub

! What we need for Vel Grid Generator output:
!NVperp1GT, NVperp2GT, NVparGT, &
!Vperp1GLT, Vperp1GHT, Vperp1GCT, dVperp1GT, Vperp2GLT, Vperp2GHT, Vperp2GCT, dVperp2GT, VparGLT, VparGHT, VparGCGT, dVparGT, d3vCT
!(VperpGLT, VperpGHT, VperpGCGT, dVperpGT, VparGLT, VparGHT, VparGCGT, dVparGT, d3vCT), &
!NVpGT, NVqGT, NVphiGT, VpGLT, VpGHT, VpGCGT, hVpCT, dVpCT, VqGLGT, VqGHGT, VqGCGT, hVqCT, dVqCT, VphiGLT, VphiGHT, VphiGCGT, hVphiCT, dVphiCT,d33vCT

    ! ---------------------------------------------

    ! SET PRELIMINARY NUMBER OF VELOCITY-SPACE GRID CELLS PER ION PARTICLE SPECIES AND
    ! FLUX TUBE:

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1)= NVperp1GpF
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1)= NVperp2GpF

      end if
      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGpT(1)= NVperpGpF

      end if

      SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)= NVparGpF

    end do

    ! ---------------------------------------------

    ! CREATE PRELIMINARY FIELD-ALIGNED ION VELOCITY SPACE GRID:

    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        ! Set range of Eulerian Vperp1 coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1)/ddVperp1)- 1d0

        ! Set range of Eulerian Vperp2 coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1)/ddVperp2)- 1d0

        ! Positive Vperp1 values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GAT(1)= 0d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GBT(1)= VperpsigmaFac(1)*Vperp12sigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(int(2d0*Vperp12NlinRange- 1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GpT(int(2d0*Vperp12NlinRange- 1)))

        ! Vperp1 values
				SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(int(2d0*Vperp12NlinRange))= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GBT(1)
        do Vperp1ind= 1, int(Vperp12NlinRange), 1
					if (Vperp1ind == 1) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind)= -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GBT(1)
					end if
					if (Vperp1ind == Vperp12NlinRange) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GAT(1)
					end if
					if ((Vperp1ind .gt. 1) .and. (Vperp1ind .lt. int(Vperp12NlinRange))) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind)= &
							-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GBT(1)+ &
							(Vperp1ind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GBT(1)- &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GAT(1))/Vperp12NlinRange)
					end if
				end do
				do Vperp1ind= 1, int(Vperp12NlinRange), 1
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp12NlinRange+ Vperp1ind)= &
						-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp12NlinRange- Vperp1ind)
        end do

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GpT(:)= &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(:)

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1)- 1d0

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1)- 1d0

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperp1GpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperp2GpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1)))

        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperp1GpT(Vperp1ind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind)
        end do
        do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperp2GpT(Vperp2ind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GpT(Vperp2ind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GpT(Vperp2ind)
        end do
      end do
    end if

    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        ! Set range of Eulerian Vperp coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGpT(1)/ddVperp)- 1d0

        ! Positive Vperp values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGAT(1)= 0d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGBT(1)= VperpsigmaFac(1)*Vperpsigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(int(2d0*VperpNlinRange- 1)))

				! Vperp values
				SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(int(2d0*VperpNlinRange))= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGBT(1)
        do Vperpind= 1, int(VperpNlinRange), 1
					if (Vperpind == 1) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind)= -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGBT(1)
					end if
					if (Vperpind == VperpNlinRange) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGAT(1)
					end if
					if ((Vperpind .gt. 1) .and. (Vperpind .lt. int(VperpNlinRange))) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind)= &
							-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGBT(1)+ &
							(Vperpind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGBT(1)- &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGAT(1))/VperpNlinRange)
					end if
				end do
				do Vperpind= 1, int(VperpNlinRange), 1
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(VperpNlinRange+ Vperpind)= &
						-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(VperpNlinRange- Vperpind)
        end do

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1)- 1d0

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperpGpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1)))

        do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVperpGpT(Vperpind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind)
        end do
      end do
    end if

    ! Construct Vpar Grid
    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
      ! Set range of Eulerian Vpar coords.
      SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)= int((SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)/ddVpar)- 1d0)

      if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 1) then
        ! Positive Vpar values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1)= 0d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)= VparsigmaFac(1)*Vparsigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(VparNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp2T(VparNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(int(2d0*VparNlinRange- 1)))

				! Vpar values
				SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(int(2d0*VparNlinRange))= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)
        do Vparind= 1, int(VparNlinRange), 1
					if (Vparind == 1) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)= -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)
					end if
					if (Vparind == VparNlinRange) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1)
					end if
					if ((Vparind .gt. 1) .and. (Vparind .lt. int(VparNlinRange))) then
						SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)= &
							-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)+ &
							(Vparind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)- &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1))/VparNlinRange)
					end if
				end do
				do Vparind= 1, int(VparNlinRange), 1
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(VparNlinRange+ Vparind)= &
						-SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(VparNlinRange- Vparind)
        end do

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)- 1d0

      end if
      if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 0) then
        ! Positive Vpar values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1)= -4d0*Vparsigma
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)= VparsigmaFac(1)*Vparsigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(VparNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(VparNlinRange))

        do Vparind= 1, VparNlinRange, 1
          if (Vparind /= VparNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(Vparind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1)+ &
							(Vparind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)- &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGAT(1))/VparNlinRange) ! Set preliminary range of Vpar-coord
          end if
          if (Vparind == VparNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(Vparind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGBT(1)
          end if
        end do

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)- 1d0

        if (qGA(1) > 0d0) then ! Northern Magnetic Hemisphere
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(:)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(:)
        end if
        if (qGA(1) <= 0d0) then ! Southern Magnetic Hemisphere
          do Vparind= 1, VparNlinRange, 1
            if (Vparind == 1) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)= &
                -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGp1T(VparNlinRange)
            end if
            if (Vparind /= 1) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)= &
                -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VparGp1T(int(VparNlinRange- (Vparind- 1d0)))
            end if
          end do
        end if
      end if

      allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVparGpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

      do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVparGpT(Vparind)= &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)
      end do
    end do

    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then

      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1)))

        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GpT(1), 1
          do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GpT(1), 1
            do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1), 1

              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GT(Vperp1ind, Vperp2ind, Vparind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GpT(Vperp1ind)

              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GT(Vperp1ind, Vperp2ind, Vparind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GpT(Vperp2ind)

              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGT(Vperp1ind, Vperp2ind, Vparind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)

            end do
          end do
        end do
      end do
    end if

    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then

      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGpT(1), 1
          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGpT(1), 1

            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGT(1)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VperpGpT(Vperpind)

            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGT(1)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGpT(Vparind)

          end do
        end do
      end do
    end if

		! ----------------------------------------------------

		! ALLOCATE VELOCITY-SPACE DERIVED DATA TYPES NESTED IN CONFIGURATION-SPACE
		! TYPE NESTED IN FLUX TUBE TYPES NESTED IN PARTICLE SPECIES TYPE:

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

			if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
				! Allocate V2PerpCellT(Vperp1ind, Vperp2ind, Vparind) derived data type nested in
				! FluxTubeT(f) nested in SpecieT(s)
				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
					V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))
			else
				! Allocate VCellT(Vperpind, Vparind) derived data type nested in
				! FluxTubeT(f) nested in SpecieT(s)
				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
					VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))
			end if

			if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
				! Allocate V3CellT(Vpind, Vqind, Vphiind) derived data type nested in
				! FluxTubeT(f) nested in SpecieT(s)
				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
					V3CellT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)))
			end if

		end do

    ! ---------------------------------------------

    ! COMPUTE VELOCITY SPACE GRID BOUNDARIES, CENTER VALUES, METRIC FACTORS, AND VOLUMES:

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then

				allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%d3vCTp( &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), &
					SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1)))

        do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
          do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
            do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

              ! Lower Vperp1 limits of Vperp1 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GT((Vperp1ind- 1)*ddVperp1+ 1, 1, 1)
              ! Lower Vperp2 limits of Vperp2 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GT(1, (Vperp2ind- 1)*ddVperp2+ 1, 1)
              ! Upper Vperp1 limits of Vperp1 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp1GT((Vperp1ind- 1)*ddVperp1+ 1+ ddVperp1, 1, 1)
              ! Upper Vperp2 limits of Vperp2 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%Vperp2GT(1, (Vperp2ind- 1)*ddVperp2+ 1+ ddVperp2, 1)
              ! Center Vperp1 values of Vperp1 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1)= &
                (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(1)+ &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1))/2d0
              ! Center Vperp2 values of Vperp2 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1)= &
                (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1)+ &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1))/2d0

              ! Lower Vpar limits of Vpar cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGT(1, 1, (Vparind- 1)*ddVpar+ 1)
              ! Upper Vpar limits of Vpar cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VparGT(1, 1, (Vparind- 1)*ddVpar+ 1+ ddVpar)
              ! Center Vpar values of Vpar cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCGT(1)= &
                (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1)+ &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1))/2d0

              ! dVperp1 across Vperp1 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1)= &
                abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1)- &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GLT(1))
              ! dVperp2 across Vperp2 cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1)= &
                abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1)- &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1))
              ! dVpar across Vpar cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1)= &
                abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1)- &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1))

              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1)= hVperp1*hVperp2*hVpar* &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1)* &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1)* &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1)

							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%d3vCTp(Vperp1ind, Vperp2ind, Vparind)= &
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1)

              if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1) == 0d0) then
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1)= 1d-9
              end if

              if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1) .lt. 0d0) then
                write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE d3vCT VALUE', &
                 ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
              end if

            end do
          end do
        end do
      end if

      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then

        do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1), 1
          do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

            ! Lower Vperp limits of Vperp cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGLT(1)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT((Vperpind- 1)*ddVperp+ 1, 1)%VperpGT(1)
            ! Upper Vperp limits of Vperp cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGHT(1)= &
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VCellT((Vperpind- 1)*ddVperp+ 1+ ddVperp, 1)%VperpGT(1)
            ! Center Vperp values of Vperp cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGCGT(1)= &
              (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGLT(1)+ &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGHT(1))/2d0

            ! Lower Vpar limits of Vpar cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGLT(1)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(1, (Vparind- 1)*ddVpar+ 1)%VparGT(1)
            ! Upper Vpar limits of Vpar cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGHT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(1, (Vparind- 1)*ddVpar+ 1+ ddVpar)%VparGT(1)
            ! Center Vpar values of Vpar cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCGT(1)= &
              (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGLT(1)+ &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGHT(1))/2d0

            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%hVthetaCT(1)= &
              abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGCGT(1))

            ! dVperp across Vperp cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVperpGT(1)= &
              abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGHT(1)- &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VperpGLT(1))

            ! dVpar across Vpar cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVparGT(1)= &
              abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGHT(1)- &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGLT(1))

            ! dVtheta across Vtheta cells
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVthetaGT(1)= 2d0*pi

            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%d3vCT(1)= &
              hVperp*hVpar*SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%hVthetaCT(1)* &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVperpGT(1)* &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVparGT(1)* &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%dVthetaGT(1)

            if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%d3vCT(1) == 0e0) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%d3vCT(1)= 1d-9
            end if

            if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%d3vCT(1) .lt. 0d0) then
              write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE d3vCT VALUE', &
               ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
            end if
          end do
        end do

      end if
    end do

    ! ---------------------------------------------

    ! SET PRELIMINARY NUMBER OF VELOCITY-SPACE GRID CELLS PER ENA PARTICLE SPECIES AND
    ! FLUX TUBE:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        ! Total number of Vp, Vq, Vphi grid values (Make odd number for 2D distrib fnc. moment integration)
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1)= NVpGpF
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)= NVqGpF
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)= NVphiGpF

      end do
    end if

    ! ---------------------------------------------

    ! CREATE PRELIMINARY DIPOLE ENA VELOCITY SPACE GRID:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        ! Set range of Eulerian Vp coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1)/ddVp)- 1d0

        ! Positive Vp values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGAT(1)= 0d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGBT(1)= VpphisigmaFac*Vpphisigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp1T(VpphiNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp2T(VpphiNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(int(2d0*VpphiNlinRange)))

        ! Postive Vp values
        do Vpind= 1, VpphiNlinRange, 1
          if (Vpind /= VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp1T(Vpind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGAT(1)+ (Vpind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGBT(1)- &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGAT(1))/VpphiNlinRange) ! Set preliminary range of Vp-coord
          end if
          if (Vpind == VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp1T(Vpind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGBT(1)
          end if
        end do

        ! Negative Vp values
        do Vpind= 1, VpphiNlinRange, 1
          if (Vpind == 1) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp2T(Vpind)= &
              -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VpGp1T(VpphiNlinRange)
          end if
          if (Vpind /= 1) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp2T(Vpind)= &
              -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VpGp1T(int(VpphiNlinRange- (Vpind- 1d0)))
          end if
        end do

        ! Put negative and positive Vp values together
        do Vpind= 1, int(2d0*VpphiNlinRange), 1
          if (Vpind .le. VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(Vpind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp2T(Vpind)
          end if
          if (Vpind .gt. VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(Vpind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGp1T(Vpind)
          end if
        end do

        ! Set range of Eulerian Vphi coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)/ddVphi)- 1d0

        ! Positive Vphi values
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGAT(1)= 0d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGBT(1)= VpphisigmaFac*Vpphisigma

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp1T(VpphiNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp2T(VpphiNlinRange))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(int(2d0*VpphiNlinRange)))

        ! Postive Vphi values
        do Vphiind= 1, VpphiNlinRange, 1
          if (Vphiind /= VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp1T(Vphiind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGAT(1)+ &
							(Vphiind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGBT(1)- &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGAT(1))/VpphiNlinRange) ! Set preliminary range of Vphi-coord
          end if
          if (Vphiind == VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp1T(Vphiind)= &
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGBT(1)
          end if
        end do

        ! Negative Vphi values
        do Vphiind= 1, VpphiNlinRange, 1
          if (Vphiind == 1) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp2T(Vphiind)= &
              -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VphiGp1T(VpphiNlinRange)
          end if
          if (Vphiind /= 1) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp2T(Vphiind)= &
              -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VphiGp1T(int(VpphiNlinRange- (Vphiind- 1d0)))
          end if
        end do

        ! Put negative and positive Vphi values together
        do Vphiind= 1, int(2d0*VpphiNlinRange), 1
          if (Vphiind .le. VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(Vphiind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp2T(Vphiind)
          end if
          if (Vphiind .gt. VpphiNlinRange) then
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(Vphiind)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGp1T(Vphiind)
          end if
        end do

        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)- 1d0
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)- 1d0

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVpGpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVphiGpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1)))

        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVpGpT(Vpind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(Vpind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(Vpind)
        end do

        do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVphiGpT(Vphiind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(Vphiind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(Vphiind)
        end do

      end do

      ! Construct Vq Grid
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
        ! Set range of Eulerian Vq coords.
        SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)= (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)/ddVq)- 1d0

        if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 1) then
          ! Positive Vq values
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1)= 0d0
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)= VqsigmaFac*Vqsigma

          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(VqNlinRange))
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp2T(VqNlinRange))
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(int(2d0*VqNlinRange)))

          do Vqind= 1, VqNlinRange, 1
            if (Vqind /= VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(Vqind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1)+ (Vqind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)- &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1))/VqNlinRange) ! Set preliminary range of Vq-coord
            end if
            if (Vqind == VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(Vqind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)
            end if
          end do

          ! Negative Vq values
          do Vqind= 1, VqNlinRange, 1
            if (Vqind == 1) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp2T(Vqind)= &
                -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(VqNlinRange)
            end if
            if (Vqind /= 1) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp2T(Vqind)= &
                -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(int(VqNlinRange- (Vqind- 1d0)))
            end if
          end do

          ! Put negative and positive Vq values together
          do Vqind= 1, int(2d0*VqNlinRange), 1
            if (Vqind .le. VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp2T(Vqind)
            end if
            if (Vqind .gt. VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(Vqind)
            end if
          end do

          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)- 1d0
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)- 1d0

        end if
        if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 0) then
          ! Positive Vq values
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1)= -4d0*Vqsigma
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)= VqsigmaFac*Vqsigma

          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(VqNlinRange))
          allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(VqNlinRange))

          do Vqind= 1, VqNlinRange, 1
            if (Vqind /= VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(Vqind)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1)+ (Vqind- 1d0)*(abs(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)- &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGAT(1))/VqNlinRange) ! Set preliminary range of Vq-coord
            end if
            if (Vqind == VqNlinRange) then
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(Vqind)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGBT(1)
            end if
          end do

          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1)- 1d0
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)= SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)- 1d0

          if (qGA(1) > 0d0) then ! Northern Magnetic Hemisphere
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(:)= &
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(:)
          end if
          if (qGA(1) <= 0d0) then ! Southern Magnetic Hemisphere
            do Vqind= 1, VqNlinRange, 1
              if (Vqind == 1) then
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)= &
                  -SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGp1T(VqNlinRange)
              end if
              if (Vqind /= 1) then
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)= &
                  -SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VqGp1T(int(VqNlinRange- (Vqind- 1d0)))
              end if
            end do
          end if
        end if

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVqGpT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1)))

        do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%dVqGpT(Vqind)= &
            SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind+ 1)- SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)
        end do
      end do

      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)))
        allocate(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGT(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1), &
          SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1), SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1)))

        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGpT(1), 1
          do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGpT(1), 1
            do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGpT(1), 1

              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGT(Vpind, Vqind, Vphiind)= &
                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGpT(Vpind)
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGT(Vpind, Vqind, Vphiind)= &
                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VqGpT(Vqind)
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGT(Vpind, Vqind, Vphiind)= &
                  SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VphiGpT(Vphiind)

            end do
          end do
        end do
      end do
    end if

    ! ---------------------------------------------

    ! COMPUTE ENA VELOCITY SPACE GRID BOUNDARIES, CENTER VALUES, METRIC FACTORS, AND VOLUMES:

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
        do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1), 1
          do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1), 1
            do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1), 1

              ! Lower Vp limits of Vp cells
              SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGLT(1)= &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VpGT((Vpind- 1)*ddVp+ 1, 1, 1)
              ! Upper Vp limits of Vp cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGHT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VpGT((Vpind- 1)*ddVp+ 1+ ddVp, 1, 1)
              ! Center Vp values of Vp cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(1)= &
                (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGLT(1)+ &
                SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGHT(1))/2d0

              ! Lower Vq limits of Vq cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGLGT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VqGT(1, (Vqind- 1)*ddVq+ 1, 1)
              ! Upper Vq limits of Vq cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGHGT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VqGT(1, (Vqind- 1)*ddVq+ 1+ ddVq, 1)
              ! Center Vq values of Vq cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(1)= &
                (SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGLGT(1)+ &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGHGT(1))/2d0

              ! Lower Vphi limits of Vphi cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VphiGT(1, 1, (Vphiind- 1)*ddVphi+ 1)
              ! Upper Vphi limits of Vphi cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%VphiGT(1, 1, (Vphiind- 1)*ddVphi+ 1+ ddVphi)
              ! Center Vphi values of Vphi cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGCGT(1)= &
                (SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(1)+ &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(1))/2d0

              ! Spherical velocity components
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VrGCGT(1)= &
                2d0*SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(1)* &
                cos(thetaGC0(Qind))/sqrt(ellGC0(Qind))+ &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(1)* &
                sin(thetaGC0(Qind))/sqrt(ellGC0(Qind))
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VthetaGCGT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(1)* &
                sin(thetaGC0(Qind))/sqrt(ellGC0(Qind))- &
                2d0*SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(1)* &
                cos(thetaGC0(Qind))/sqrt(ellGC0(Qind))
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VellGCGT(1)= &
                (1d0+ 3d0*((cos(thetaGC0(Qind)))**2d0))

              ! Test dipole velocity components
              VpGCGTest= SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VrGCGT(1)* &
                sin(thetaGC0(Qind))/sqrt(ellGC0(Qind))- &
                2d0*SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VthetaGCGT(1)* &
                cos(thetaGC0(Qind))/sqrt(ellGC0(Qind))
              VqGCGTest= 2d0*SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VrGCGT(1)* &
                cos(thetaGC0(Qind))/sqrt(ellGC0(Qind))+ &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VthetaGCGT(1)* &
                sin(thetaGC0(Qind))/sqrt(ellGC0(Qind))

              if (abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(1)- VpGCGTest) .gt. 1d-6) then
                write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD VpGCGTest VALUE', &
	                ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
              end if

              if (abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(1)- VqGCGTest) .gt. 1d-6) then
                write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD VqGCGTest VALUE', &
	                ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
              end if

              ! ENA Vel-space metric factors
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1)= &
                abs(RE*sin(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VthetaGCGT(1))**3d0/ &
                (sqrt(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VellGCGT(1))))
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1)= &
                abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VrGCGT(1)**3d0/ &
                (RE**2d0*sqrt(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VellGCGT(1))))
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)= &
                abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VrGCGT(1)* &
                sin(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VthetaGCGT(1)))

              ! dVp across Vpar cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1)= &
                abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGHT(1)- &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VpGLT(1))
              ! dVq across Vpar cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVqCT(1)= &
                  abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGHGT(1)- &
                  SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VqGLGT(1))
              ! dVphi across Vpar cells
              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(1)= &
                abs(SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGHT(1)- &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%VphiGLT(1))

              SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1)= &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVpCT(1)* &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVqCT(1)* &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1)* &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1)* &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVqCT(1)* &
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(1)

              if (SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1) == 0d0) then
                SpecieT(s)%FluxtubeT(f)%QCellT(Qind)%V3CellT(Vpind, Vqind, Vphiind)%d33vCT(1)= 1d-9
              end if

            end do
          end do
        end do
      end do
    end if

		! ----------------------------------------------------

		! ION VELOCITY-SPACE GRID DIMENSIONS:

		if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp1GT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperp2GT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do

			! ----------------------------------------------------

		else

			! ----------------------------------------------------

			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperpGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVperpGT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVparGT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		! ENA VELOCITY-SPACE GRID DIMENSIONS:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

				! ----------------------------------------------------

				! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVpGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVpGT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVqGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVqGT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(:)) /= 1) .or. &
					(isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVphiGT(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NVphiGT HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN VELOCITY-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

			end do
		end if

		! ----------------------------------------------------

		! ION EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES:

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
			if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then

				! ----------------------------------------------------

				do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp1GT(1), 1
					do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVperp2GT(1), 1
						do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%NVparGT(1), 1

							! ----------------------------------------------------

							if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCGT(1) == 0d0) then
								SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCGT(1)= 1d-15
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
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GHT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GHT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp1GCT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp1GCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp1GT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GLT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GLT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GHT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GHT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%Vperp2GCT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Vperp2GCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVperp2GT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGLT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGHT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%VparGCGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%d3vCT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3vCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperp1ind= ', Vperp1ind, &
									', Vperp2ind= ', Vperp2ind, ', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
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

						if (SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCGT(1) == 0d0) then
							SpecieT(s)%FluxTubeT(f)%QCellT(Qind)%VCellT(Vperpind, Vparind)%VparGCGT(1)= 1d-15
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
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%dVparGT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%dVparGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVparGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGLT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGLT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGLT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGHT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGHT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGHT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGCGT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VperpGCGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VperpGCGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGLT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGLT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGLT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGHT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGHT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGHT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGCGT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%VparGCGT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VparGCGT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%d3vCT(:)) &
							/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
							VCellT(Vperpind, Vparind)%d3vCT(1))) .eqv. .true.)) then
							write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3vCT HAS', &
								' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
								', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vperpind= ', Vperpind, &
								', AND Vparind= ', Vparind, ' IN VELOCITY-SPACE GRID GENERATOR', &
								' SUBROUTINE' // achar(27) // '[0m.'
						end if

						! ----------------------------------------------------

					end do
				end do

				! ----------------------------------------------------

			end if
		end do

		! ----------------------------------------------------

		! ENA EULERIAN VELOCITY-SPACE GRID LIMITS AND VOLUMES:

		if (SpecieT(1)%FluxTubeT(1)%QEXCHANGEflagT(1) == 1) then
			do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

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
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VpGCGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VpGCGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGLGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGLGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGLGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGHGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGHGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGHGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VqGCGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VqGCGT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									Vphiind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VphiGCGT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%VphiGCGT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
									' VphiGCGT HAS BAD INVERSION, SIZE, OR HAS NaN VALUE', &
									' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', &
									Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%hVphiCT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hVphiCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							if ((size(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%dVpCT(:)) &
								/= 1) .or. (isnan(real(SpecieT(s)%FluxTubeT(f)%QCellT(Qind)% &
								V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1))) .eqv. .true.)) then
								write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dVpCT HAS', &
									' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
									', FLUX TUBE= ', f, ', Qind= ', Qind, ', Vpind= ', Vpind, &
									', Vqind= ', Vqind, ', AND Vphiind= ', Vphiind, &
									' IN VELOCITY-SPACE GRID GENERATOR', &
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
									Vphiind, ' IN VELOCITY-SPACE GRID GENERATOR', &
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
									Vphiind, ' IN VELOCITY-SPACE GRID GENERATOR', &
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
									Vphiind, ' IN VELOCITY-SPACE GRID GENERATOR', &
									' SUBROUTINE' // achar(27) // '[0m.'
							end if

							! ----------------------------------------------------

						end do
					end do
				end do

				! ----------------------------------------------------

			end do
		end if

    ! ---------------------------------------------

  end subroutine VelGridGeneratorSub

end module VelGridGenerator
