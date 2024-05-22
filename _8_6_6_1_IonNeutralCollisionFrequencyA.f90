module IonNeutralCollisionFrequencyA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.6.6.1 ION-NEUTRAL COLLISION FREQUENCY A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use MomentFilter
use GaussianRNG
use PoissonRNG

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE ION-NEUTRAL COLLISION FREQUENCIES:

	subroutine IonNeutralCollisionFrequencyASub

  ! ----------------------------------------------------

  if (rank == 0) then

		! ----------------------------------------------------

		! FILTER ION MOMENTS ACCORDING TO MOVING AVERAGE POINT:

		if (SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1) then

			! Filter Parallel Velocity Moment
			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M1ParMAfilterPt- 1d0)/2d0)) then
				MAfilterPt(1)= M1ParMAfilterPt
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					call MomentFilterSub
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
				end do
			end if

			! Filter Perp1 Velocity Moment
			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M1Perp1MAfilterPt- 1d0)/2d0)) then
				MAfilterPt(1)= M1Perp1MAfilterPt
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					call MomentFilterSub
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
				end do
			end if

			! Filter Perp2 Velocity Moment
			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M1Perp2MAfilterPt- 1d0)/2d0)) then
				MAfilterPt(1)= M1Perp2MAfilterPt
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind)
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					call MomentFilterSub
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%M1Perp2FiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
				end do
			end if

			! Filter Parallel Energy Moment
			if (SpecieT(s)%FluxTubeT(f)%NqICT(1) >= nint((M2ParMAfilterPt- 1d0)/2d0)) then
				MAfilterPt(1)= M2ParMAfilterPt
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, Qind)
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					call MomentFilterSub
				end do
				do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
					SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)
				end do
			end if

		end if

		! ----------------------------------------------------

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1

      ! ----------------------------------------------------

			! COMPUTE ION-NEUTRAL COLLISIONAL CROSS-SECTIONS:
			! Note: Use collisional cross-sections from Lindsay and Stebbings 2005.

			if (SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 0) then
	      if (SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, Qind) /= 0d0) then
	        SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind)= &
	          ((4.07d0- 0.269d0*log((6.242d18)*(1d-3)*SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, Qind)))**2d0)* &
	          ((1d0- exp(-415d0/((6.242d18)*(1d-3)*SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, Qind))))**0.8d0)
				else if (SpecieT(s)%FluxTubeT(f)%M2ParphRT(nn, Qind) == 0d0) then
	        SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind)= 0d0
				end if
      end if

			if (SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1) then
				if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then
		      if (SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, Qind) /= 0d0) then
		        SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind)= &
		          ((4.07d0- 0.269d0*log((6.242d18)*(1d-3)*SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, Qind)))**2d0)* &
		          ((1d0- exp(-415d0/((6.242d18)*(1d-3)*SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, Qind))))**0.8d0)
		      else if (SpecieT(s)%FluxTubeT(f)%M2ParFiltAvrgRT(nn, Qind) == 0d0) then
		        SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind)= 0d0
					end if
				end if
      end if

      ! ----------------------------------------------------

      ! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

      if ((isnan(real(SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind))) &
        .eqv. .true.) .or. (size(sigmaIonNeut(:)) /= 1d0)) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
          ' sigmaIonNeut HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
          s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
          ', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
          ' SUBROUTINE' // achar(27) // '[0m.'
      end if

      ! ----------------------------------------------------

      ! DIAGNOSTIC FLAGS FOR CONSISTENT ION-NEUTRAL COLLISIONAL CROSS-SECTION:

      if (((SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) /= 0d0)) .or. &
        ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) == 0d0))) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NON-ZERO VALUE ', &
          '	M0phRT= ', SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind), &
          ' AND sigmaIonNeut= ', SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind), &
          ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
          ', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
          ' SUBROUTINE' // achar(27) // '[0m.'
      end if

			! ----------------------------------------------------

			! COMPUTE BACKGROUND MAXWELLIAN TEST NEUTRAL VELOCITY:

			! ----------------------------------------------------

      if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) then

				! Compute Gaussian random variates
				UniformRN1(1)= VxNeutrandn1(1)
				UniformRN2(1)= VxNeutrandn2(1)
				call GaussianRNGSub
				VxNeutp(1)= GaussianRN(1)

				UniformRN1(1)= VyNeutrandn1(1)
				UniformRN2(1)= VyNeutrandn2(1)
				call GaussianRNGSub
				VyNeutp(1)= GaussianRN(1)

				UniformRN1(1)= VzNeutrandn1(1)
				UniformRN2(1)= VzNeutrandn2(1)
				call GaussianRNGSub
				VzNeutp(1)= GaussianRN(1)

				! MB distribution standard deviations (Global Cartesian coords.)
  			sigmaVxNeut(1)= sqrt(kB*TNeut(1)/(mNeut)) ! Where TparX= TparY= TparZ= TNeut
  			sigmaVyNeut(1)= sigmaVxNeut(1)
  			sigmaVzNeut(1)= sigmaVxNeut(1)

				! MB Cartesian velocity components
  			VxNeut(1)= sigmaVxNeut(1)*VxNeutp(1)
				VyNeut(1)= sigmaVyNeut(1)*VyNeutp(1)
				VzNeut(1)= sigmaVzNeut(1)*VzNeutp(1)

  			! ----------------------------------------------------

        ! Ion population translational (parallel) Cartesian velocity components [m/s]
        !MionVx(1)= (cos(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))/ &
        !  (sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))))* &
        !  (3d0*SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
        !  cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
        !  sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind)))
        !MionVy(1)= (sin(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))/ &
        !  (sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))))* &
        !  (3d0*SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
        !  cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
        !  sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind)))
        !MionVz(1)= SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
        !  ((3d0*((cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind)))**2d0)- 1d0)/ &
        !  sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind)))

				if (SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 0) then
					MionVx(1)= cos(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))* &
						(3d0*SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)* &
						(1d0- 3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))- &
						SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind)* &
						sin(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))
					MionVy(1)= sin(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))* &
						(3d0*SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)* &
						(1d0- 3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp2phRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))
					MionVz(1)= SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)* &
						(3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)- 1d0)/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))+ &
						3d0*SpecieT(s)%FluxTubeT(f)%M1Perp1phRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))
				end if
				if (SpecieT(s)%FluxTubeT(f)%MOMENTFILTERflagT(1) == 1) then
					MionVx(1)= cos(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))* &
						(3d0*SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(nn, Qind)* &
						(1d0- 3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))- &
						SpecieT(s)%FluxTubeT(f)%M1Perp2FiltAvrgRT(nn, Qind)* &
						sin(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))
					MionVy(1)= sin(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))* &
						(3d0*SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(nn, Qind)* &
						(1d0- 3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))+ &
						SpecieT(s)%FluxTubeT(f)%M1Perp2FiltAvrgRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%phiGCGT(nn, Qind))
					MionVz(1)= SpecieT(s)%FluxTubeT(f)%M1ParFiltAvrgRT(nn, Qind)* &
						(3d0*(cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))**2d0)- 1d0)/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))+ &
						3d0*SpecieT(s)%FluxTubeT(f)%M1Perp1FiltAvrgRT(nn, Qind)* &
						cos(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))* &
						sin(SpecieT(s)%FluxTubeT(f)%thetaGCGT(nn, Qind))/ &
						sqrt(SpecieT(s)%FluxTubeT(f)%ellGCGT(nn, Qind))
				end if

  			! Ion-neutral relative translational velocities [m/s]
  			VelIonNeut(1)= abs(sqrt(abs(MionVx(1)- VxNeut(1))**2d0+ abs(MionVy(1)- VyNeut(1))**2d0+ &
  				abs(MionVz(1)- VzNeut(1))**2d0)) ! Relative O+ O velocities

  			! ----------------------------------------------------

  			! COMPUTE NUMBER ION ION-NEUTRAL COLLISIONS PER TIME-STEP
  			! AS POISSON DISTRIBUTED EVENTS (BY ACCEPTANCE-REJECTION TECHNIQUE):

  			! ----------------------------------------------------

  			if (SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) /= 0d0) then

  				! Mean number of ion-neutral collisions per time-step interval (unitless)
					nuIonNeutM(1)= SpecieT(s)%FluxTubeT(f)%nsnormCNeutGT(nn, Qind)* &
  					SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind)*VelIonNeut(1)*SpecieT(s)%hT

  				! Poisson distributed number of total ion-neutral collisions per time-step interval (unitless)
					MeanPoisson(1)= nuIonNeutM(1)
					call PoissonRNGSub
  				SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)= RandPoisson(1)

					! ----------------------------------------------------

					! Ensure number of collisions does not exceed number of available ions in config-grid cell
					if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > &
						SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
						(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) then
							SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)= SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
							(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))
					end if

					! ----------------------------------------------------

  			end if
  			if (SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) == 0d0) then

  				SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)= -2d0

  			end if

  			! ----------------------------------------------------

  			! ENSURE COMPUTATIONAL TIME-STEP RESOLVES ION-NEUTRAL COLLISION FREQUENCIES:

  			!if ((1d0/(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)) > 1d-10) .and. &
        !  (SpecieT(s)%hT >= 1d0/(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)))) then
  			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' COMPUTATIONAL ', &
  			!		' TIME-STEP IS GREATER OR EQUAL TO 1/nuIonNeutRT= ', &
  			!		1d0/(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
  			!		', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
  			!		' IN ION-NEUTRAL COLLISION FREQUENCY A SUBROUTINE' // achar(27) // '[0m.'
  			!end if

        ! ----------------------------------------------------

      else

        SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind)= -2d0

      end if

      ! ----------------------------------------------------

      ! DIAGNOSTIC FLAGS FOR CONSISTENT ION-NEUTRAL COLLISIONAL FREQUENCY:

			if (((SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) == 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) /= -2d0)) .or. &
        ((SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind) /= 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) == -2d0))) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NON-ZERO VALUE', &
          '	sigmaIonNeut= ', SpecieT(s)%FluxTubeT(f)%sigmaIonNeutRT(nn, Qind), &
          ' AND nuIonNeutRT= ', SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind), &
          ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
          ', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
          ' SUBROUTINE' // achar(27) // '[0m.'
      end if

      if (((SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) /= -2d0)) .or. &
        ((SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) /= 0d0) .and. &
        (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) == -2d0))) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' INCONSISTENT NON-ZERO VALUE', &
          '	M0phRT= ', SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind), &
          ' AND nuIonNeutRT= ', SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind), &
          ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
          ', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
          ' SUBROUTINE' // achar(27) // '[0m.'
      end if

			! ----------------------------------------------------

    end do
  end if

	! ----------------------------------------------------

  ! DIAGNOSTIC FLAGS FOR CONSISTENT TOTAL NUMBER OF ION-NEUTRAL COLLISIONS:

	do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
		if (rank == 0) then
			if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > &
				SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
				(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NUMBER OF TOTAL ION-NEUTRAL COLLISIONS', &
					' PER TIME-STEP INTERVAL nuIonNeutRT= ', SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind), &
					' IS GREATER THAN NUMBER OF PARTICLES= ', &
					SpecieT(s)%FluxTubeT(f)%NqRTp(nn, Qind)* &
					(SpecieT(s)%FluxTubeT(f)%d3xCGT(nn, Qind)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)), &
					' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
					', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

  		if (SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)) then
  			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NUMBER OF ION-NEUTRAL COLLISIONS', &
  				' PER TIME-STEP INTERVAL nuIonNeutRT= ', SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind), &
  				' IS GREATER THAN TOTAL PARTICLE NUMBER OVER ALL RANKS= ', SpecieT(s)%FluxTubeT(f)%NsnRRT(nn), &
  				' FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
  				', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
  				' SUBROUTINE' // achar(27) // '[0m.'
  		end if
    end if
  end do

	! ----------------------------------------------------

	! BROADCAST POISSON DISTRIBUTED NUMBER OF ION-NEUTRAL COLLISIONS PER TIME-STEP INTERVAL ACROSS ALL MPI RANKS:

	call mpi_barrier(MPI_COMM_WORLD, ierr)
	call mpi_bcast(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, :), &
		(((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)), &
		MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

	! ----------------------------------------------------

  ! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

  do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1), 1
    if ((isnan(real(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(nn, Qind))) &
      .eqv. .true.) .or. (size(SpecieT(s)%FluxTubeT(f)%nuIonNeutRT(:, :)) /= &
      (SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1)* &
      (((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1)))) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
        ' nuIonNeutRT HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
        s, ', FLUX TUBE= ', f, ', Qind= ', Qind, &
        ', AND STATISTICAL TIME-STEP= ', nn, ' IN ION NEUTRAL COLLISION FREQUENCY A', &
        ' SUBROUTINE' // achar(27) // '[0m.'
    end if
  end do

  ! ----------------------------------------------------

	end subroutine IonNeutralCollisionFrequencyASub

end module IonNeutralCollisionFrequencyA
