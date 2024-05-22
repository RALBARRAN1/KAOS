module ConfigGridGenerator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	0 - CONFIGURATION-SPACE GRID GENERATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use GridDipolePolynomialSolver

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! CONSTRUCT CONFIGURATION-SPACE GRID:

	subroutine ConfigGridGeneratorSub

		! What we need for Grid Generator output:
		! qGC0T, hqC0T, dpC0T, dqC0T, dphiC0T, rGC0T, phiGC0T, thetaGC0T, ellGC0T, qGL0T, qGH0T, pGC0T, d3xC0T, TsPerp0T, TsPar0T, Ts0T, &
		! NqG0T, NqGT, Te0T

		! ----------------------------------------------------

		if ((INITIALGRIDflag == 1) .and. (nn /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' INCONSISTENT INITIAL GRID FLAG FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((INITIALGRIDflag == 0) .and. (nn == 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' INCONSISTENT INITIAL GRID FLAG FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ---------------------------------------------

		! SET CONFIGURATION-SPACE GRID DIMENSIONS:

		SpecieT(s)%FluxTubeT(f)%NqGpT(1)= &
			(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4)
		SpecieT(s)%FluxTubeT(f)%NqICT(1)= &
			abs(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1))+ 1d0 ! Q cell range for density initialization
		SpecieT(s)%FluxTubeT(f)%NqGTp(1)= (SpecieT(s)%FluxTubeT(f)%NqGpT(1)/dq)- 1d0 ! Number of config cells
		SpecieT(s)%FluxTubeT(f)%NqG0T(1)= SpecieT(s)%FluxTubeT(f)%NqGTp(1)
		SpecieT(s)%FluxTubeT(f)%NqGT(1)= SpecieT(s)%FluxTubeT(f)%NqGTp(1)- 2d0

		! ---------------------------------------------

		! DIAGNOSTIC FLAGS FOR CORRECT CONFIGURATION-SPACE GRID DIMENSIONS:

		if ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
			SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4d0) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' (NqUBT- NqLBT+ 4d0)= ', &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 4d0), &
				', NqGpT= ', SpecieT(s)%FluxTubeT(f)%NqGpT(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((SpecieT(s)%FluxTubeT(f)%NqUBT(1)- &
			SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3d0) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' (NqUBT- NqLBT+ 3d0)= ', &
				(SpecieT(s)%FluxTubeT(f)%NqUBT(1)- SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 3d0), &
				', NqG0T= ', SpecieT(s)%FluxTubeT(f)%NqG0T(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqLBT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqLBT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqLBT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqUBT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqUBT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqUBT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqGpT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGpT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqGpT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqICT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqGTp(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGTp(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqGTp HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

	  ! ----------------------------------------------------

    ! SET PRELIMINARY VELOCITY-SPACE GRID PARAMETERS:

    ! Set Vperp parameters
    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
      NVperp1GpF= NVperpGpF
      NVperp2GpF= NVperpGpF
      Vperp12NlinRange= (NVperp1GpF)/2d0 ! Number of grid cells spanning linear grid
    end if
    if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then
      VperpNlinRange= (NVperpGpF)/2d0 ! Number of grid cells spanning linear grid
    end if

    ! Set Vpar parameters
    if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 1) then
        VparNlinRange= (NVparGpF)/2d0
    end if
    if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 0) then
        VparNlinRange= NVparGpF
    end if

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 1) then
        NVpGpF= NVperp1GpF ! Number of Vp Grid Cells (even (div by 2 odd) for +/- log10 Vp grid) (+ 3)
        NVqGpF= NVparGpF ! Number of Vq Grid Cells (even (div by 2 odd) for +/- log10 Vq grid) (+ 3)
        NVphiGpF= NVperp2GpF ! Number of Vphi Grid Cells (even (div by 2 odd) for +/- log10 Vphi grid) (+ 3)
      end if

      if (SpecieT(s)%FluxTubeT(f)%ION2VPERPflagT(1) == 0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NO ENA 1D VPERP GRID CONSTRUCTED', &
          ' IN CONFIG. GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

      VpphisigmaFac= VperpsigmaFac(1) ! Sigma factor with linear grid to resolve thermal core of MB distrib.
      VqsigmaFac= VparsigmaFac(1)

      VpphiNlinRange= (NVpGpF)/2d0
      if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 1) then
          VqNlinRange= (NVqGpF)/2d0
      end if
      if (SpecieT(s)%FluxTubeT(f)%SYMVPARflagT(1) == 0) then
          VqNlinRange= NVqGpF
      end if
    end if

    ! Set Velocity space grid spacing by MB sigma value for initial temperature only.
		! FIXME update velocity grid with new base altitude Ti, Te. in time

		TsPerp= Ti(1) ! Initialize as isotropic (in pitch-angle), where Ti= (1/3)*TsPar+ (2/3)*TsPerp
    TsPar= Ti(1)
    Vperp12sigma= sqrt(kB*TsPerp/SpecieT(s)%msT(1)) ! MB sigma for Vperp12
    Vperpsigma= Vperp12sigma
    Vparsigma= sqrt(kB*TsPar/SpecieT(s)%msT(1)) ! MB sigma for Vpar

    if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
      mNeut= SpecieT(s)%msT(1) ! O neutral mass [kg]
      Vpphisigma= sqrt(kB*TNeut(1)/mNeut) ! MB sigma for Vp, Vphi
      Vqsigma= sqrt(kB*TNeut(1)/mNeut) ! MB sigma for Vq
    end if

		! ---------------------------------------------

    SpecieT(s)%Qindns0GT(1)= 1d0

		! ---------------------------------------------

		! DIAGNOSTIC FLAGS FOR NAN VALUES:

		if ((isnan(real(SpecieT(s)%Qindns0GT(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%Qindns0GT(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0GT= ', &
				SpecieT(s)%Qindns0GT(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ---------------------------------------------

		! SET CONFIGURATION-SPACE GRID TEMPERATURES:

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

			! Note: Isotropic temperature from grid data (Ti= Tperp= Tpar= Tperp1= Tperp2)
      TsPerp0(Qind)= TsPerp
      TsPar0(Qind)= TsPar
			Ts0(Qind)= (1d0/3d0)*(TsPar0(Qind))+ (2d0/3d0)*(TsPerp0(Qind))
      Te0(Qind)= Te(1)

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((size(TsPerp0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(TsPerp0(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerp0 HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(TsPar0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(TsPar0(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPar0 HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(Ts0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(Ts0(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Ts0 HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(Te0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(Te0(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0 HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (Te0(Qind) > dNTeEND) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' INITIAL Te0T IS GREATER THAN Te CAP FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

		! ---------------------------------------------

		if (INITIALGRIDflag == 1) then
	    allocate(qGp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
			allocate(pGp(1))
			allocate(phiGp(1))
	    allocate(rGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(thetaGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(xGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(yGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
	    allocate(zGridOutp(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
		end if

		! ---------------------------------------------

		! SET CONFIGURATION-SPACE LSHELL AND INVARIANT LONGITUDE:

		pGp(1)= Lshell(1) ! Set grid L-shell
		phiGp(1)= phiLshell(1) ! Set grid longitude

		! ---------------------------------------------

		! SET PRELIMINARY FIELD-ALIGNED GRID VALUES:

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1
      if (Qind /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
				if (qGA(1) <= 0d0) then ! SMH
	        qGp(Qind)= qGA(1)+ (Qind- 1d0)*(abs(qGB(1)- qGA(1))/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
				if (qGA(1) > 0d0) then ! NMH
	        qGp(Qind)= qGA(1)- (Qind- 1d0)*(abs(qGB(1)- qGA(1))/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
      end if
      if (Qind == SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
        qGp(Qind)= qGB(1)
      end if

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(qGp(Qind)) .eqv. .true.) .or. (size(qGp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGp= ', &
					qGp(Qind), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

		! ---------------------------------------------

		! DIAGNOSTIC FLAGS FOR NAN VALUES:

		if ((isnan(pGp(1)) .eqv. .true.) .or. (size(pGp(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGp= ', &
				pGp(1), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(phiGp(1)) .eqv. .true.) .or. (size(phiGp(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGp= ', &
				pGp(1), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		! ---------------------------------------------

		! COMPUTE SPHERICAL AND CARTESIAN VALUES OF PRELIMINARY GRID:

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1

      pGridIn(1)= pGp(1)
      qGridIn(1)= qGp(Qind)
			phiGridIn(1)= phiGp(1)

			! ---------------------------------------------

      call GridDipolePolynomialSolverSub

			! ---------------------------------------------

      rGridOutp(Qind)= rGridOut(1)
      thetaGridOutp(Qind)= thetaGridOut(1)
      xGridOutp(Qind)= xGridOut(1)
      yGridOutp(Qind)= yGridOut(1)
      zGridOutp(Qind)= zGridOut(1)

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(rGridOutp(Qind))) .eqv. .true.) .or. &
				(size(rGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGridOutp= ', &
					rGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(thetaGridOutp(Qind))) .eqv. .true.) .or. &
				(size(thetaGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGridOutp= ', &
					thetaGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(xGridOutp(Qind))) .eqv. .true.) .or. &
				(size(xGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGridOutp= ', &
					xGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(yGridOutp(Qind))) .eqv. .true.) .or. &
				(size(yGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGridOutp= ', &
					yGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(zGridOutp(Qind))) .eqv. .true.) .or. &
				(size(zGridOutp(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGridOutp= ', &
					zGridOutp(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    ! ---------------------------------------------

    ! SET CONFIGURATION SPACE GRID LIMITS:

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

      qGL0(Qind)= qGp((Qind- 1)*dq+ 1) ! lower q limits of FA cells
      qGH0(Qind)= qGp((Qind- 1)*dq+ 1+ dq) ! upper q limits of FA cells
      qGC0(Qind)= (qGL0(Qind)+ qGH0(Qind))/2d0 ! center q values of FA cells

      pGC0(Qind)= pGp(1) ! center p values of FA cells
			phiGC0(Qind)= phiGp(1) ! center phi values of FA cells

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(qGL0(Qind))) .eqv. .true.) .or. &
				(size(qGL0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGL0= ', qGL0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(qGH0(Qind))) .eqv. .true.) .or. &
				(size(qGH0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGH0= ', qGH0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(qGC0(Qind))) .eqv. .true.) .or. &
				(size(qGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGC0= ', qGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(pGC0(Qind))) .eqv. .true.) .or. &
				(size(pGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGC0= ', pGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(phiGC0(Qind))) .eqv. .true.) .or. &
				(size(phiGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGC0= ', phiGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if (phiGC0(Qind) /= phiGC0(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGC0= ', phiGC0(Qind), &
					' DOES NOT EQUAL FLUX-TUBE FOOT PRINT phiCG0= ', phiGC0(1), ' FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if (pGC0(Qind) /= pGC0(1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGC0= ', pGC0(Qind), &
					' DOES NOT EQUAL FLUX-TUBE FOOT PRINT phiCG0= ', pGC0(1), ' FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

		! ---------------------------------------------

		! COMPUTE GRID CENTER (r, theta) VALUES, METRIC FACTORS, AND CONFIGURATION-SPACE VOLUMES:

    do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

      pGridIn(1)= pGC0(Qind)
      qGridIn(1)= qGC0(Qind)
			phiGridIn(1)= phiGC0(Qind)

			! ---------------------------------------------

      call GridDipolePolynomialSolverSub

			! ---------------------------------------------

      ! SET CONFIGURATION SPACE GRID CENTER VALUES, METRIC FACTORS, AND VOLUMES:

      rGC0(Qind)= rGridOut(1)
      thetaGC0(Qind)= thetaGridOut(1)

      ellGC0(Qind)= (1d0+ 3d0*cos(thetaGC0(Qind))**2d0)

      hqC0(Qind)= abs(rGC0(Qind)**3d0/(RE**2d0*sqrt(ellGC0(Qind))))
      hpC0(Qind)= abs(RE*sin(thetaGC0(Qind))**3d0/(sqrt(ellGC0(Qind))))
      hphiC0(Qind)= abs(rGC0(Qind)*sin(thetaGC0(Qind)))

      dqC0(Qind)= abs(qGH0(Qind)- qGL0(Qind)) ! dq across FA cells
      dpC0(Qind)= 1d-3 ! dp across FA cells [RE]
      dphiC0(Qind)= pi/180d0 ! dphi across FA cells [1 deg= pi/180 rads]

      d3xC0(Qind)= hqC0(Qind)*hpC0(Qind)*hphiC0(Qind)* &
        dqC0(Qind)*dpC0(Qind)*dphiC0(Qind)

      if (d3xC0(Qind) == 0d0) then
        d3xC0(Qind)= 1d-9
      end if

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR GRID VALUES:

			if (d3xC0(Qind) < 0d0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
          ' NEGATIVE d3xC0T VALUE IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

			if ((isnan(real(rGC0(Qind))) .eqv. .true.) .or. &
				(size(rGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGC0= ', rGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(thetaGC0(Qind))) .eqv. .true.) .or. &
				(size(thetaGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGC0= ', thetaGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(ellGC0(Qind))) .eqv. .true.) .or. &
				(size(ellGC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGC0= ', ellGC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(hqC0(Qind))) .eqv. .true.) .or. &
				(size(hqC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqC0= ', hqC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(hpC0(Qind))) .eqv. .true.) .or. &
				(size(hpC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hpC0= ', hpC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(hphiC0(Qind))) .eqv. .true.) .or. &
				(size(hphiC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hphiC0= ', hphiC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(dqC0(Qind))) .eqv. .true.) .or. &
				(size(dqC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqC0= ', dqC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(dpC0(Qind))) .eqv. .true.) .or. &
				(size(dpC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpC0= ', dpC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(dphiC0(Qind))) .eqv. .true.) .or. &
				(size(dphiC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiC0= ', dphiC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(d3xC0(Qind))) .eqv. .true.) .or. &
				(size(d3xC0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xC0= ', d3xC0(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if (d3xC0(Qind) < 0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE ', &
					' d3xCGT FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

		! ----------------------------------------------------

		! SET STATISTICAL TIME-STEP WITH NEW GRID:

		! Injection time-step [s] (must be > 2 and excludes initial time-step)
		! Mean transit time through LB ghost cell

		SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)= hqC0(SpecieT(s)%FluxTubeT(f)%NqLBT(1))* &
			dqC0(SpecieT(s)%FluxTubeT(f)%NqLBT(1))*(abs(sqrt(SpecieT(s)%msT(1)/ &
			(kB*TsPar0(SpecieT(s)%FluxTubeT(f)%NqLBT(1))))))/SpecieT(s)%hT

		! ----------------------------------------------------

  end subroutine ConfigGridGeneratorSub

end module ConfigGridGenerator
