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

! CONSTRUCT PHASE-SPACE GRID:

	subroutine ConfigGridGeneratorSub

		! What we need for Grid Generator output:
		! qGC0T, hqC0T, dpC0T, dqC0T, dphiC0T, rGC0T, phiGC0T, thetaGC0T, ellGC0T, qGL0T, qGH0T, pGC0T, d3xC0T, TsPerp0T, TsPar0T, Ts0T, &
		! NqG0T, NqGT, Te0T

		! ----------------------------------------------------

		! SET GRID RANGE ALONG FLUX TUBES:

		if (s == 1) then
			! Starting and ending Q cell for density initialization
			! (including selected cells and does not include lower boundary ghost cell)

			NqLB(1)= NqICA
			NqUB(1)= NqICB

			NqGpF= (NqUB(1)- NqLB(1)+ 4)

		end if

		NqIC(1)= abs(NqICB- NqICA)+ 1d0 ! Q cell range for density initialization

		SpecieT(s)%FluxTubeT(f)%NqICAT= NqICA ! Create nested derived data types
		SpecieT(s)%FluxTubeT(f)%NqICBT= NqICB
		SpecieT(s)%FluxTubeT(f)%NqICT= NqIC

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

		if ((NqICA /= SpecieT(s)%FluxTubeT(f)%NqICAT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqICAT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICAT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICAT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((NqICB /= SpecieT(s)%FluxTubeT(f)%NqICBT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqICBT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICBT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICBT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((NqIC(1) /= SpecieT(s)%FluxTubeT(f)%NqICT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqICT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%NqICT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqICT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR', &
				' SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! ALLOCATE INITIALIZATION CONFIG-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE
		! TYPES NESTED IN PARTICLE SPECIES TYPE:

		! Allocate QCellICT(Qind) derived data type nested in FluxTubeT(f)
		! nested in SpecieT(s)
		allocate(SpecieT(s)%FluxTubeT(f)%QCellICT(SpecieT(s)%FluxTubeT(f)%NqICT(1)))

		! ----------------------------------------------------

		! SET PARTICLE SPECIES INITIAL DENSITY PARAMETERS:

		! ----------------------------------------------------

		SpecieT(s)%FluxTubeT(f)%zns0T= zns0 ! Create nested derived data types
		SpecieT(s)%FluxTubeT(f)%ns0T= ns0
		SpecieT(s)%FluxTubeT(f)%nsnormfacT= nsnormfac

		SpecieT(s)%FluxTubeT(f)%zns0NeutT= zns0Neut ! Create nested derived data types
		SpecieT(s)%FluxTubeT(f)%ns0NeutT= ns0Neut

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

		if ((zns0 /= SpecieT(s)%FluxTubeT(f)%zns0T(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%zns0T(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%zns0T(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zns0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN SIMULATION', &
				' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((ns0 /= SpecieT(s)%FluxTubeT(f)%ns0T(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%ns0T(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%ns0T(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ns0T HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN SIMULATION', &
				' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((nsnormfac /= SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%nsnormfacT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormfacT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN SIMULATION', &
				' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((zns0Neut /= SpecieT(s)%FluxTubeT(f)%zns0NeutT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%zns0NeutT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%zns0NeutT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zns0NeutT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN SIMULATION', &
				' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		if ((ns0Neut /= SpecieT(s)%FluxTubeT(f)%ns0NeutT(1)) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%ns0NeutT(:)) /= 1) .or. &
			(isnan(real(SpecieT(s)%FluxTubeT(f)%ns0NeutT(1))) .eqv. .true.)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ns0NeutT HAS', &
				' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
				' AND FLUX TUBE= ', f, ' IN SIMULATION', &
				' PARAMETERIZATION SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! ALLOCATE CONFIGURATION-SPACE DERIVED DATA TYPES NESTED IN FLUX TUBE TYPES NESTED IN
		! PARTICLE SPECIES TYPE:

		! Allocate QCellT(Qind) derived data type nested in FluxTubeT(f)
		! nested in SpecieT(s)
		allocate(SpecieT(s)%FluxTubeT(f)%QCellT(((NqUB(1)- NqLB(1))+ 1)))
		allocate(SpecieT(s)%FluxTubeT(f)%QCell0T(((NqUB(1)- NqLB(1))+ 3)))

	  ! ----------------------------------------------------

    ! SET PRELIMINARY NUMBER OF CONFIGURATION-SPACE GRID CELLS PER PARTICLE SPECIES AND
    ! FLUX TUBE:

    if (qGA <= 0d0) then ! Southern Magnetic Hemisphere
      SMagHemFlag= 1
    end if
    if (qGA > 0d0) then ! Northern Magnetic Hemisphere
      SMagHemFlag= 0
    end if

    ! Set Vperp parameters
    if (ION2VPERPflag == 1) then
      NVperp1GpF= NVperpGpF
      NVperp2GpF= NVperpGpF
      Vperp12sigmaFac= VperpsigmaFac
      Vperp12NlinRange= (NVperp1GpF)/2d0 ! Number of grid cells spanning linear grid
    end if
    if (ION2VPERPflag == 0) then
      VperpNlinRange= (NVperpGpF)/2d0 ! Number of grid cells spanning linear grid
    end if

    ! Set Vpar parameters
    if (SYMVPARflag == 1) then
        VparNlinRange= (NVparGpF)/2d0
    end if
    if (SYMVPARflag == 0) then
        VparNlinRange= NVparGpF
    end if

    if (QEXCHANGEflag == 1) then
      if (ION2VPERPflag == 1) then
        NVpGpF= NVperp1GpF ! Number of Vp Grid Cells (even (div by 2 odd) for +/- log10 Vp grid) (+ 3)
        NVqGpF= NVparGpF ! Number of Vq Grid Cells (even (div by 2 odd) for +/- log10 Vq grid) (+ 3)
        NVphiGpF= NVperp2GpF ! Number of Vphi Grid Cells (even (div by 2 odd) for +/- log10 Vphi grid) (+ 3)
      end if

      if (ION2VPERPflag == 0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NO ENA 1D VPERP GRID CONSTRUCTED', &
          ' IN CONFIG. GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

      if (ICRflag == 0) then
        VpphisigmaFac= Vperp12sigmaFac ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac

        VpphiNlinRange= (NVpGpF)/2d0
        if (SYMVPARflag == 1) then
            VqNlinRange= (NVqGpF)/2d0
        end if
        if (SYMVPARflag == 0) then
            VqNlinRange= NVqGpF
        end if
      end if
      if (ICRflag == 1) then
        VpphisigmaFac= Vperp12sigmaFac ! Sigma factor with linear grid to resolve thermal core of MB distrib.
        VqsigmaFac= VparsigmaFac

        VpphiNlinRange= (NVpGpF)/2d0
        if (SYMVPARflag == 1) then
            VqNlinRange= (NVqGpF)/2d0
        end if
        if (SYMVPARflag == 0) then
            VqNlinRange= NVqGpF
        end if
      end if
    end if

    ! Set Velocity space grid spacing by MB sigma value for initial temperature only. FIXME update velocity grid with new base altitude Ti, Te. in time
    TsPerp= Ti ! Initialize as isotropic (in pitch-angle), where Ti= (1/3)*TsPar+ (2/3)*TsPerp
    TsPar= Ti
    Vperp12sigma= sqrt(kB*TsPerp/mO) ! MB sigma for Vperp12
    Vperpsigma= Vperp12sigma
    Vparsigma= sqrt(kB*TsPar/mO) ! MB sigma for Vpar

    if (QEXCHANGEflag == 1) then
      mNeut= mO ! O neutral mass [kg]
      Vpphisigma= sqrt(kB*TNeut/mNeut) ! MB sigma for Vp, Vphi
      Vqsigma= sqrt(kB*TNeut/mNeut) ! MB sigma for Vq
    end if

    ! ----------------------------------------------------

    ! Total number of q grid values (Make even number to avoid equatorial
    ! grid boundary):
    if ((s == 1) .and. (f .lt. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 1) .and. (f .ge. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 2) .and. (f .lt. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if
    if ((s == 2) .and. (f .ge. SpecieT(s)%NfT(1)/2d0)) then
      SpecieT(s)%FluxTubeT(f)%NqGpT(1)= NqGpF
    end if

		! ---------------------------------------------
		! ENSURE CORRECT CONFIGURATION-SPACE GRID CELL DIMENSION:

		if ((NqUB(1)- NqLB(1)+ 4d0) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' (NqUB(1)- NqLB(1)+ 4d0)= ', (NqUB(1)- NqLB(1)+ 3d0), &
				', NqGpT= ', SpecieT(s)%FluxTubeT(f)%NqGpT(1), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ---------------------------------------------

    ! SET MAGNETIC L-SHELL VALUES AND PRELIMINARY FIELD-ALIGNED GRID CELLS FOR EACH PARTICLE SPECIES
    ! AND FLUX-TUBE:

    SpecieT(s)%FluxTubeT(f)%NqGT(1)= (SpecieT(s)%FluxTubeT(f)%NqGpT(1)/dq)- 1d0 ! Number of config cells

    allocate(SpecieT(s)%FluxTubeT(f)%qGT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
		allocate(SpecieT(s)%FluxTubeT(f)%pGT(1))
    allocate(SpecieT(s)%FluxTubeT(f)%phiGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGridOutT(SpecieT(s)%FluxTubeT(f)%NqGpT(1)))

		SpecieT(s)%FluxTubeT(f)%pGT(1)= Lshell ! Set grid L-shell

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1
      if (Qind /= SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
				if (qGA <= 0d0) then ! NMH
	        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGA+ (Qind- 1d0)*(abs(qGB- qGA)/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
				if (qGA > 0d0) then ! SMH
	        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGA- (Qind- 1d0)*(abs(qGB- qGA)/SpecieT(s)%FluxTubeT(f)%NqGpT(1))
				end if
      end if
      if (Qind == SpecieT(s)%FluxTubeT(f)%NqGpT(1)) then
        SpecieT(s)%FluxTubeT(f)%qGT(Qind)= qGB
      end if
    end do

    do Qind= 1, SpecieT(s)%FluxTubeT(f)%NqGpT(1), 1

      pGridIn(1)= SpecieT(s)%FluxTubeT(f)%pGT(1)
      qGridIn(1)= SpecieT(s)%FluxTubeT(f)%qGT(Qind)

      call GridDipolePolynomialSolverSub

      SpecieT(s)%FluxTubeT(f)%rGridOutT(Qind)= rGridOut(1)
      SpecieT(s)%FluxTubeT(f)%thetaGridOutT(Qind)= thetaGridOut(1)
			SpecieT(s)%FluxTubeT(f)%phiGridOutT(Qind)= phiGridOut(1)
      SpecieT(s)%FluxTubeT(f)%xGridOutT(Qind)= xGridOut(1)
      SpecieT(s)%FluxTubeT(f)%yGridOutT(Qind)= yGridOut(1)
      SpecieT(s)%FluxTubeT(f)%zGridOutT(Qind)= zGridOut(1)

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%rGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%rGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%rGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%thetaGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%thetaGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%thetaGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%phiGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%phiGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%phiGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%xGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%yGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGridOutT(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGridOutT(:)) /= SpecieT(s)%FluxTubeT(f)%NqGpT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGridOutT= ', &
					SpecieT(s)%FluxTubeT(f)%zGridOutT(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    ! ---------------------------------------------

    ! SET CONFIGURATION SPACE GRID LIMITS:

    allocate(SpecieT(s)%FluxTubeT(f)%phiGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%phiGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%rGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%thetaGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%xGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%yGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%zGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%ellGL0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%ellGH0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%hpC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))
    allocate(SpecieT(s)%FluxTubeT(f)%hphiC0T(SpecieT(s)%FluxTubeT(f)%NqGT(1)))

    do Qind= NqLB(1), NqUB(1)+ 2, 1

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)= &
        SpecieT(s)%FluxTubeT(f)%qGT((Qind- 1)*dq+ 1) ! lower q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)= &
        SpecieT(s)%FluxTubeT(f)%qGT((Qind- 1)*dq+ 1+ dq) ! upper q limits of FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)= &
        (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)+ &
	      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))/2d0 ! center q values of FA cells

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)= SpecieT(s)%FluxTubeT(f)%pGT(1) ! center p values of FA cells

      SpecieT(s)%FluxTubeT(f)%phiGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%phiGridOutT((Qind- 1)*dq+ 1) ! lower phid limits of FA cells
      SpecieT(s)%FluxTubeT(f)%phiGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%phiGridOutT((Qind- 1)*dq+ 1+ dq) ! upper phid limits of FA cells

      SpecieT(s)%FluxTubeT(f)%rGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%rGridOutT((Qind- 1)*dq+ 1) ! lower r limits of FA cells
      SpecieT(s)%FluxTubeT(f)%rGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%rGridOutT((Qind- 1)*dq+ 1+ dq) ! upper r limits of FA cells

      SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%thetaGridOutT((Qind- 1)*dq+ 1) ! lower theta limits of FA cells
      SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%thetaGridOutT((Qind- 1)*dq+ 1+ dq) ! upper theta limits of FA cells

      SpecieT(s)%FluxTubeT(f)%yGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%yGridOutT((Qind- 1)*dq+ 1) ! lower y limits of FA cells
      SpecieT(s)%FluxTubeT(f)%yGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%yGridOutT((Qind- 1)*dq+ 1+ dq) ! upper y limits of FA cells

      SpecieT(s)%FluxTubeT(f)%zGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%zGridOutT((Qind- 1)* dq+ 1) ! lower z limits of FA cells
      SpecieT(s)%FluxTubeT(f)%zGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%zGridOutT((Qind- 1)*dq+ 1+ dq) ! upper z limits of FA cells

			SpecieT(s)%FluxTubeT(f)%xGL0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%xGridOutT((Qind- 1)* dq+ 1) ! lower x limits of FA cells
      SpecieT(s)%FluxTubeT(f)%xGH0T(Qind)= &
        SpecieT(s)%FluxTubeT(f)%xGridOutT((Qind- 1)*dq+ 1+ dq) ! upper x limits of FA cells

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%phiGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%phiGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%phiGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%phiGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%phiGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%phiGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%rGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%rGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%rGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%rGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%rGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%rGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%thetaGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%thetaGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    do Qind= NqLB(1), NqUB(1)+ 2, 1

      pGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1)
      qGridIn(1)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1)

      call GridDipolePolynomialSolverSub

      ! SET CONFIGURATION SPACE GRID CENTER VALUES, METRIC FACTORS, AND VOLUMES:

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)= rGridOut(1)
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1)= thetaGridOut(1)
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1)= phiGridOut(1)
      SpecieT(s)%FluxTubeT(f)%xGC0T(Qind)= xGridOut(1)
      SpecieT(s)%FluxTubeT(f)%yGC0T(Qind)= yGridOut(1)
      SpecieT(s)%FluxTubeT(f)%zGC0T(Qind)= zGridOut(1)

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))**2d0)
      SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%thetaGL0T(Qind))**2d0)
      SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind)= &
        (1d0+ 3d0*cos(SpecieT(s)%FluxTubeT(f)%thetaGH0T(Qind))**2d0)

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)**3d0/ &
        (RE**2d0*sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))) ! q metric factor centered in FA cells
      SpecieT(s)%FluxTubeT(f)%hpC0T(Qind)= &
        abs(RE*sin(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))**3d0/ &
        (sqrt(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1)))) ! p metric factor centered in FA cells
      SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1)* &
        sin(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) ! phi metric factor centered in FA cells

      ! Note: Make dpC equal to the Lshell drift limit over entire simulation and dphidC equal to the
      ! phid drift limit in kinetic solver f90 subroutine.

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)= &
        abs(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1)- SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1)) ! dq across FA cells
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1)= 1d-3 ! dp across FA cells [RE]
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1)= pi/180d0 ! dphi across FA cells [1 deg= pi/180 rads]

      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1)= &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%hpC0T(Qind)* &
        SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1)* &
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1)

      if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) == 0d0) then
        SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1)= 1d-9
      end if

      if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) < 0d0) then
        write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
          ' NEGATIVE d3xC0T VALUE IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
      end if

			! ---------------------------------------------

			! DIAGNOSTIC FLAGS FOR NAN VALUES:

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%xGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%xGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' xGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%xGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%yGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%yGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' yGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%yGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%zGC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%zGC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' zGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%zGC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%ellGL0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGL0T= ', &
					SpecieT(s)%FluxTubeT(f)%ellGL0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%ellGH0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGH0T= ', &
					SpecieT(s)%FluxTubeT(f)%ellGH0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%hpC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%hpC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hpC0T= ', &
					SpecieT(s)%FluxTubeT(f)%hpC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%hphiC0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqGT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hphiC0T= ', &
					SpecieT(s)%FluxTubeT(f)%hphiC0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xC0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			! ---------------------------------------------

    end do

    SpecieT(s)%Qindns0T(1)= 1d0

    ! ---------------------------------------------

    ! Re-index to account for lower and upper ghost cells (non-computational domain)
    SpecieT(s)%FluxTubeT(f)%NqG0T(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)
    SpecieT(s)%FluxTubeT(f)%NqGT(1)= SpecieT(s)%FluxTubeT(f)%NqGT(1)- 2d0

		allocate(SpecieT(s)%FluxTubeT(f)%Te0T(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))

    do Qind= NqLB(1), NqUB(1)+ 2, 1
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1)= TsPerp
      SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1)= TsPar
      SpecieT(s)%FluxTubeT(f)%Te0T(Qind)= Te
    end do

    if (SpecieT(s)%FluxTubeT(f)%NqG0T(1) /= (NqUB(1)- NqLB(1))+ 3) then
      write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' BAD NqG0T VALUE', &
        ' IN GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
    end if

		! ---------------------------------------------

		! DIAGNOSTIC FLAGS FOR NAN VALUES:

		if ((isnan(real(SpecieT(s)%Qindns0T(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%Qindns0T(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Qindns0T= ', &
				SpecieT(s)%Qindns0T(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(SpecieT(s)%FluxTubeT(f)%NqG0T(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqG0T(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqG0T= ', &
				SpecieT(s)%FluxTubeT(f)%NqG0T(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		if ((isnan(real(SpecieT(s)%FluxTubeT(f)%NqGT(1))) .eqv. .true.) .or. &
			(size(SpecieT(s)%FluxTubeT(f)%NqGT(:)) /= 1)) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NqGT= ', &
				SpecieT(s)%FluxTubeT(f)%NqGT(1), &
				' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
				', FLUX TUBE= ', f, ' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
				// achar(27) // '[0m.'
		end if

		do Qind= NqLB(1), NqUB(1)+ 2, 1

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerp0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPar0T= ', &
					SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

			if ((isnan(real(SpecieT(s)%FluxTubeT(f)%Te0T(Qind))) .eqv. .true.) .or. &
				(size(SpecieT(s)%FluxTubeT(f)%Te0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0T= ', &
					SpecieT(s)%FluxTubeT(f)%Te0T(Qind), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' &
					// achar(27) // '[0m.'
			end if

		end do

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

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

		do Qind= NqLB(1), NqUB(1)+ 2, 1
			if ((size(SpecieT(s)%FluxTubeT(f)%Te0T(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%Te0T(Qind))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' Te0T HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%Te0T(Qind) > dNTeEND) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' INITIAL Te0T IS GREATER THAN Te CAP FOR SPECIE= ', s, &
					' AND FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if
		end do

		! ----------------------------------------------------

		! COMPUTE TOTAL FIELD-LINE ARC LENGTH:

		allocate(SpecieT(s)%FluxTubeT(f)%ellqCT(SpecieT(s)%FluxTubeT(f)%NqG0T(1)))
		SpecieT(s)%FluxTubeT(f)%ellqCT(:)= SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1)*SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1)
		SpecieT(s)%FluxTubeT(f)%SUMellqCT(1)= sum(SpecieT(s)%FluxTubeT(f)%ellqCT(:))

		! ----------------------------------------------------

		do Qind= NqLB(1), NqUB(1)+ 2, 1

			! ----------------------------------------------------

			! Total initial ion temperature [K]
			! Note: Isotropic temperature from grid data (Ti= Tperp= Tpar= Tperp1= Tperp2)
			SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%Ts0T(1)= (1d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))+ &
				(2d0/3d0)*(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%hqC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' hqCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dpC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dpCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dqC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dqCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%dphiC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' dphiCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPerp0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsPerpT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%TsPar0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' TsParT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%rGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' rGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%phiGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' phiGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%thetaGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' thetaGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%ellGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ellGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGL0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGLT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%qGH0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' qGHT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%pGC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' pGCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(:)) /= 1) .or. &
				(isnan(real(SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' d3xCT HAS', &
					' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%QCell0T(Qind)%d3xC0T(1) < 0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' NEGATIVE ', &
					' d3xCT FOR SPECIE= ', s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

  end subroutine ConfigGridGeneratorSub

end module ConfigGridGenerator
