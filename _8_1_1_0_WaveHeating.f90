module WaveHeating

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.1.0 WAVE HEATING:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use GaussianRNG

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine WaveHeatingSub

		! ----------------------------------------------------

		! PERFORM ICR HEATING:

		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.) &
			.and. (Qindk1(j) /= 0) .and. (Qindk1(j) /= -1)) then

			! ----------------------------------------------------

			! SET POSITION-DEPENDENT ELF WAVE POWER DIRECTIONS AND SPECTRAL INDICES FOR ICR HEATING:

			lambdaPerp(1)= SpecieT(s)%FluxTubeT(f)%lambdaPerppT(nnind, Qindk1(j))
			EtaLH(1)= SpecieT(s)%FluxTubeT(f)%EtaLHpT(nnind, Qindk1(j))
			XiPerp1(1)= SpecieT(s)%FluxTubeT(f)%XiPerp1pT(nnind, Qindk1(j))
			XiPerp2(1)= SpecieT(s)%FluxTubeT(f)%XiPerp2pT(nnind, Qindk1(j))
			S0(1)= SpecieT(s)%FluxTubeT(f)%S0pT(nnind, Qindk1(j))
			OmegaG0(1)= SpecieT(s)%FluxTubeT(f)%OmegaG0pT(nnind, Qindk1(j))
			ChiPerp1(1)= SpecieT(s)%FluxTubeT(f)%ChiPerp1pT(nnind, Qindk1(j))
			ChiPerp2(1)= SpecieT(s)%FluxTubeT(f)%ChiPerp2pT(nnind, Qindk1(j))

		end if

		if ((SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 1) .and. (ENAflag(j) .eqv. .false.) &
			.and. (Qindk1(j) /= 0) .and. (Qindk1(j) /= -1)) then

			! ----------------------------------------------------

			! Compute Gaussian random variates
			UniformRN1(1)= U1GammaPerp1randn(1)
			UniformRN2(1)= U2GammaPerp1randn(1)
			call GaussianRNGSub
			GammaPerp1(1)= GaussianRN(1)

			UniformRN1(1)= U1GammaPerp2randn(1)
			UniformRN2(1)= U2GammaPerp2randn(1)
			call GaussianRNGSub
			GammaPerp2(1)= GaussianRN(1)

			! ----------------------------------------------------

			! COMPUTE ICR HEATING TRANSVERSE VELOCITY SPACE DIFFUSION COEFFICIENTS
			! FOR INFINITE OR FINITE PERPENDICULAR BBELF WAVELENGTHS:
			! Note: Refer to diffusion coefficients per (Barghouthi '06) for altitude dependence

			if ((lambdaPerp(1) == 0d0) .or. ((lambdaPerp(1) /= 0d0) .and. & ! Long wavelength limit
				(((2d0*pi*(sqrt(Vperp1(j)**2d0+ Vperp2(j)**2d0)))/(OmegaGk1(1))) < lambdaPerp(1)))) then
				SigmaPerp(1)= 1d0
			end if

			if ((lambdaPerp(1) /= 0d0) .and. & ! Short wavelength limit
				(((2d0*pi*(sqrt(Vperp1(j)**2d0+ Vperp2(j)**2d0)))/(OmegaGk1(1))) >= lambdaPerp(1))) then
				SigmaPerp(1)= (2d0*pi*(sqrt(Vperp1(j)**2d0+ Vperp2(j)**2d0)))/(lambdaPerp(1)*OmegaGk1(1))
			end if

			! Diffusion coefficients [m^2 s^-3]
			DPerp1(1)= ((SpecieT(s)%qsT(1)**2d0/(2d0*SpecieT(s)%msT(1)**2d0))* &
				(XiPerp1(1)*EtaLH(1))*(S0(1))* &
				(OmegaGk1(1)/OmegaG0(1))**(-ChiPerp1(1)))*(SigmaPerp(1)**(-3d0))
			DPerp2(1)= ((SpecieT(s)%qsT(1)**2d0/(2d0*SpecieT(s)%msT(1)**2d0))* &
				(XiPerp2(1)*EtaLH(1))*(S0(1))* &
				(OmegaGk1(1)/OmegaG0(1))**(-ChiPerp2(1)))*(SigmaPerp(1)**(-3d0))

			! ----------------------------------------------------

			! COMPUTE WAVE-PARTICLE INTERACTION TIME FROM (Schulz and Lanzerotti '74):

			! Compute field line arc length of heating region
			SpecieT(s)%FluxTubeT(f)%dsICRT(nnind, Qindk1(j))= &
				SpecieT(s)%FluxTubeT(f)%hqCT(nnind, Qindk1(j))*SpecieT(s)%FluxTubeT(f)%dqCT(nnind, Qindk1(j))
			SUMdsICR(1)= sum(SpecieT(s)%FluxTubeT(f)%dsICRT(nnind, :))

			! Compute epsilon factor (unitless)
			epsPerp(1)= abs(sqrt((VxR(1)**(2d0))+ (VyR(1)**(2d0))+ (VzR(1)**(2d0)))/ &
				(OmegaGk1(1)*SUMdsICR(1)))

			! ICR interaction time [s]
			tauPerp(1)= 1d0/(((OmegaGk1(1)/(2d0*pi)))*abs(sqrt(epsPerp(1))))

			! ----------------------------------------------------

			! ENSURE ICR HEATING EPSILON FACTOR IS LESS THAN UNITY:

			!if (epsPerp(1) >= 1d0) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ICR heating Epsilon factor= ', epsPerp(1), &
			!		' IS GREATER THAN UNITY FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
			!		', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
			!		// achar(27) // '[0m.'
			!end if

			!if (SpecieT(s)%FluxTubeT(f)%hT(1) > tauPerp(1)) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ICR interaction time= ', tauPerp(1), &
			!		' IS LESS THAN COMPUTATIONAL TIME-STEP FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
			!		', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
			!		// achar(27) // '[0m.'
			!end if

			! ----------------------------------------------------

			! COMPUTE TRANSVERSE VELOCITY KICKS:
			! Note: Heat ions by ICR only if simulation time-step is less than ICR interaction time
			! Scale Vperp speed kicks [m/s] to ICR interaction time

			! Resonant cyclotron heating
			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 1) then
				!if (SpecieT(s)%FluxTubeT(f)%hT(1) <= tauPerp(1)) then
	        DVperpicr(1)= (SpecieT(s)%FluxTubeT(f)%hT(1)/tauPerp(1))* &
						sqrt(2d0*(DPerp1(1)+ DPerp2(1))*tauPerp(1))*GammaPerp1(1)
				!end if
				!if (SpecieT(s)%FluxTubeT(f)%hT(1) > tauPerp(1)) then
				!	DVperpicr(1)= 0d0
				!end if
			end if

			! Stochastic cyclotron heating
			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 0) then
				if (SpecieT(s)%FluxTubeT(f)%hT(1) <= tauPerp(1)) then
					DVperp1icr(1)= (SpecieT(s)%FluxTubeT(f)%hT(1)/tauPerp(1))* &
						sqrt(2d0*DPerp1(1)*tauPerp(1))*GammaPerp1(1)
					DVperp2icr(1)= (SpecieT(s)%FluxTubeT(f)%hT(1)/tauPerp(1))* &
						sqrt(2d0*DPerp2(1)*tauPerp(1))*GammaPerp2(1)
				end if
				if (SpecieT(s)%FluxTubeT(f)%hT(1) > tauPerp(1)) then
					DVperp1icr(1)= 0d0
					DVperp2icr(1)= 0d0
				end if
			end if

			! ----------------------------------------------------

			! UPDATE TRANSVERSE VELOCITIES FOR COHERENT HEATING:

			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 1) then

				! Update Vperp values
				VperpN(j)= Vperp(j)+ DVperpicr(1)

				! Compute uniformly distributed random gyro-angle
				UniformGA(1)= GammaGArandn(1)
				call random_number(UniformGA) ! Uniformly distributed RN from [0, 1)
				URNgyroangle(1)= UniformGA(1)*2d0*pi ! Uniformly distributed gyro-angle from [0, 2*pi) [rads]

				! ----------------------------------------------------

				! ENSURE GYRO-ANGLE IS WITHIN BOUNDS:

				if ((URNgyroangle(1) < 0d0) .or. (URNgyroangle(1) >= 2d0*pi)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
						' COHERENT GYRO-ANGLE= ', URNgyroangle(1), &
						' IS OUT OF BOUNDS FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
						', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
						// achar(27) // '[0m.'
				end if

				! ----------------------------------------------------

				! RECONSTRUCT Vperp1 AND Vperp2 VALUES WITH UNIFORMLY DISTRIBUTED GYRO-ANGLE FOR COHERENT BBELF WAVES:

				if ((0d0 .lt. URNgyroangle(1)) .and. (URNgyroangle(1) .le. pi/2d0)) then
					Vperp1N(j)= abs(VperpN(j)*cos(URNgyroangle(1)))
					Vperp2N(j)= abs(VperpN(j)*sin(URNgyroangle(1)))

					! ----------------------------------------------------

					! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

					if (Vperp1N(j) < 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp1N= ', Vperp1N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (Vperp2N(j) < 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp2N= ', Vperp2N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if
				if ((pi/2d0 .lt. URNgyroangle(1)) .and. (URNgyroangle(1) .le. pi)) then
					Vperp1N(j)= -abs(VperpN(j)*cos(pi- URNgyroangle(1)))
					Vperp2N(j)= abs(VperpN(j)*sin(pi- URNgyroangle(1)))

					! ----------------------------------------------------

					! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

					if (Vperp1N(j) > 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp1N= ', Vperp1N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (Vperp2N(j) < 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp2N= ', Vperp2N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if
				if ((pi .lt. URNgyroangle(1)) .and. (URNgyroangle(1) .le. 3d0*pi/2d0)) then
					Vperp1N(j)= -abs(VperpN(j)*cos(URNgyroangle(1)- pi))
					Vperp2N(j)= -abs(VperpN(j)*sin(URNgyroangle(1)- pi))

					! ----------------------------------------------------

					! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

					if (Vperp1N(j) > 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp1N= ', Vperp1N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (Vperp2N(j) > 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp2N= ', Vperp2N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if
				if ((3d0*pi/2d0 .lt. URNgyroangle(1)) .and. (URNgyroangle(1) .le. 2d0*pi)) then
					Vperp1N(j)= abs(VperpN(j)*cos(2d0*pi- URNgyroangle(1)))
					Vperp2N(j)= -abs(VperpN(j)*sin(2d0*pi- URNgyroangle(1)))

					! ----------------------------------------------------

					! ENSURE CORRECT SIGNS OF Vperp1 AND Vperp2 VALUES:

					if (Vperp1N(j) < 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp1N= ', Vperp1N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					if (Vperp2N(j) > 0d0) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
							' INCORRECT COHERENT ICR SIGN OF Vperp2N= ', Vperp2N(j), &
							' FOR GYRO-ANGLE [rads]= ', URNgyroangle(1), &
							' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
							', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

					! ----------------------------------------------------

				end if

				! ----------------------------------------------------

			end if

			! ----------------------------------------------------

			! UPDATE TRANSVERSE VELOCITIES FOR STOCHASTIC HEATING:

			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 0) then

				Vperp1N(j)= Vperp1(j) + DVperp1icr(1)
				Vperp2N(j)= Vperp2(j) + DVperp2icr(1)
				VperpN(j)= abs(sqrt(Vperp1N(j)**2d0+ Vperp2N(j)**2d0))

			end if

			! ----------------------------------------------------

			! ENSURE LARMOR RADIUS IS SMALLER THAN BBELF WAVELENGTH (Barghouthi '94) (Bouhram '03a) (Zeng '06) (Wu '02)
			! AND FLUX-TUBE ARC LENGTH:

			!if (VperpN(j)/OmegaGk1(1) >= lambdaPerp(1)) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION LARMOR RADIUS [m]= ', VperpN(j)/OmegaGk1(1), &
			!		' IS GREATER THAN VLF WAVELENGTH [m]= ', lambdaPerp(1), ' AT LOCAL GYROFREQUENCY FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
			!		', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
			!		// achar(27) // '[0m.'
			!end if

			!if (VperpN(j)/OmegaGk1(1) >= &
			!	sum(SpecieT(s)%FluxTubeT(f)%hqCT(nn, :)*SpecieT(s)%FluxTubeT(f)%dqCT(nn, :))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ION LARMOR RADIUS [m]= ', VperpN(j)/OmegaGk1(1), &
			!		' IS GREATER THAN TOTAL FLUX-TUBE ARC LENGTH [m]= ', &
			!		sum(SpecieT(s)%FluxTubeT(f)%hqCT(nn, :)*SpecieT(s)%FluxTubeT(f)%dqCT(nn, :)), &
			!		' AT LOCAL GYROFREQUENCY FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
			!		', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
			!		// achar(27) // '[0m.'
			!end if

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR ALL Vperp VALUES WITHIN VELOCITY-SPACE GRID:

			!if ((VperpN(j) <= 0d0) .or. (VperpN(j) > SpecieT(s)%FluxTubeT(f)% &
			!	QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
			!	SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VperpGHT(1))) then
			!	write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
			!		' Vperp= ', sqrt(Vperp1(j)**2d0+ Vperp2(j)**2d0), &
			!		' VperpN= ', VperpN(j), ' ION VALUE OUT OF GRID WITH VperpGL= ', &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%VCellT(1, 1)%VperpGLT(1), &
			!		' AND VperpGH= ', SpecieT(s)%FluxTubeT(f)% &
			!		QCellT(1)%VCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), &
			!		SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VperpGHT(1), &
			!		'FOR SPECIE= ', s, &
			!		', FLUX TUBE= ', f, ', TIME-STEP= ', n, &
			!		', AND PARTICLE= ', j, ' IN WAVE HEATING SUBROUTINE' &
			!		// achar(27) // '[0m.'
				!call MPI_FINALIZE(ierr)
				!stop
			!end if

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

			if ((isnan(real(lambdaPerp(1))) .eqv. .true.) .or. &
				(size(lambdaPerp(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' lambdaPerp= ', lambdaPerp(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(EtaLH(1))) .eqv. .true.) .or. &
				(size(EtaLH(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' EtaLH= ', EtaLH(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(XiPerp1(1))) .eqv. .true.) .or. &
				(size(XiPerp1(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp1= ', XiPerp1(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(XiPerp2(1))) .eqv. .true.) .or. &
				(size(XiPerp2(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' XiPerp2= ', XiPerp2(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(S0(1))) .eqv. .true.) .or. &
				(size(S0(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' S0= ', S0(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(OmegaG0(1))) .eqv. .true.) .or. &
				(size(OmegaG0(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' OmegaG0= ', OmegaG0(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(ChiPerp1(1))) .eqv. .true.) .or. &
				(size(ChiPerp1(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp1= ', ChiPerp1(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(ChiPerp2(1))) .eqv. .true.) .or. &
				(size(ChiPerp2(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ChiPerp2= ', ChiPerp2(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(GammaPerp1(1))) .eqv. .true.) .or. &
				(size(GammaPerp1(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' GammaPerp1= ', GammaPerp1(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(GammaPerp2(1))) .eqv. .true.) .or. &
				(size(GammaPerp2(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' GammaPerp2= ', GammaPerp2(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(DPerp1(1))) .eqv. .true.) .or. &
				(size(DPerp1(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DPerp1= ', DPerp1(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(DPerp2(1))) .eqv. .true.) .or. &
				(size(DPerp2(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DPerp2= ', DPerp2(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(epsPerp(1))) .eqv. .true.) .or. &
				(size(epsPerp(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' epsPerp= ', epsPerp(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(tauPerp(1))) .eqv. .true.) .or. &
				(size(epsPerp(:)) /= 1)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' tauPerp= ', tauPerp(1), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 1) then
				if ((isnan(real(DVperpicr(1))) .eqv. .true.) .or. &
					(size(DVperpicr(:)) /= 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DVperpicr= ', DVperpicr(1), &
						' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
						', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if
			end if

			if (SpecieT(s)%FluxTubeT(f)%ICRCOHERENCEflagT(1) == 0) then
				if ((isnan(real(DVperp1icr(1))) .eqv. .true.) .or. &
					(size(DVperp1icr(:)) /= 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DVperp1icr= ', DVperp1icr(1), &
						' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
						', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((isnan(real(DVperp2icr(1))) .eqv. .true.) .or. &
					(size(DVperp2icr(:)) /= 1)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' DVperp2icr= ', DVperp2icr(1), &
						' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
						', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
						' SUBROUTINE' // achar(27) // '[0m.'
				end if
			end if

			if ((isnan(real(Vperp1N(j))) .eqv. .true.) .or. &
				(size(Vperp1N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerp1N= ', VPerp1N(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(Vperp2N(j))) .eqv. .true.) .or. &
				(size(Vperp2N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerp2N= ', VPerp2N(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(VperpN(j))) .eqv. .true.) .or. &
				(size(VperpN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerpN= ', VPerpN(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		else if ((((SpecieT(s)%FluxTubeT(f)%ICRflagT(1) == 0) .and. (ENAflag(j) .eqv. .false.)) .or. &
			((SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) .and. &
			(ENAflag(j) .eqv. .true.))) .and. (Qindk1(j) /= 0) .and. (Qindk1(j) /= -1)) then

			Vperp1N(j)= Vperp1(j) ! Set drifting MB distrib for no ICR heating
			Vperp2N(j)= Vperp2(j)
			VperpN(j)= Vperp(j)

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:

			if ((isnan(real(Vperp1N(j))) .eqv. .true.) .or. &
				(size(Vperp1N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerp1N= ', VPerp1N(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(Vperp2N(j))) .eqv. .true.) .or. &
				(size(Vperp2N(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerp2N= ', VPerp2N(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((isnan(real(VperpN(j))) .eqv. .true.) .or. &
				(size(VperpN(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VPerpN= ', VPerpN(j), &
					' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
					', TIME-STEP= ', n, ', AND PARTICLE= ', j, ' IN WAVE HEATING', &
					' SUBROUTINE' // achar(27) // '[0m.'
			end if

			! ----------------------------------------------------

		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine WaveHeatingSub

end module WaveHeating
