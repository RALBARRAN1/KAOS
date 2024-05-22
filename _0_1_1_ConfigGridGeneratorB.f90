module ConfigGridGeneratorB

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	0.1.1 - CONFIGURATION-SPACE GRID GENERATOR B:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine ConfigGridGeneratorSubB

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR STATISTICAL TIME-STEP SIZE:

		if (ndatfacG(1) < 3d0) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' ndatfacG= ', ndatfacG(1), &
				' IS NOT COMPATIBLE WITH KINETIC SOLVER TIME FORMAT FOR SPECIE= ', s, &
				' , FLUX TUBE= ', f, ', AND STATISTICAL TIME-STEP= ', nn, &
				' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! COMPUTE ION DENSITY PROFILE:

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

			! ----------------------------------------------------

			! COMPUTE STEADY-STATE ION DENSITY PROFILE WITH NORMALIZATION:
			! NOTE: SET ninormfac SUCH THAT NO FA CELL IS EMPTY.

			ICbbp(Qind)= 2d0*GG*ME*((hqC0(Qind)*dqC0(Qind)*cos(thetaGC0(Qind)))/ &
				((rGC0(Qind)**2d0)*(sqrt(ellGC0(Qind)))))

		end do

		IC0bb(1)= sum(ICbbp(1:SpecieT(s)%Qindns0GT(1)))

		do Qind= SpecieT(s)%FluxTubeT(f)%NqLBT(1), SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2, 1

			!Altitude dependent gravitational acceleration, integration over ds
			ICbb(1)= sum(ICbbp(1:Qind))
      gC(1)= ICbb(1)- IC0bb(1)

			if (qGA(1) <= 0d0) then ! SMH
        argC(1)= ((SpecieT(s)%msT(1)*gC(1))/(kB*(TsPar0(Qind)+ Te0(Qind))))
			end if

			if (qGA(1) > 0d0) then ! NMH
        argC(1)= -((SpecieT(s)%msT(1)*gC(1))/(kB*(TsPar0(Qind)+ Te0(Qind))))
			end if

			nsC(1)= ns0(1)*exp(argC(1)) ! Number of ions per flux-tube grid cell [m^-3]

			! Normalized number of macroparticles per flux-tube grid cell [unitless]
			if (INITIALGRIDflag == 1) then
				if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) then
					nsnormC0(Qind)= (nsC(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))*(d3xC0(Qind))
				end if
				if (SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 0) then
					if (Qind == SpecieT(s)%FluxTubeT(f)%NqLBT(1)) then
						nsnormC0(Qind)= SpecieT(s)%FluxTubeT(f)%nsnormCLBInputT(1)
					end if
					if (Qind == SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2) then
						nsnormC0(Qind)= SpecieT(s)%FluxTubeT(f)%nsnormCUBInputT(1)
					end if
					if ((Qind /= SpecieT(s)%FluxTubeT(f)%NqLBT(1)) .and. (Qind /= SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 2)) then
						nsnormC0(Qind)= &
							(SpecieT(s)%FluxTubeT(f)%QCellT(Qind- 1)%DensityInputT(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
							(d3xC0(Qind))
					end if
				end if
			end if

			if (INITIALGRIDflag == 0) then
				nsnormC0(Qind)= (nsC(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
					(d3xC0(Qind))
			end if

			! ----------------------------------------------------

			! COMPUTE STEADY-STATE O NEUTRAL DENSITY PROFILE:

			if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

				gCNeut(1)= (GG*ME)/(kB*TNeut(1))

				if (qGA(1) <= 0d0) then ! SMH
					argCNeut(1)= mNeut*gCNeut(1)*((1d0/zns0Neut(1))- &
						(1d0/rGC0(Qind)))
				end if
				if (qGA(1) > 0d0) then ! NMH
					argCNeut(1)= -mNeut*gCNeut(1)*((1d0/zns0Neut(1))- &
						(1d0/rGC0(Qind)))
				end if

				nsCNeut(1)= ns0Neut(1)*exp(argCNeut(1)) ! Number of H neutrals per flux-tube grid cell [m^-3]

				! Normalized number of macroparticles per flux-tube grid cell [unitless]
				nsnormCNeut(1)= (nsCNeut(1)/SpecieT(s)%FluxTubeT(f)%nsnormfacT(1))* &
					(d3xC0(Qind))

				nsnormCNeut0(Qind)= nsnormCNeut(1)

			end if

			! ----------------------------------------------------

			! DIAGNOSTIC FLAGS FOR PROPER ARRAY INVERSIONS, SIZES, AND FINITE VALUES:

			if ((size(gC(:)) /= 1) .or. (isnan(real(gC(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gC HAS', &
					' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(argC(:)) /= 1) .or. (isnan(real(argC(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argC HAS', &
					' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
			end if

			if ((size(nsC(:)) /= 1) .or. (isnan(real(nsC(1))) .eqv. .true.)) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsC HAS', &
					' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
					', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
					' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
			end if

			if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then

				if ((size(gCNeut(:)) /= 1) .or. (isnan(real(gCNeut(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' gCNeut HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(argCNeut(:)) /= 1) .or. (isnan(real(argCNeut(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' argCNeut HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((size(nsnormCNeut(:)) /= 1) .or. (isnan(real(nsnormCNeut(1))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeut HAS', &
						' BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
				end if

				if ((nsnormCNeut(1) /= nsnormCNeut0(Qind)) .or. &
					(size(nsnormCNeut0(:)) /= SpecieT(s)%FluxTubeT(f)%NqG0T(1)) .or. &
					(isnan(real(nsnormCNeut0(Qind))) .eqv. .true.)) then
					write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nsnormCNeutG HAS', &
						' BAD INVERSION, SIZE, OR HAS NaN VALUE FOR SPECIE= ', s, &
						', FLUX TUBE= ', f, ', AND Qind= ', Qind, &
						' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
				end if

			end if

			! ----------------------------------------------------

		end do

		! ----------------------------------------------------

		! SET LB AND UB INJECTION DENSITIES:

		nsnormCLBG(1)= nsnormC0(SpecieT(s)%FluxTubeT(f)%NqLBT(1)+ 1)
		nsnormCUBG(1)= nsnormC0(SpecieT(s)%FluxTubeT(f)%NqUBT(1)+ 1)

		! ----------------------------------------------------

		! DIAGNOSTIC FLAG FOR CONSISTENT LOWER GHOST CELL DENSITIY:

		if (SpecieT(s)%Qindns0GT(1) /= 1d0) then
			write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
				' LOWER GHOST CELL DOES NOT MATCH ALTITUDE OF REFERENCE ', &
				' DENSITY FOR SPECIE= ', s, ', AND FLUX TUBE= ', f, &
				' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
		end if

		! ----------------------------------------------------

		! DIAGNOSTIC FLAGS FOR PROPER GHOST CELL DENSITIES:

		if (SpecieT(1)%FluxTubeT(1)%LBCONDITIONflagT(1) == 1) then
			if (nsnormCLBG(1) == 0d0) then
				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
					' ZERO LOWER BOUNDARY DENSITY= ', nsnormCLBG(1), &
					' FOR SPECIE= ', s, ', AND STATISTICAL TIME-STEP= ', nn, ', FLUX TUBE= ', f, &
					' IN CONFIGURATION-SPACE GRID GENERATOR B SUBROUTINE' // achar(27) // '[0m.'
			end if
		end if

		! ----------------------------------------------------

  end subroutine ConfigGridGeneratorSubB

end module ConfigGridGeneratorB
