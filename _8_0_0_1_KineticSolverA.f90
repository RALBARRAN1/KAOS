module KineticSolverA

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.0.1 3D KINETIC SOLVER A:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams
use KineticRK4Update
use ParticleCounts
use DistributionFunctions
use ZerothIonMoment
use FirstPerpIonMoment
use FirstParIonMoment
use SecondIonMoment
use SecondPerpIonMoment
use SecondParIonMoment
use AmbipolarEfield
use Gravfield
use PotentialStructure
use DataExport1
use DataExport2

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

	subroutine KineticSolverASub

  	! ----------------------------------------------------

  	! SET ALL KINETIC VARIABLES AT INITIAL TIME:

  	if (n == 1) then

			nnind= 1

  		! ----------------------------------------------------

  		do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

  			! Set all initial particles as ions
  			ENAflag(j)= .false.
  			ENAflagN0ind(j)= n

  			if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) <= 0) then ! N. Magnetic Hemisphere

					! ----------------------------------------------------

					if (SpecieT(s)%FluxTubeT(f)%q0T(j) <= SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1)) then
						SpecieT(s)%FluxTubeT(f)%q0T(j)= SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1)
					else if (SpecieT(s)%FluxTubeT(f)%q0T(j) > SpecieT(s)%FluxTubeT(f)%qGHT(nn, (NqUB(1)- NqLB(1))+ 1)) then
						SpecieT(s)%FluxTubeT(f)%q0T(j)= SpecieT(s)%FluxTubeT(f)%qGHT(nn, ((NqUB(1)- NqLB(1))+ 1))
					end if

  				! ----------------------------------------------------

					QloopKSA1: do Qind= NqLB(1), NqUB(1), 1
						if ((Qind == NqLB(1)) .and. &
							(SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) <= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then
							Qindk0(j)= Qind

							exit QloopKSA1

						else if ((Qind /= NqLB(1)) .and. &
							(SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) < SpecieT(s)%FluxTubeT(f)%q0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%q0T(j) <= SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then
							Qindk0(j)= Qind

							exit QloopKSA1

						end if
					end do QloopKSA1

  				! ----------------------------------------------------

  				! DIAGNOSTIC FLAGS FOR ALL q0T VALUES WITHIN CONFIGURATION-SPACE GRID:

					if (((Qindk0(j) == 1) .and. &
						(SpecieT(s)%FluxTubeT(f)%q0T(j) < SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)))) &
						.or. ((Qindk0(j) == 1) .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) > &
						SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j))))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' q0T= ', &
							SpecieT(s)%FluxTubeT(f)%q0T(j), ' VALUE OUT OF Q GRID WHERE qGLT= ', &
							SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)), ' AND qGHT= ', &
							SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j)), 'FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', Qindk0= ', Qindk0(j), ', STATISTICAL TIME-STEP= ', nn, &
							', AND PARTICLE= ', j, ' IN KINETIC SOLVER A SUBROUTINE' &
							// achar(27) // '[0m.'
					end if
					if (((Qindk0(j) /= 1) .and. &
						(SpecieT(s)%FluxTubeT(f)%q0T(j) <= SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)))) &
						.or. ((Qindk0(j) /= 1) .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) > &
						SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j))))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' q0T= ', &
							SpecieT(s)%FluxTubeT(f)%q0T(j), ' VALUE OUT OF Q GRID WHERE qGLT= ', &
							SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)), ' AND qGHT= ', &
							SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j)), 'FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', Qindk0= ', Qindk0(j), ', STATISTICAL TIME-STEP= ', nn, &
							', AND PARTICLE= ', j, ' IN KINETIC SOLVER A SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

  				! ----------------------------------------------------

  			end if
  			if (SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1) > 0) then ! S. Magnetic Hemisphere

  				! ----------------------------------------------------

					if (SpecieT(s)%FluxTubeT(f)%q0T(j) > SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1)) then
						SpecieT(s)%FluxTubeT(f)%q0T(j)= SpecieT(s)%FluxTubeT(f)%qGLT(nn, 1)
					else if (SpecieT(s)%FluxTubeT(f)%q0T(j) < SpecieT(s)%FluxTubeT(f)%qGHT(nn, (NqUB(1)- NqLB(1))+ 1)) then
						SpecieT(s)%FluxTubeT(f)%q0T(j)= SpecieT(s)%FluxTubeT(f)%qGHT(nn, ((NqUB(1)- NqLB(1))+ 1))
					end if

  				! ----------------------------------------------------

					QloopKSA2: do Qind= NqLB(1), NqUB(1), 1
						if ((Qind == NqLB(1)) .and. &
							(SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) >= SpecieT(s)%FluxTubeT(f)%q0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then
							Qindk0(j)= Qind

							exit QloopKSA2

						else if ((Qind /= NqLB(1)) .and. &
							(SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qind) > SpecieT(s)%FluxTubeT(f)%q0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%q0T(j) >= SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qind))) then
							Qindk0(j)= Qind

							exit QloopKSA2

						end if
					end do QloopKSA2

  				! ----------------------------------------------------

  				! DIAGNOSTIC FLAGS FOR ALL q0T VALUES WITHIN CONFIGURATION-SPACE GRID:

					if (((Qindk0(j) == 1) .and. &
						(SpecieT(s)%FluxTubeT(f)%q0T(j) > SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)))) &
						.or. ((Qindk0(j) == 1) .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) < &
						SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j))))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' q0T= ', &
							SpecieT(s)%FluxTubeT(f)%q0T(j), ' VALUE OUT OF Q GRID WHERE qGLT= ', &
							SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)), ' AND qGHT= ', &
							SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j)), 'FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', Qindk0= ', Qindk0(j), ', STATISTICAL TIME-STEP= ', nn, &
							', AND PARTICLE= ', j, ' IN KINETIC SOLVER A SUBROUTINE' &
							// achar(27) // '[0m.'
					end if
					if (((Qindk0(j) /= 1) .and. &
						(SpecieT(s)%FluxTubeT(f)%q0T(j) >= SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)))) &
						.or. ((Qindk0(j) /= 1) .and. (SpecieT(s)%FluxTubeT(f)%q0T(j) < &
						SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j))))) then
						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' q0T= ', &
							SpecieT(s)%FluxTubeT(f)%q0T(j), ' VALUE OUT OF Q GRID WHERE qGLT= ', &
							SpecieT(s)%FluxTubeT(f)%qGLT(nn, Qindk0(j)), ' AND qGHT= ', &
							SpecieT(s)%FluxTubeT(f)%qGHT(nn, Qindk0(j)), ' FOR SPECIE= ', s, &
							', FLUX TUBE= ', f, ', Qindk0= ', Qindk0(j), ', STATISTICAL TIME-STEP= ', nn, &
							', AND PARTICLE= ', j, ' IN KINETIC SOLVER A SUBROUTINE' &
							// achar(27) // '[0m.'
					end if

  				! ----------------------------------------------------

  			end if

  			! ----------------------------------------------------

				if (SpecieT(1)%FluxTubeT(1)%ION2VPERPflagT(1) == 1) then

	  			VparloopKSA1: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
						if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(1, 1, 1)%VparGLT(1) > SpecieT(s)%FluxTubeT(f)%Vpar0T(j)) &
							.or. (SpecieT(s)%FluxTubeT(f)%Vpar0T(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(1, 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1))%VparGHT(1))) then
							Vparindk0(j)= 0d0

							exit VparloopKSA1

						end if
						if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(1, 1, Vparind)%VparGLT(1) <= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vpar0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
	  					Vparindk0(j)= Vparind

							exit VparloopKSA1

	  				end if
	  				if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(1, 1, Vparind)%VparGLT(1) < SpecieT(s)%FluxTubeT(f)%Vpar0T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vpar0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(1, 1, Vparind)%VparGHT(1))) then
	  					Vparindk0(j)= Vparind

							exit VparloopKSA1

	  				end if
	  			end do VparloopKSA1

	  			Vperp1loopKSA: do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
						if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(1, 1, 1)%Vperp1GLT(1) > SpecieT(s)%FluxTubeT(f)%Vperp10T(j)) &
							.or. (SpecieT(s)%FluxTubeT(f)%Vperp10T(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1, 1)%Vperp1GHT(1))) then
							Vperp1indk0(j)= 0d0

							exit Vperp1loopKSA

						end if
	  				if ((Vperp1ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) <= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vperp10T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
	  					Vperp1indk0(j)= Vperp1ind

							exit Vperp1loopKSA

	  				end if
	  				if ((Vperp1ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GLT(1) < SpecieT(s)%FluxTubeT(f)%Vperp10T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vperp10T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(Vperp1ind, 1, 1)%Vperp1GHT(1))) then
	  					Vperp1indk0(j)= Vperp1ind

							exit Vperp1loopKSA

	  				end if
	  			end do Vperp1loopKSA

					Vperp2loopKSA: do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
						if ((SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(1, 1, 1)%Vperp2GLT(1) > SpecieT(s)%FluxTubeT(f)%Vperp20T(j)) &
							.or. (SpecieT(s)%FluxTubeT(f)%Vperp20T(j) > SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							V2PerpCellT(1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1)%Vperp2GHT(1))) then
							Vperp2indk0(j)= 0d0

							exit Vperp2loopKSA

						end if
						if ((Vperp2ind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) <= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vperp20T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
	  					Vperp2indk0(j)= Vperp2ind

							exit Vperp2loopKSA

	  				end if
	  				if ((Vperp2ind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
	  					V2PerpCellT(1, Vperp2ind, 1)%Vperp2GLT(1) < SpecieT(s)%FluxTubeT(f)%Vperp20T(j)) &
	  					.and. (SpecieT(s)%FluxTubeT(f)%Vperp20T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
  						V2PerpCellT(1, Vperp2ind, 1)%Vperp2GHT(1))) then
	  					Vperp2indk0(j)= Vperp2ind

							exit Vperp2loopKSA

	  				end if
	  			end do Vperp2loopKSA

				else

					VparloopKSA2: do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1
						if ((Vparind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							VCellT(1, Vparind)%VparGLT(1) <= SpecieT(s)%FluxTubeT(f)%Vpar0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%Vpar0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							VCellT(1, Vparind)%VparGHT(1))) then
							Vparindk0(j)= Vparind

							exit VparloopKSA2

						end if
						if ((Vparind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							VCellT(1, Vparind)%VparGLT(1) < SpecieT(s)%FluxTubeT(f)%Vpar0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%Vpar0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
							VCellT(1, Vparind)%VparGHT(1))) then
							Vparindk0(j)= Vparind

							exit VparloopKSA2

						end if
					end do VparloopKSA2

					VperploopKSA: do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
						if ((Vperpind == 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qindk0(j))% &
							VCellT(Vperpind, 1)%VperpGLT(1) <= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%Vperp0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qindk0(j))% &
							VCellT(Vperpind, 1)%VperpGHT(1))) then
							Vperpindk0(j)= Vperpind

							exit VperploopKSA

						end if
						if ((Vperpind /= 1) .and. (SpecieT(s)%FluxTubeT(f)%QCellT(Qindk0(j))% &
							VCellT(Vperpind, 1)%VperpGLT(1) < SpecieT(s)%FluxTubeT(f)%Vperp0T(j)) &
							.and. (SpecieT(s)%FluxTubeT(f)%Vperp0T(j) <= SpecieT(s)%FluxTubeT(f)%QCellT(Qindk0(j))% &
							VCellT(Vperpind, 1)%VperpGHT(1))) then
							Vperpindk0(j)= Vperpind

							exit VperploopKSA

						end if
					end do VperploopKSA

				end if

  			! ----------------------------------------------------

  		end do

			! ----------------------------------------------------

			! GATHER STATISTICS AND FORM MACROSCOPIC PARAMETERS:

  		call ParticleCountsSub
  		call DistributionFunctionsSub
  		call ZerothIonMomentSub
  		call FirstPerpIonMomentSub
  		call FirstParIonMomentSub
  		call SecondIonMomentSub
  		call SecondPerpIonMomentSub
  		call SecondParIonMomentSub

			! ----------------------------------------------------

			! DIAGNOSTIC FLAG FOR NON-ZERO GRID CELLS IN SPIN-UP:

			if ((SpecieT(s)%FluxTubeT(f)%SPINUPflagT(1) == 1) .and. (rank == 0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
					if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
						(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then
						if (nn == SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1) then
							do Qind= NqLB(1), NqUB(1), 1
								if (SpecieT(s)%FluxTubeT(f)%M0phRT(nn, Qind) == 0) then
									write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, &
										' SPIN-UP SIMULATION HAS ZERO DENSITY AT FINAL TIME FOR SPECIE= ', &
										s, ', FLUX TUBE= ', f, ', AND Qind= ', Qind, ' IN KINETIC SOLVER A SUBROUTINE' &
										// achar(27) // '[0m.'
								end if
							end do
						end if
					end if
				end do
			end if

			! ----------------------------------------------------

			! Compute ambipolar and parallel electric fields
			call AmbipolarEfieldSub
			call GravfieldSub
			call PotentialStructureSub

  		! ----------------------------------------------------

  		jloopKS0: do j= 1, SpecieT(s)%FluxTubeT(f)%NsT(1), 1

  			! ----------------------------------------------------

				! Reset Ambipolar Electric Field Accelerations
				if (ENAflag(j) .eqv. .false.) then
					AEAmag(j)= AEAmagN(j)
					AGmag(j)= AGmagN(j)
				end if

				! Reset Parallel Electric Field Accelerations
				if (ENAflag(j) .eqv. .false.) then
					AEPmag(j)= AEPmagN(j)
				end if

  			! ----------------------------------------------------

  			x(j)= SpecieT(s)%FluxTubeT(f)%x0T(j)
  			y(j)= SpecieT(s)%FluxTubeT(f)%y0T(j)
  			z(j)= SpecieT(s)%FluxTubeT(f)%z0T(j)
  			Vperp1(j)= SpecieT(s)%FluxTubeT(f)%Vperp10T(j)
  			Vperp2(j)= SpecieT(s)%FluxTubeT(f)%Vperp20T(j)
  			Vperp(j)= SpecieT(s)%FluxTubeT(f)%Vperp0T(j)
  			Vx(j)= SpecieT(s)%FluxTubeT(f)%Vx0T(j)
  			Vy(j)= SpecieT(s)%FluxTubeT(f)%Vy0T(j)
  			Vz(j)= SpecieT(s)%FluxTubeT(f)%Vz0T(j)

  			! ----------------------------------------------------

  			call KineticRK4UpdateSub

				cycle jloopKS0

  			! ----------------------------------------------------

  		end do jloopKS0

  		! ----------------------------------------------------

  		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if (((n == 1) .and. (nn == 1)) .or. ((n /= 1) .and. (nn /= 1) .and. &
  				(n == (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)))) then

  				! ----------------------------------------------------

  				SpecieT(s)%FluxTubeT(f)%TimeT(nn)= Time(1)

  				! ----------------------------------------------------

  			end if
  		end do

  		! ----------------------------------------------------

  		! EXPORT ALL INITIAL KINETIC DATA:

  		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if (((n == 1) .and. (nn == 1)) .and. &
  				(n /= (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) .and. &
  				((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
  				(nn /= SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1))) then
  				call DataExport1Sub
  				call DataExport2Sub
  			end if
  		end do

			! Export injection data on injection time-step
			if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1)  .and. &
				(rank == 0)) then
				do nn= 1, SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1, 1
	  			if (((n == 1) .and. (nn == 1)) .and. &
	  				(n /= (nn- 1)*SpecieT(s)%FluxTubeT(f)%Q0ndatfacT(1)) .and. &
	  				((n /= SpecieT(s)%FluxTubeT(f)%NtT(1)) .and. &
	  				(nn /= SpecieT(s)%FluxTubeT(f)%Q0NNtT(1)+ 1))) then

						write(nnstring, '(I5)') nn
						write(sstring, '(I5)') s
						write(fstring, '(I5)') f

						expstring= adjustl(adjustr(rankstring) &
							// '_' // adjustl(adjustr(nnstring) &
							// '_' // adjustl(adjustr(sstring) &
							// '_' // adjustl(adjustr(fstring) // '_'))))

						NsnTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnT(nn)
						close(expint)

						NsnRRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NsnRRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NsnRRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NsnRRT(nn)
						close(expint)

						NqLBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqLBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxIonRT(nn)
						close(expint)

						LBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(LBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxIonRT(nn)
						close(expint)

						NqUBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(NqUBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxIonRT(nn)
						close(expint)

						UBoutfluxIonRTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxIonRTfort.bin')))
						open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
							adjustl(adjustr(UBoutfluxIonRTfile))), &
							status= 'replace', form= 'unformatted', access= 'stream')
						write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxIonRT(nn)
						close(expint)

						if (SpecieT(s)%FluxTubeT(f)%QEXCHANGEflagT(1) == 1) then
							NqLBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqLBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqLBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqLBoutfluxENART(nn)
							close(expint)

							LBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBoutfluxENART(nn)
							close(expint)

							NqUBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('NqUBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(NqUBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%NqUBoutfluxENART(nn)
							close(expint)

							UBoutfluxENARTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBoutfluxENARTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBoutfluxENARTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBoutfluxENART(nn)
							close(expint)
						end if

						if (SpecieT(s)%FluxTubeT(f)%LBCONDITIONflagT(1) == 1) then
							LBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('LBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(LBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%LBNetDensityT(nn)
							close(expint)

							UBNetDensityTfile= adjustl(adjustr(expstring) // adjustl(adjustr('UBNetDensityTfort.bin')))
							open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
								adjustl(adjustr(UBNetDensityTfile))), &
								status= 'replace', form= 'unformatted', access= 'stream')
							write(expint) SpecieT(s)%FluxTubeT(f)%UBNetDensityT(nn)
							close(expint)

						end if

					end if
  			end do
  		end if

  		! ----------------------------------------------------

  		! PRINT OUT STATISTICAL TIME-STEP:

  		do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1)+ 1, 1
  			if ((n == 1) .and. (nn == 1)) then
  				if (rank == 0) then
  					call cpu_time(KS0End)
  					write(nnstring, '(I5)') nn
  					write(sstring, '(I5)') s
  					write(fstring, '(I5)') f
  					write(Timestring, '(I5)') nint(Time(1))
  					write(KS0string, '(i10)')  nint(KS0End)
  					write(*, *) trim('** COMPLETE: STATISTICAL TIME-STEP= ' &
  						// adjustl(nnstring)) // &
  						trim(', RANK= ' // adjustl(rankstring)) // &
  						trim(', PARTICLE SPECIE= ' // adjustl(sstring)) // &
  						trim(', FLUX-TUBE= ' // adjustl(fstring)) // &
  						trim(', SIM-TIME= ' // adjustl(Timestring)) // &
  						trim(' s., REAL-TIME= ' // adjustl(KS0string)) // &
  						trim(' s.')
  				end if
  			end if
  		end do

  		! ----------------------------------------------------

  	end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine KineticSolverASub

end module KineticSolverA
