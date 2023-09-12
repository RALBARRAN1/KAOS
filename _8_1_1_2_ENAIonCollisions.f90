module ENAIonCollisions

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.1.1.2 ENA-ION COLLISIONS:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! UPDATE KINETIC VARIABLE ENA FLAGS BY ENA-ION COLLISIONS:

	subroutine ENAIonCollisionsSub

		! ----------------------------------------------------

		! if (Qindk1(j) /= 0d0) then
!
! 			! ----------------------------------------------------
!
! 			! COMPUTE ENA-ION COLLISION FREQUENCIES:
! 			! Note: From ENA to ion use zeroth ion moment (plasma density).
!
! 			! ----------------------------------------------------
!
! 			sigmaENAIon(1)= pi*(2d0*SpecieT(s)%radsT(1))**2d0
!
! 				M1pT(1)= SpecieT(s)%FluxTubeT(f)%M1PerpphT(nn, Qind)/2d0
! 				M1qT(1)= SpecieT(s)%FluxTubeT(f)%M1ParphT(nn, Qind)
! 				M1phiT(1)= SpecieT(s)%FluxTubeT(f)%M1PerpphT(nn, Qind)/2d0
!
! 			if (rank == 0) then
! 				M1pT(1)= SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)/2d0
! 				M1qT(1)= SpecieT(s)%FluxTubeT(f)%M1ParphRT(nn, Qind)
! 				M1phiT(1)= SpecieT(s)%FluxTubeT(f)%M1PerpphRT(nn, Qind)/2d0
! 			end if
!
! 			M1xT(1)= (cos(phik1(j))/sqrt(ellk1(j)))*(3d0*M1qT(1)* &
! 				cos(thetak1(j))*sin(thetak1(j))+ M1pT(1)*(1d0- 3d0*(cos(thetak1(j)))**2d0))- &
! 				M1phiT(1)*sin(phik1(j))
!
! 			M1yT(1)= (sin(phik1(j))/sqrt(ellk1(j)))*(3d0*M1qT(1)* &
! 				cos(thetak1(j))*sin(thetak1(j))+ M1pT(1)*(1d0- 3d0*(cos(thetak1(j)))**2d0))+ &
! 				M1phiT(1)*cos(phik1(j))
!
! 			M1zT(1)= (M1qT(1)*(3d0*(cos(thetak1(j)))**2d0- 1d0)/sqrt(ellk1(j)))+ &
! 				(3d0*M1pT(1)*cos(thetak1(j))*sin(thetak1(j))/sqrt(ellk1(j)))
!
! 			VelENAIon(j)= abs(sqrt(abs(Vx(j)- M1xT(1))**2d0+ abs(Vy(j)- M1yT(1))**2d0+ &
! 				abs(Vz(j)- M1zT(1))**2d0))
!
! 			do nn= 1, SpecieT(s)%FluxTubeT(f)%NNtT(1), 1
! 				if ((n > (nn- 1)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1)) .and. &
! 					(n <= (nn)*SpecieT(s)%FluxTubeT(f)%ndatfacT(1))) then
! 					if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qindk1(j)) /= 0) then
!
! 						nuENAIon(j)= SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qindk1(j))* &
! 							SpecieT(s)%FluxTubeT(f)%NqReNormT(nn, Qindk1(j))* &
! 							sigmaENAIon(1)*VelENAIon(j)
!
! 					else if (SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qindk1(j)) == 0) then
!
! 						nuENAIon(j)= SpecieT(s)%FluxTubeT(f)%NqReNormENAT(nn, Qindk1(j)- 1)* &
! 							SpecieT(s)%FluxTubeT(f)%QCellT(Qindk1(j))%nsCNeutT(1)* &
! 							sigmaIonNeut(1)*VelIonNeut(j)
!
! 					end if
! 				end if
! 			end do
!
! 			tauENAIon(j)= 1d0/nuENAIon(j)
!
! 			! ----------------------------------------------------
!
! 			! DIAGNOSTIC FLAGS FOR PROPER ARRAY SIZES AND FINITE VALUES:
!
! 			if ((isnan(real(sigmaENAIon(1))) .eqv. .true.) .or. &
! 				(size(sigmaENAIon(:)) /= 1d0)) then
! 				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' sigmaENAIon= ', &
! 					sigmaENAIon(1), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
! 					s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
! 					' IN ENA-ION COLLISIONS SUBROUTINE' // achar(27) // '[0m.'
! 			end if
!
! 			if ((isnan(real(VelENAIon(j))) .eqv. .true.) .or. &
! 				(size(VelENAIon(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
! 				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' VelIonNeut= ', &
! 					VelENAIon(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
! 					s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
! 					' IN ENA-ION COLLISIONS SUBROUTINE' // achar(27) // '[0m.'
! 			end if
!
! 			if ((isnan(real(nuENAIon(j))) .eqv. .true.) .or. &
! 				(size(nuENAIon(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
! 				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' nuIonNeut= ', &
! 					nuENAIon(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
! 					s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
! 					' IN ENA-ION COLLISIONS SUBROUTINE' // achar(27) // '[0m.'
! 			end if
!
! 			if ((isnan(real(tauENAIon(j))) .eqv. .true.) .or. &
! 				(size(tauENAIon(:)) /= SpecieT(s)%FluxTubeT(f)%NsT(1))) then
! 				write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' tauENAIon= ', &
! 					tauENAIon(j), ' HAS BAD SIZE OR HAS NaN VALUE FOR SPECIE= ', &
! 					s, ', FLUX TUBE= ', f, ', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
! 					' IN ENA-ION COLLISIONS SUBROUTINE' // achar(27) // '[0m.'
! 			end if
!
! 			! ----------------------------------------------------
!
! 			! ENSURE COMPUTATIONAL TIME-STEP RESOLVES ENA-ION COLLISION FREQUENCIES:
!
! 			! if (SpecieT(s)%FluxTubeT(f)%hT(1) >= tauENAIon(j)) then
! ! 						write(*, *) achar(27) // '[33m ERROR: RANK= ', rank, ' COMPUTATIONAL ', &
! ! 							' TIME-STEP IS GREATER OR EQUAL TO tauIonNeut= ', &
! ! 							tauENAIon(j), ' FOR SPECIE= ', s, ', FLUX TUBE= ', f, &
! ! 							', TIME-STEP= ', n, ', AND PARTICLE= ', j, &
! ! 							' IN ENA-ION COLLISIONS SUBROUTINE' // achar(27) // '[0m.'
! ! 					end if
!
! 			! ----------------------------------------------------
!
! 			! CHANGE ENAflag VALUES FOR ENAs TO IONs ON ENA-ION COLLISION FREQUENCIES:
!
! 			! ----------------------------------------------------
!
! 			if (abs(n- ENAflagN1ind(j))*SpecieT(s)%FluxTubeT(f)%hT(1) &
! 				>= tauENAIon(j)) then
!
! 				ENAflag(j)= 0
! 				ENAflagN0ind(j)= n
!
! 				! Fold ENA translational velocities to ion perp velocity components
! 	1			! and 3D Cartesian components.
!
! 				Vperp1(j)= abs((1d0- 3d0*(cos(thetak1(j))**2d0))*(Vx(j)*cos(phik1(j))+ &
! 					Vy(j)*sin(phik1(j)))/sqrt(ellk1(j))+ 3d0*Vz(j)*cos(thetak1(j))* &
! 					sin(thetak1(j))/sqrt(ellk1(j)))
! 				Vperp2(j)= abs(Vy(j)*cos(phik1(j))- Vx(j)*sin(phik1(j)))
!
! 				Vparp(j)= 3d0*cos(thetak1(j))*sin(thetak1(j))*(Vx(j)*cos(phik1(j))+ &
! 					Vy(j)*sin(phik1(j)))/sqrt(ellk1(j))+ Vz(j)*(3d0*(cos(thetak1(j))**2d0)- 1d0)/ &
! 					sqrt(ellk1(j))
!
! 				! Note: Translational velocities are from parallel velocity component only.
! 				Vx(j)= 3d0*Vparp(j)*cos(thetak1(j))*sin(thetak1(j))*cos(phik1(j))/ &
! 					(sqrt(ellk1(j)))
!
! 				Vy(j)= 3d0*Vparp(j)*cos(thetak1(j))*sin(thetak1(j))*sin(phik1(j))/ &
! 					(sqrt(ellk1(j)))
!
! 				Vz(j)= Vparp(j)*(3d0*(cos(thetak1(j))**2d0)- 1d0)/ &
! 					(sqrt(ellk1(j)))
!
! 				if (abs(Vx(j)) < 1d-6) then ! Set Cartesian coords according to B longitude
! 					Vx(j)= 0d0
! 				end if
!
! 				if (abs(Vy(j)) < 1d-6) then
! 					Vy(j)= 0d0
! 				end if
!
! 				if (abs(Vz(j)) < 1d-6) then
! 					Vz(j)= 0d0
! 				end if
!
! 				write(*, *) 'Secondary Ion produced'
! 				write(*, *) '(n-ENAflagN1ind)*h test= ', &
! 					abs(n- ENAflagN1ind(j))*SpecieT(s)%FluxTubeT(f)%hT(1)
! 				write(*, *) 'tauENAIon test= ', &
! 					tauENAIon(j)
!
! 			end if
!
! 			! ----------------------------------------------------
!
! 		else if (Qindk1(j) == 0d0) then
! 			tauENAIon(j)= 0d0
! 		end if

		! ----------------------------------------------------

	end subroutine ENAIonCollisionsSub

end module ENAIonCollisions
