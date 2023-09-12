module NonUniform3DIntegrator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.9 NON-UNIFORM 3D INTEGRATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE 3D NUMERICAL INTEGRATION:

	subroutine NonUniform3DIntegratorSub

		! ----------------------------------------------------

		do Vpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVpGT(1), 1
			do Vqind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVqGT(1), 1
				do Vphiind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVphiGT(1), 1

					Sum1ENA(Vpind, Vqind, Vphiind, nn)= &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(Vpind, Vqind, Vphiind)%dVpCT(1))* &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(Vpind, Vqind, Vphiind)%dVqCT(1))* &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
						V3CellT(Vpind, Vqind, Vphiind)%dVphiCT(1))* &
						gggENA(Vpind, Vqind, Vphiind, nn)

				end do
			end do
		end do

		IIENA(1)= sum(Sum1ENA(:, :, :, nn))
		MMENA(1)= IIENA(1)

		! ----------------------------------------------------

	end subroutine NonUniform3DIntegratorSub

end module NonUniform3DIntegrator
