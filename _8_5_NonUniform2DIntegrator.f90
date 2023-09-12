module NonUniform2DIntegrator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.5 NON-UNIFORM 2D INTEGRATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE 2D NUMERICAL INTEGRATION:

	subroutine NonUniform2DIntegratorSub

		! ----------------------------------------------------

		do Vperpind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperpGT(1), 1
			do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

				Sum1(Vperpind, Vparind, nn)= &
					(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					VCellT(Vperpind, Vparind)%dVperpGT(1))* &
					(SpecieT(s)%FluxTubeT(f)%QCellT(1)% &
					VCellT(Vperpind, Vparind)%dVparGT(1))* &
					ggg(Vperpind, Vparind, nn)

			end do
		end do

		II(1)= sum(Sum1(:, :, nn))
		MM(1)= 2d0*pi*II(1)

		! ----------------------------------------------------

	end subroutine NonUniform2DIntegratorSub

end module NonUniform2DIntegrator
