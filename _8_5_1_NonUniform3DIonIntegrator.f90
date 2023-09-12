module NonUniform3DIonIntegrator

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.5.1 NON-UNIFORM 3D ION INTEGRATOR:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! COMPUTE 3D NUMERICAL INTEGRATION:

	subroutine NonUniform3DIonIntegratorSub

		! ----------------------------------------------------

    do Vperp1ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp1GT(1), 1
      do Vperp2ind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVperp2GT(1), 1
        do Vparind= 1, SpecieT(s)%FluxTubeT(f)%QCellT(1)%NVparGT(1), 1

					Sum12Perp(Vperp1ind, Vperp2ind, Vparind, nn)= &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp1GT(1))* &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVperp2GT(1))* &
						(SpecieT(s)%FluxTubeT(f)%QCellT(1)%V2PerpCellT(Vperp1ind, Vperp2ind, Vparind)%dVparGT(1))* &
						ggg2Perp(Vperp1ind, Vperp2ind, Vparind, nn)

				end do
			end do
		end do

		II(1)= sum(Sum12Perp(:, :, :, nn))
		MM(1)= II(1)

		! ----------------------------------------------------

	end subroutine NonUniform3DIonIntegratorSub

end module NonUniform3DIonIntegrator
