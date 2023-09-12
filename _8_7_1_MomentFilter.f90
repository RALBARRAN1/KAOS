module MomentFilter

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.7.1 MomentFilter:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! FILTER MOMENTS BY SLIDING AVERAGE:

	subroutine MomentFilterSub

		! ----------------------------------------------------
		! ----------------------------------------------------

	  if ((Qind >= NqLB(1)) .and. (Qind <= NqUB(1))) then

			! ----------------------------------------------------

			! Filter low altitude moments
			if (Qind <= nint(NqLB(1)- 1d0+ (MAfilterPt(1)- 1d0)/2d0)) then

				! Average moments above (and including) given cell
				FiltSum1(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)+ &
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind+ 1)
				do MAfilterQind1= 2, nint((MAfilterPt(1)- 1d0)/2d0), 1
					FiltSum1(1)= FiltSum1(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind+ MAfilterQind1))
				end do

				FiltAvrg1(1)= FiltSum1(1)/(nint((MAfilterPt(1)- 1d0)/2d0)+ 1)

				! Average moments below given cell
				if (Qind == NqLB(1)) then
					FiltSum2(1)= 0d0

					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= FiltAvrg1(1)

				else if (Qind == NqLB(1)+ 1d0) then
					FiltSum2(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, NqLB(1))

					FiltAvrg2(1)= FiltSum2(1)
					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= (FiltAvrg1(1)+ FiltAvrg2(1))/2d0

				else
					FiltSum2(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind- 1)
					do MAfilterQind2= 2, nint(real(Qind- NqLB(1))), 1
						FiltSum2(1)= FiltSum2(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind- MAfilterQind2))
					end do

					FiltAvrg2(1)= FiltSum2(1)/nint(real(Qind- NqLB(1)))
					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= (FiltAvrg1(1)+ FiltAvrg2(1))/2d0

				end if

			end if

			! ----------------------------------------------------

			! Filter mid altitude moments
			if ((Qind > nint(NqLB(1)- 1d0+ (MAfilterPt(1)- 1d0)/2d0)) .and. &
				(Qind < (NqUB(1)- nint((MAfilterPt(1)- 1d0)/2d0)))) then

				! Average moments above (and including) given cell
				FiltSum1(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)+ &
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind+ 1)
				do MAfilterQind1= 2, nint((MAfilterPt(1)- 1d0)/2d0), 1
					FiltSum1(1)= FiltSum1(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind+ MAfilterQind1))
				end do

				! Average moments below given cell
				FiltSum2(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind- 1)
				do MAfilterQind2= 2, nint((MAfilterPt(1)- 1d0)/2d0), 1
					FiltSum2(1)= FiltSum2(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind- MAfilterQind2))
				end do

				FiltAvrg1(1)= FiltSum1(1)/(nint((MAfilterPt(1)- 1d0)/2d0)+ 1)
				FiltAvrg2(1)= FiltSum2(1)/nint((MAfilterPt(1)- 1d0)/2d0)

				SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= (FiltAvrg1(1)+ FiltAvrg2(1))/2d0

			end if

			! ----------------------------------------------------

			! Filter high altitude moments
			if (Qind >= (NqUB(1)- nint((MAfilterPt(1)- 1d0)/2d0))) then

				! Average moments below (and including) given cell
				FiltSum2(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind)+ &
					SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind- 1)
				do MAfilterQind2= 2, nint((MAfilterPt(1)- 1d0)/2d0), 1
					FiltSum2(1)= FiltSum2(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind- MAfilterQind2))
				end do

				FiltAvrg2(1)= FiltSum2(1)/(nint((MAfilterPt(1)- 1d0)/2d0)+ 1)

				! Average moments above given cell
				if (Qind == NqUB(1)) then
					FiltSum1(1)= 0d0

					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= FiltAvrg2(1)

				else if (Qind == NqUB(1)- 1d0) then
					FiltSum1(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, NqUB(1))

					FiltAvrg1(1)= FiltSum1(1)
					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= (FiltAvrg1(1)+ FiltAvrg2(1))/2d0

				else
					FiltSum1(1)= SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, Qind+ 1)
					do MAfilterQind1= 2, nint(real(NqUB(1)- Qind)), 1
						FiltSum1(1)= FiltSum1(1)+ SpecieT(s)%FluxTubeT(f)%MomentFiltInT(nn, (Qind+ MAfilterQind1))
					end do

					FiltAvrg1(1)= FiltSum1(1)/nint(real(NqUB(1)- Qind))
					SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= (FiltAvrg1(1)+ FiltAvrg2(1))/2d0

				end if

			end if

			! ----------------------------------------------------

    else

      SpecieT(s)%FluxTubeT(f)%MomentFiltOutT(nn, Qind)= 0d0

    end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine MomentFilterSub

end module MomentFilter
