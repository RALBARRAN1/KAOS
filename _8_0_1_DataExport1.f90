module DataExport1

! ----------------------------------------------------
! ----------------------------------------------------

! %%%%%%	8.0.1- DATA EXPORT 1:	%%%%%%

! ----------------------------------------------------
! ----------------------------------------------------

use KineticMainParams

! ----------------------------------------------------

implicit none

! ----------------------------------------------------

contains

! ----------------------------------------------------

! EXPORT DATA AS BINARY FILES:

	subroutine DataExport1Sub

		! ----------------------------------------------------
		! ----------------------------------------------------

		if ((SpecieT(s)%FluxTubeT(f)%FLUIDIONEXPORTflagT(1) == 1) .and. &
			(rank == 0)) then

			write(sstring, '(I5)') s
			write(fstring, '(I5)') f

			expstring= adjustl(adjustr(rankstring) // '_' // adjustl(adjustr(sstring) &
				// '_' // adjustl(adjustr(fstring) // '_')))

			nsnormfacTfile= adjustl(adjustr(expstring) // adjustl(adjustr('nsnormfacTfort.bin')))
			open(unit= expint, file= adjustl(adjustr(dataexportdir) // &
				adjustl(adjustr(nsnormfacTfile))), &
				status= 'replace', form= 'unformatted', access= 'stream')
			write(expint) SpecieT(s)%FluxTubeT(f)%nsnormfacT(1)
			close(expint)

		end if

		! ----------------------------------------------------
		! ----------------------------------------------------

	end subroutine DataExport1Sub

end module DataExport1
