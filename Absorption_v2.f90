! ********************************************************************************************
! Program to calculate the integral and FM absorption from the transverse speed distributions
! ********************************************************************************************

	include "Maths/fmDetection_v2.f90"

	Program Absorption

	use fmDetection_v2

	Implicit None

! ****************************************
! Define variables
! ****************************************

	Integer:: i, j, nfiles, npad, nAbs, nFMPoints, numberOfSpeedPoints, npaddedAbs
	Double Precision :: lowlimit, uplimit, stepSize, modFreq, modStep, speedStep, min, maxAbs, maxFM, minX
	Double Precision, allocatable :: freq(:), amplitude(:)
	Double Precision, Allocatable :: transSpeed(:), FMabs(:), absLine(:), sp(:), spX(:), paddedAbs(:)
	Character(8), Parameter:: inpath = 'Outputs/'
	Character(8), Parameter:: outpath = 'Outputs/'
	Character(15):: absFile, FMFile
	Character(18):: infile
	Logical :: norm


! ****************************************
! Read input files
! ****************************************

	Open(unit=9,file='Inputs/fmInputs.inp')
	Read(9,*) lowlimit							! Minimum Speed
	Read(9,*) uplimit							! Maximum Speed
	Read(9,*) ! bin								!
	Read(9,*) stepSize							! Step size in MHz
	Read(9,*) ! kTransition				
	Read(9,*) modFreq							! Modulation frequency in MHz
	Close(9)

	Open(unit=10,file='Inputs/binParameters.inp')
	Read(10,*)
	Read(10,*)
	Read(10,*) norm								! Normalize output T/F
	Read(10,*) nfiles							! Number of files


	Write(*,*) "Generating Absorption and FM Lineshape"


	Do i = 1, nfiles							! For each input file

		write(infile,'("SpeedDistr_",I3.3,".txt")')i			! Generate input filename
		Open(unit=1000+i,file= inpath//infile)				! Open file
		Read(1000+i,*) speedStep, numberOfSpeedPoints			! Read speed step and number of speed points


!*****************************************
! Calculate parameters for array sizes
!*****************************************

		modStep = modFreq / stepSize					! Calculate ratio of modulation frequency to step size	
		npad = 3 * modStep						! Calculate number of padding points
		nAbs = 2 * (npad + numberOfSpeedPoints -1)			! Calculate number of points in integrated absorption lineshape
!		nFMPoints = nAbs - (2*npad)					! Calculate number of points in FM lineshape

                nPaddedAbs = 2 * (npad + numberOfSpeedPoints)
		nFMPoints = nPaddedAbs - (6 * modStep)
		
		If (Allocated(absLine)) then					! Deallocate all arrays from previous iteration (if one is allocated all were)
			Deallocate(absLine)					! Memory issues are created if these aren't deallocated at the start of the loop
			Deallocate(transSpeed)
			Deallocate(freq)
			Deallocate(amplitude)
			Deallocate(FMabs)
			Deallocate(Sp)
			Deallocate(SpX)
			Deallocate(paddedAbs)
		End if

	

!**************************************
! Allocate array sizes
!**************************************

		allocate(transSpeed(numberOfSpeedPoints))			! Allocate arrays for this iteration
		allocate(FMabs(nFMPoints))
		allocate(Sp(nFMPoints))
		allocate(absLine(nAbs))
		allocate(SpX(nAbs))
		allocate(freq(numberOfSpeedPoints))
		allocate(amplitude(numberOfSpeedPoints))
		allocate(paddedAbs(npaddedAbs))


		Do j = 1, numberOfSpeedPoints					! For each speed
			Read(1000+i,*) transSpeed(j), amplitude(j)		! Read in transverse speed / absorption line 

		End do
		Close(1000+i)

!******************************************
! Flip Speed Array
!******************************************
!	As we are reading in speeds (rather than velocities) the speed distribution is flipped about zero
!	to give the negative component we are assuning that the speed distribution is symmetric about zero
!	NOTE that because of the geometry of the G-L experiment the process of obtaining a FM lineshape from
!	speed distribution is different from the one used in photodissociation


	Write(absFile,'("Abs_",I3.3,".txt")')i					! Create output filename for 
	Open(unit=2000+i, file=outpath//absFile)				! Open file

	maxAbs=0.0

		Do j = 1, 2*numberOfSpeedPoints						! Create an array twice the length of the speed distribution
			If (j .le. numberOfSpeedPoints) then				! The first (negative) half is the reverse of the input speed distrib
				Sp(j) = -1*transSpeed(numberOfSpeedPoints-j+1)		! The speeds used are -1 time that in the positive half
				absLine(j) = amplitude(numberOfSpeedPoints-j+1) 	! The amplitudes are the same
			else 
				Sp(j) = transSpeed(j-numberOfSpeedPoints)		! Create the second (positive) half - these are the speeds 
				absLine(j) = amplitude(j-numberOfSpeedPoints)		! and amplitudes as read in
			End if

			If(absLine(j) .gt. maxAbs) then
				maxAbs = absLine(j)
			End if
		End do

		

		Do j = 1, 2*numberOfSpeedPoints
			If(norm .eqv. .true.) then
				absLine(j) = absLine(j)/maxAbs
			End if

			Write(2000+i,*) Sp(j), absLine(j)				! Write absorption lineshape to file (there isn't as much need for this
		End do									! as there would be for photodissociation data but it is useful to have
											! all the same)

!******************************************************************
! Create a padded absorption array prior to doing the FM conversion
!******************************************************************

		paddedAbs = 0.0								! Create an array of zeros
		Do j = 1, npaddedAbs
			if(j .le. npad) then
				paddedAbs(j) = 0.0
			else if(j .le. ((2*numberOfSpeedPoints)+npad)) then
				paddedAbs(j) = absLine(j-npad)
			else
				paddedAbs(j) = 0.0
			end if

			If(paddedAbs(j) .gt. maxAbs) then                         ! Normalize max absorption to 1
				maxAbs = paddedAbs(j)
			End if

		End do

!		If(norm .eqv. .true.) then
!	                Do j = 1, npaddedAbs
!				paddedAbs(j) = paddedAbs(j)/maxAbs
!			End do
!               End if

		Do j = 1, npaddedAbs
!                        Write(*,*) i, paddedAbs(j), maxAbs
!			Write(2000+i,*) Sp(j), absLine(j) 
		End do

!******************************************
! Convert from absorption to FM absorption
!******************************************

                call intAbsToFM(FMabs, nFMPoints, modStep, paddedAbs)

		Write(FMfile,'("FMabs_",I3.3,".txt")')i				! Generate filename
		Open(unit=3000+i, file=outpath//FMfile)				! Open file

		maxFM = 0.0							! Set max initially to zero

                Do j = 1, nFMPoints
                        if (FMabs(j) .gt. maxFM) then
                                maxFM = FMabs(j)				! Find max
                        endif
                End do



		if(norm .eqv. .true.) then
			Do j = 1, nFMPoints
				FMabs(j) = FMabs(j)/maxFM			! Normalize FM output to max
			End do
		end if
	
		Do j = 1, nFMPoints
	 		Write(3000+i,*) Sp(j), FMabs(j)				! Write output
		End do
		Close(3000+i)

	

	End do
	

	End Program Absorption
