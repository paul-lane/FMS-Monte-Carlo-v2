! ********************************************************************************************
! Program to calculate the integral and FM absorption from the transverse speed distributions
! ********************************************************************************************

	include "Maths/fmDetection.f90"

	Program Absorption

	use fmDetection

	Implicit None

! ****************************************
! Define variables
! ****************************************

	Integer:: i, j, nfiles, npad, nAbs, nFMPoints, numberOfSpeedPoints
	Double Precision :: lowlimit, uplimit, stepSize, modFreq, modStep, speedStep, min, maxAbs, maxFM, minX
	Double Precision, allocatable :: freq(:), amplitude(:)
	Double Precision, Allocatable :: transSpeed(:), FMabs(:), absLine(:), sp(:), spX(:)
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
		nFMPoints = nAbs - (6 * modStep)				! Calculate number of points in FM lineshape

		If (Allocated(absLine)) then					! Deallocate all arrays from previous iteration (if one is allocated all were)
			Deallocate(absLine)					! Memory issues are created if these aren't deallocated at the start of the loop
			Deallocate(transSpeed)
			Deallocate(freq)
			Deallocate(amplitude)
			Deallocate(FMabs)
			Deallocate(Sp)
			Deallocate(SpX)
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



		Do j = 1, numberOfSpeedPoints					! For each speed
			Read(1000+i,*) transSpeed(j), amplitude(j)		! Read in transverse speed (midpoint of bin) and amplitude of f(v)v^2/2v

		End do
		Close(1000+i)

!*****************************************
! Convert from speed to absoprtion line
!*****************************************

		call speedToAbsorptionLine(amplitude, numberOfSpeedPoints, nPad, nAbs, speedStep, absLine)

		Write(absFile,'("Abs_",I3.3,".txt")')i				! Generate filename
		Open(unit=2000+i, file=outpath//absFile)			! Open file

		minX = -1* speedStep * (numberOfSpeedPoints + nPad) 	
		maxAbs = 0.0							! Initialize max to zero

		Do j = 1, nAbs
			if (absLine(j) .gt. maxAbs) then	
				maxAbs = absLine(j)				! Find max value
			endif
		End do


		if (norm .eqv. .true.) then 
			Do j = 1, nAbs
				absLine(j) =  absLine(j)/maxAbs			! if normalizing output modify absLine
			End do
		end if

	
		Do j = 1, nAbs							! for each absorption point
			SpX(j) = minX + ((j+1)*speedStep) - (speedStep/2) 	! Generate x axis
			Write(2000+i,*) Spx(j), absLine(j)			! write output
		End do
		
		Close (2000+i)


!******************************************
! Convert from absorption to FM absorption
!******************************************
		call intAbsToFM(FMabs, nFMPoints, modStep, absLine)	

		Write(FMfile,'("FMabs_",I3.3,".txt")')i				! Generate filename
		Open(unit=3000+i, file=outpath//FMfile)				! Open file

		min = -1 * speedStep * numberOfSpeedPoints			! Generate x axis
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
			Sp(j) = min + ((j+1)*speedStep) - (speedStep/2)		! Generate x axis values
	 		Write(3000+i,*) Sp(j), FMabs(j)				! Write output
		End do
		Close(3000+i)

	

	End do
	

	End Program Absorption
