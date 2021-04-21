	Program TimeSliceBinning

! Program to process speed data as a function of time
! Sums the transverse speeds up over a given time interval
! Assumes input data is in 1 us bins
! Step 1 in data processing

	Implicit none

!*********************************
! DEFINE VARIABLES
!*********************************
	
	Integer:: nfiles, noutputs
	Integer:: dataStart, dataEnd, time
	Integer:: outStart, outEnd
	Integer:: nSpeeds, nOutput
	Integer:: i, j, k


	Double Precision, dimension(:), allocatable :: start, end 	
	Double Precision :: lowlimit, uplimit, step
	Double Precision, dimension(:), allocatable :: speed, amplitude, summedAmplitude
	Character*30 filename
	Character*15 outfile
	Character(20), PARAMETER :: filepath = "Outputs/Processing1/"


!*********************************
! READ INPUT FILES
!*********************************

! Read fmInputs file

	Open(unit=9,file="Inputs/fmInputs.inp")
	Read(9,*) lowlimit							! Read lower limit of Speed distribution
	Read(9,*) uplimit							! Read upper limit of Speed distribution
	Read(9,*) step								! Read step size of Speed distribution

	nSpeeds = int((uplimit-lowlimit)/step)+1				! Calculate number of speeds

! Read bin parameters file

	Open(unit=10,file="Inputs/binParameters.inp")
	Read(10,*) 
	Read(10,*) dataStart, dataEnd						! Start and end speeds in filenames
	Read(10,*)
	Read(10,*) noutputs							! Number of outputs required

	Write(*,*) "Binning Time Slices"

!**********************************
! Allocate arrays
!**********************************

	Allocate(start(noutput), end(noutput))
	Allocate(speed(nSpeeds), amplitude(nSpeeds),summedAmplitude(nSpeeds))
	
	nfiles = dataEnd - dataStart + 1					! Calculate total no of files

!***********************************
! GENERATE OUPUT FILENAMES
!***********************************

! For each output listed generate a file to write to

	Do i = 1, noutputs
		write(outfile, '("Output_",I3.3,".txt")')i
		Open(unit=1000+i,file=filepath//outfile)		


		summedAmplitude = 0.0						! Initialize value to zero at start of each iteration

! Read it's start and end time value

		Read(10,*) outStart, outEnd					! These are the start and end times of the binning 	
		Write(1000+i,*) outStart, outEnd				! Write these as a header for the output files


!***********************************
! CONSTRUCT SPEED INPUT FILENAMES
!***********************************

		Do j = 1, (outEnd-outStart+1)					! Make the appropriate filename/path for the Speed file for each time to be read in 
			time = outStart + (j-1)
			write(filename, '("Outputs/Speed/Speed_",I3.3,".txt")')time

!******************************************
! Open Transverse Speed file for this time
!******************************************

			Open(unit=100+j,file=filename)	


			Do k = 1, nSpeeds					! For each speed create a summed amplitude over the time range
				read(100+j,*) speed(k), amplitude(k)		! containing all the times up to this point in the bin
				summedAmplitude(k) = summedAmplitude(k) + amplitude(k)
			End do 	

			Close(100+j)
		End do

!*****************************************************************************
! Once the loop has finished all the files, evaluate write the output to file		
!*****************************************************************************

		Do k = 1, nSpeeds
			Write(1000+i,*) speed(k), summedAmplitude(k) 
		End do
	Close(1000+i)
	End do



	End Program TimeSliceBinning
