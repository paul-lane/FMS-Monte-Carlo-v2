include "Constants/mathConstants.f90"


	Program PickingSteps

! A program to combine multiple transverse speed bins in order to make them
! the appropriate size for conversion into FM lineshapes
! Step 2 in data processing

        use mathConstants

	Implicit none

!**************************************
! Define variables
!**************************************
	
	Double Precision :: kTransition, transFreq 	
	Double Precision :: stepSize, DopplerSize
	Double Precision :: lowlimit, uplimit, bin
	Double Precision, dimension(:), allocatable :: speed, amplitude, midPointSpeed, freq, fvv2_2v
	

	Double Precision, dimension(:), allocatable :: sampledSpeed, outputLine
	Integer :: i, k, roundedDoppler, nspeeds, binsize, nbins, remove, ninputs
	Integer :: sampled, j, n

	Character(15) :: filename
	Character(18) :: fileout
	Character(20), Parameter :: filepath1 = "Outputs/Processing1/" 	
	Character(8), Parameter :: filepath2 = "Outputs/"

!**********************************
! READ INPUT PARAMETERS
!**********************************

	Open(unit=9,file='Inputs/fmInputs.inp')

	Read(9,*) lowlimit								! first bin
	Read(9,*) uplimit								! last bin 
	Read(9,*) bin									! binsize
	Read(9,*) stepsize								! laser step size in MHz
	Read(9,*) kTransition								! Transition frequency in cm-1

	Close(9)

	Open(unit=10,file='Inputs/binParameters.inp')
	Read(10,*)
	Read(10,*)
	Read(10,*)
	Read(10,*) ninputs								! Number of input files
	Close(10)

	Write(*,*) "Binning in Frequency Steps"
!*******************************************************
! GET DOPPLER SHIFT SIZES FOR GIVEN FREQUENCY STEPSIZE
!*******************************************************

	binsize = int(bin)								! Convert binsize to integer

	nspeeds =int((uplimit-lowlimit)/binsize)+1					! Calculate number of speeds based on bin width

	transFreq = kTransition*100*speedLight						! Convert transition frquency from cm-1 to Hz

	DopplerSize = (stepSize*1e6*speedLight)/transFreq				! Convert laser step size into Doppler shift in ms-1

!	roundedDoppler = floor(DopplerSize)						! Round Doppler shift to nearest integers (data is in 1ms-1 bins)
	roundedDOppler = nint(DopplerSize)

	nbins = int((1.0*nspeeds)/DopplerSize)+1						! Calculate number of bins


!**************************************
! Allocate arrays
!**************************************

	allocate(speed(nspeeds), amplitude(nspeeds), midPointSpeed(nSpeeds),freq(nSpeeds),fvv2_2v(nSpeeds))

! Calculate the remainder to ensure we don't have a partially filled bin at the end 

	remove = mod(nspeeds,nbins)							! Calculate the remainder when speed data is broken up in bins 

        Write(*,*) uplimit, roundedDoppler, nbins


	Allocate(sampledSpeed(nbins))
	Allocate(outputLine(nbins))
											! This remainder is removed so that we don't have a partially filled bin at the end
!**************************************
! Read File
!**************************************

	Do i = 1, ninputs								! For every input file...
		Write(filename,'("Output_",I3.3,".txt")')i             			! Generate file names to read
		Open(unit=1000+i,file=filepath1//filename)            
		Read(1000+i,*)								! Read header line
		
		Do k = 1, nspeeds 							! For every speed point in the file...
			Read(1000+i,*) speed(k), amplitude(k)				! Read the speed point and it's amplitude
		End do

                Write(fileout,'("SpeedDistr_",I3.3,".txt")')i
                Open(unit=2000+i, file=filepath2//fileout)
                Write(2000+i,*) roundedDoppler, int(nSpeeds/roundedDoppler)

		sampled=0								! Set sampled = 0
		j=0									! Set j = 0 each time the loop starts
		Do While (sampled .le. (uplimit-roundedDoppler))			! For the range of Doppler limits...
			n = j+1								! n is the index of the array where the data is to be stored
			sampled = nint((j+0.5)*DopplerSize)				! Calculate the sampled speeds (first point is 0.5 step for zero)
!		Write(*,*) i, j, sampled, uplimit-roundedDoppler
				sampledSpeed(n) = speed(sampled)			
				outputLine(n) = amplitude(sampled)
!!!!!				Write(2000+i,*) speed(sampled), amplitude(sampled)
				Write(2000+i,*) speed(sampled), outputLine(n)
			j=j+1
		End do 
		Close(1000+i)	


	End do
	End Program PickingSteps
	
 	
