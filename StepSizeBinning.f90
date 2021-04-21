include "Constants/mathConstants.f90"


	Program StepSizeBinning

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
	Integer :: i, k, roundedDoppler, nspeeds, binsize, nbins, remove, ninputs
	Double Precision :: fv 						
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

	roundedDoppler = floor(DopplerSize)						! Round Doppler shift to nearest integers (data is in 1ms-1 bins)

	nbins = roundedDoppler/binsize							! Calculate number of bins

!**************************************
! Allocate arrays
!**************************************

	allocate(speed(nspeeds), amplitude(nspeeds), midPointSpeed(nSpeeds),freq(nSpeeds),fvv2_2v(nSpeeds))

! Calculate the remainder to ensure we don't have a partially filled bin at the end 

	remove = mod(nspeeds,nbins)							! Calculate the remainder when speed data is broken up in bins 
											! This remainder is removed so that we don't have a partially filled bin at the end

!**************************************
! Read File
!**************************************
	Do i = 1, ninputs								! For every input file
		
		Write(filename,'("Output_",I3.3,".txt")')i				! Generate it's filename
		Open(unit=1000+i,file=filepath1//filename) 				! Open it
		Read(1000+i,*)								! Ignore header
		Do k = 1, nspeeds  							! For each speed point read in speed and amplitude of transverse speed distrib
			Read(1000+i,*) speed(k), amplitude(k)
			midPointSpeed(k) = speed(k) + (roundedDoppler/2)		! Generate the midpoint speed of the bin 

!!			fvv2_2v(k) = (amplitude(k)*(midPointSpeed(k)**2))/(2*midPointSpeed(k))		! Calculate f(v)v^2/2v for each point

			fvv2_2v(k) = amplitude(k)/(2*midPointSpeed(k))			! Speed distribution is f(v)v^2, 1/2v is the Jacobian
		End do
		Close(1000+i)	

	
!**************************************
! Do binning and Write output
!**************************************
			
		Write(fileout,'("SpeedDistr_",I3.3,".txt")')i
		Open(unit=2000+i, file=filepath2//fileout)				! Create filename
		Write(2000+i,*) roundedDoppler, ((nSpeeds-remove)/nbins)		! Write header containing doppler stepsize and number of speed points

	
		Do k = 1, (nspeeds-remove), nbins					
			fv = sum(fvv2_2v(k:k+nbins-1))					! Integration of f(v)v^2/2v (as with of bin is 1ms-1 dt = 1ms-1 therefore *1 ignored)
			Write(2000+i,*) midPointSpeed(k), fv, amplitude(k)		! Write to file 
		End do
		Close(2000+i)
	End do
	
	End Program StepSizeBinning
	
 	
