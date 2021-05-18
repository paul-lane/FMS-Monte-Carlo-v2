! **************************************************************************************************************
! Program to calculate the FM signals for gas-liquid scattering experiments
! Calculates the appearance profiles and transverse velocities as a function of time for a multipass
! laser beam experiment
! ************************************************************************************************************** 

! List of included module paths

include "Constants/mathConstants.f90"
include "Maths/speeds_v2.f90"
include "Maths/directions.f90"
include "Maths/fmDetection.f90"

program MCScatteringFM

! Tell the code to use the specified modules
    use speeds_v2
    use directions
    use mathConstants
    use fmDetection
 
    implicit none

! *********************************************************************************************************************
! DEFINE VARIABLES
! *********************************************************************************************************************

! Variables from directions input (defined in input file)
    Integer :: ncyc, cosinePowerTD, cosinePowerIS
    Double Precision :: incidenceAngle, x0, aMax, aMin, h, s, dist, massMol, energyTrans, surfaceMass, exitAngle, temp
    Double Precision :: maxSpeed, scatterFraction, mass
    Double Precision :: gMean, gStdDev, lgFraction, lGamma
    Double Precision :: weight1, weight2, weight3, sumweight, w1, w2, w3
    Real :: m1, m2, m3, std1, std2, std3
    Logical :: scattering, outputIngoing, xyOutput, blurOn, GaussSpeeds

! Variables from chamber dimensions input (defined in input file)
    Double Precision :: skimPos, valvePos, colPos, skimRad, valveRad, colRad, probeCentre 
    Double Precision :: wheelCentrex, wheelCentrey, wheelRadius, wheelBathy	! All hit wheel subroutine parameters
    Double Precision :: bathTop							! height (y) of top of bath

! Variables from timing inputs (defined in input file)
    Double Precision :: pulseLength, probeStart, probeEnd, tStep

! Variables from FM inputs (defined in input file)
    Double Precision :: speedStart, speedEnd, speedStep

! Variables from Heriott cell inputs
    Integer :: npass								! number of passes
    Double Precision :: rLaser							! radius of laser beam
    Double Precision, dimension(:), allocatable:: yOffset, zOffset, th1, th2 	! offset positions from cell centre in y and z, Euler angles     

! Variables from Overwrite input (defined in file)
    logical :: owIngoingSp, owIngoingPos, owIngoingTraj
    logical :: owOutgoingSp, owOutgoingPos, owOutgoingTraj
    Double Precision :: owSpIn, owPosxIn, owPosyIn, owTrajxIn, owTrajyIn, owTrajzIn
    Double Precision :: owSpOut, owPosxOut, owPosyOut, owTrajxOut, owTrajyOut, owTrajzOut

! Vectors/Arrays Produced in code
    Double Precision, dimension(2,3) :: particleVector, particleStartPos        ! particle trajectory vector and start position in the lab frame
    Double Precision, dimension(2,3) :: EulerVector, EulerStartPos              ! particle trajectiry vector and start position in the Euler frame
    Double Precision, dimension(3) :: probeCentreVector, EulerCellCentre        ! vector position of the centre of the probe beam in the lab and Euler frames
    Double Precision, dimension(2) :: particleSpeed, particleTime               ! particle (scalar) speed and creation time
    Double Precision, dimension(2,3) :: Pos, LabPos                             ! positions where the trajectory intersects the probe beam in the Euler and Lab frames

! Output Arrays
    Double Precision, dimension(:), allocatable:: AppProfile			! Appearance Profile array
    Double Precision, dimension(:,:), allocatable:: TransSpeed			! 2D array for storing transverse speeds as a function of time

! Loop counters and array sizes
    Integer:: i, j, k								! Loop iterators
    Integer:: NumberOfTimePoints, NumberOfSpeedPoints				! No of time points in App Profile, No of speed points in transverse speed
    Integer:: vectorsPerParticle, startParticle						! No of vectors run for a particular particle (1 if ingoing only, 2 if scattering)	

! Variables for counter print out
    Integer:: acceptedCounter							! Count of accepted particle trajectories  
    Integer:: acceptedIngoing, acceptedScattered				! Count of accepted ingoing trajectories and scattered trajectories
    Integer:: totalTraj								! Total number possible trajectories given each particle could be probed by each beam
    Double Precision:: acceptanceRatio						! Ratio of accepted trajectories to total
    Double Precision:: startTime, endTime, runTime				! CPU start time, end time and calculation time   
    Integer, dimension(:,:), allocatable :: arrayCount				! Counts trajectories that intersect each laser beam
    Integer:: missWheel								! Counts the trajcetories that miss the wheel

! Variables for filename production
    Integer:: filetime								! Variable which provides the number (time) in the Speed filename
    Character*50 :: speedFilename						! Speed filename

! Variables used for profiling the beam
    Double Precision, dimension(2) :: profile					! Beam profile at a given z distance
    Double Precision :: z							! z position beam profile is taken at

! Probe/Intersection parameters
    Logical:: intersect, hit							! Particle intersects beam, Particle hits wheel
    Double Precision:: entryTime, exitTime					! Time particle enters and exits laser beam
    Integer:: startTimePoint, endTimePoint                                      ! The time bin that the particle enters and exits the beam (entry and exit times mapped
										! onto the App Profile grid of points)

! Scattering parameters 
    Double Precision:: mostLikelyProbability					! Most likely probability for Maxwell-Boltzmann
    Double Precision:: tWheel							! Time taken for ingoing particle to reach the wheel (time of flight)
    Double Precision:: deflectionAngle						! Deflection angle in soft sphere scattering (IS)
    Double Precision:: rand1							! Random number for IS/TD ratios
    Logical:: correctDirection							! correct direction (IS)

! Random number generation
   Double Precision :: rand
   Integer :: i_seed
   Integer, Dimension(:), Allocatable :: a_seed
   Integer, Dimension(1:8) :: dt_seed

! ****************************************
! INITIALIZE COUNTERS
! ****************************************

    acceptedCounter = 0
    acceptedIngoing = 0
    acceptedScattered = 0
    missWheel = 0

!*****************************************
! GENERATE RANDOM SEED
!*****************************************

    ! Without calling random seed, random number sequences can often be repeated

    call random_seed(size=i_seed)
    allocate(a_seed(1:i_seed))
    call date_and_time(values=dt_seed)
    a_seed(i_seed)=dt_seed(8); a_seed=dt_seed(8)*dt_seed(7)*dt_seed(6)
    call random_seed(put=a_seed)
    deallocate(a_seed)

    call cpu_time(startTime)

!******************************************
! LOAD INPUT PARAMETERS
! General input parameters (no of particles, chamber geometery etc)
! Note that to avoid any confusion/compatibility issues this is the 
! same input file as for LIF imaging so some parameters may be redundant
!******************************************

!******************************************
! READ IN DIRECTIONS INPUTS
!******************************************
     Open(unit=11,file='Inputs/directionInputs_v2.inp')
     Read(11,*) incidenceAngle						! Angle of incidence (not used)
     Read(11,*) ncyc                                                    ! Number of trajectories
     Read(11,*) GaussSpeeds						! Use Gaussian Speeds T/F 
     Read(11,*) weight1, weight2, weight3				! Weighting of 3 Gaussians
     Read(11,*) m1, m2, m3						! Mean of the 3 Gaussians
     Read(11,*) std1, std2, std3					! Standard deviation of the 3 Gaussians
     Read(11,*) dist							! Valve-Probe Distance at which the characterised apperance profile was recorded	
     Read(11,*) x0							! }
     Read(11,*) aMax							! }
     Read(11,*) aMin							! } Parameters for speed distribution from Origin function of ingoing app profile
     Read(11,*) h							! }
     Read(11,*) s							! }
     Read(11,*) outputIngoing                                           ! Output ingoing profile? T/F
     Read(11,*) xyOutput                                                ! Output x-y co-ordinates at various points
     Read(11,*) scattering						! Perform scattering? T/F
     Read(11,*) scatterFraction                                         ! Fraction of scattering which is IS and TD
     Read(11,*) massMol							! Molar mass of incident particle
     Read(11,*) energyTrans						! Fraction of energy transfered to the surface
     Read(11,*) surfaceMass						! Effective surface mass
     Read(11,*) exitAngle						! Exit angle
     Read(11,*) cosinePowerIS                                           ! IS power of cosine of scattering angle (cos^n)
     Read(11,*) temp							! Surface temperature
     Read(11,*) maxSpeed						! Maximum speed for TD
     Read(11,*) cosinePowerTD						! TD power of cosine of scattering angle (cos^n)
     Read(11,*) blurOn							! Turn on Gaussian Blurring
     Read(11,*) lgFraction                                              ! Fractional weighting of Gaussian/Loretzian in transverse temp blurring
     Read(11,*) gMean							! Mean of Gaussian for ingoing beam transverse temp blurring	
     Read(11,*) gStdDev							! Std Dev of Gaussian for ingoing beam transverse temp blurring
     Read(11,*) lGamma							! Gamma for Lorentzian in transverse temp blurring
     Close(11)

     sumWeight= weight1 + weight2 + weight3
     w1 = weight1/sumWeight
     w2 = weight2/sumWeight
     w3 = weight3/sumWeight



     mass = massMol/avagadro

!******************************************
! READ IN CHAMBER DIMENSIONS
!******************************************
     Open(unit=12, file='Inputs/chamberDimensions.inp')
     Read(12,*) skimPos							! skimmer position 
     Read(12,*) valvePos						! valve position
     Read(12,*) colPos							! collimator position
     Read(12,*) skimRad							! skimmer radius
     Read(12,*) valveRad						! valve radius
     Read(12,*) colRad							! collimator radius
     Read(12,*) probeCentre						! centre of Herriot cell
     Read(12,*) wheelCentrex						! centre of wheel in x direction
     Read(12,*) wheelCentrey						! centre of wheel in y direction (wheel axle is offset with respect to chamber axis)
     Read(12,*) wheelRadius						! radius of wheel
     Read(12,*) wheelBathy						! distance between wheel axle and top of bath 
     Close(12)

     bathTop = wheelCentrey - wheelBathy				! position of the top of the bath (used for hit wheel subroutine

!****************************************
! READ IN TIMING INPUTS
!****************************************
     Open(unit=13, file='Inputs/timingInputs.inp')
     Read(13,*) pulseLength						! Disharge pulse length (used??)
     Read(13,*) probeStart						! Probe start time (earliest point in App profile)
     Read(13,*) probeEnd						! Probe end time (latest point in App Profile)
     Read(13,*) tStep							! App Profile time step size
     Close(13)

!****************************************
! READ IN FM PARAMETERS
!****************************************
 
! Read in Doppler and FM parameters
     Open(unit=14,file='Inputs/fmInputs.inp')
     Read(14,*) speedStart						! Lowest speed in Transverse speed profile
     Read(14,*) speedEnd						! Fastest speed in Transverse speed profile
     Read(14,*) speedStep						! Speed step size of transverse speed profile

! Note that there are more lines to tis input file containing laser step size, transition frequency and modulation frequency
! but these are not required for this code so are not called...

     Close(14)

!****************************************
! Read in Herriott Cell parameters
!****************************************

     Open(unit=15,file='Inputs/herriottCellGeom.inp')
     Read(15,*) npass							! number of laser passes
     Read(15,*) rLaser							! radius of laser beam
     Read(15,*)
     allocate(yOffset(npass), zOffset(npass), th1(npass), th2(npass))	! Allocate arrays of the appropriate size

     allocate(arrayCount(2,npass))
     arrayCount=0

     Do k = 1, npass					
	Read(15,*) yOffset(k), zOffset(k), th1(k), th2(k)		! Read in internal dimensions of Herriott Cell
     Enddo	 
     Close(15)

!****************************************
! Read in Overwrite parameters
!****************************************

     Open(unit=16,file='Inputs/overwrite.inp')
     Read(16,*) owIngoingSp						! Overwrite Ingoing Speed? T/F
     Read(16,*) owSpIn							! Overwritten Ingoing Speed value
     Read(16,*) owIngoingPos						! Overwrite Ingoing start position? T/F
     Read(16,*) owPosxIn, owPosyIn					! Overwritten start positions in x & y
     Read(16,*) owIngoingTraj						! Overwrite Ingoing trajectory vector? T/F
     Read(16,*) owTrajxIn, owTrajyIn, owTrajzIn				! Overwritten values of trajectory vector in x, y & z
     Read(16,*) owOutgoingSp						! Overwrite Scattered speed? T/F
     Read(16,*) owSpOut							! Overwritten Scattered speed value
     Read(16,*) owOutgoingPos						! Overwrite Scattered start position? T/F
     Read(16,*) owPosxOut, owPosyOut					! Overwritten values of scattered start position in x & y
     Read(16,*) owOutgoingTraj						! Overwrite Scattered trajectory vector? T/F
     Read(16,*) owTrajxOut, owTrajyOut, owTrajzOut			! Overwritten values of scattered trajectory vector in x, y & z
     Close(16)

!**************************************
! Allocate arrays & initialize arrays
!**************************************

    NumberOfTimePoints = ((probeEnd - probeStart) / tStep) + 1		! Calculate number of time points 
    NumberOfSpeedPoints = ((speedEnd - speedStart) / speedStep) + 1	! Calculate number of speed points

    allocate(AppProfile(NumberOfTimePoints))                            ! This is an array for appearance profile intensity
    allocate(TransSpeed(NumberOfTimePoints,NumberOfSpeedPoints))        ! This is a 2D array holding transverse velocity info as a function of time

    AppProfile = 0.0							! Initialize arrays to zero
    TransSpeed = 0.0
    

!************************************************
! SIMULATION DETAILS
! Scattered trajectories are calculated if the input flag is set to true
! Also have the option to not output ingoing beam if required
!  
!************************************************

    if ((outputIngoing .eqv. .false.) .and. (scattering .eqv. .false.)) then
        Write(*,*) "You have chosen to output neither ingoing or scattering information"
        Write(*,*) "DON'T WASTE MY TIME"
 	STOP
    end if


    if (outputIngoing .eqv. .true.) then                                ! Have the ability to output only the scattered profile
        startParticle = 1
    else
        startParticle = 2
    end if



    if (scattering .eqv. .TRUE.) then

	mostLikelyProbability = MB_most_likely(temp, mass)		! Calculate probability of most probable speed at a given temp for TD subroutines
        vectorsPerParticle = 2						! If doing scattering sets the number of vectors per particle to 2 - used for loops later
	if(outputIngoing .eqv. .false.) then
		Write(*,*) "Simulating Scattered Beam Only"
	else
		Write(*,*) "Simulating Ingoing and Scattered Beams"
	end if 
    else
        vectorsPerParticle = 1						! Otherwise number of vectors per particle is set to 1 (ingoing simulation only)
	Write(*,*) "Simulating Ingoing Beam Only"
    end if

!********************************************************************************
! OPEN FILES TO WRITE STEP-BY-STEP DATA TO
!********************************************************************************

    Open(unit=22,file='Outputs/MissedWheel.txt')
    Write(22,*) "List of x-y co-ordinates of trajectories that miss the wheel"		! Most of the time everything should hit, shows that file has been written

    Open(unit=23,file='Outputs/IngoingCollimatorProfile.txt')
    Write(23,*) "List of x-y co-ordinates of trajectories exiting the collimator"

    Open(unit=24,file='Outputs/IngoingCellProfile.txt')
    Write(24,*) "List of x-y co-ordinates of ingoing trajectories at the centre of the probe plane (not all will be detected)"

    Open(unit=25,file='Outputs/ScatteredCellProfile.txt')
    Write(25,*) "List of x-y co-ordinates of scattered trajectories at the centre of the probe plane (not all will be detected)"

    Open(unit=26,file='Outputs/WheelProfile.txt') 
    Write(26,*) "List of x-y co-ordinates of trajectories at the wheel plane (not all will hit the wheel)"

    Open(unit=27,file='Outputs/IngoingProbed.txt')
    Write(27,*) "List of x-y co-ordinates of ingoing probed trajectories at the centre of the Herriott Cell)"

    Open(unit=28,file='Outputs/ScatteredProbed.txt')
    Write(28,*) "List of x-y co-ordinates of ingoing probed trajectories at the centre of the Herriott Cell)"

    if(xyOutput .eqv. .false.) then
	Do i = 22, 28
 	 	Write(i,*) "This file is blank as you have chosen not to output beam profiles"
	End do 
    endif


!********************************************************************************
! CREATES TRAJECTORIES 
! for ingoing and scattered trajectories as appropriate
!********************************************************************************

    do i = 1, ncyc							! For 1 to number of trajectories

	! Generate ingoing particle trajectories and speeds

	if(GaussSpeeds .eqv. .TRUE.) then				! Use sum of Gaussians
		call ingoing_speed_from_Gauss(w1, w2, w3, m1, m2, m3, std1, std2, std3, dist, pulseLength, particlespeed(1), particleTime(1))
	else								! Use Origin function
        	call ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, particleSpeed(1), particleTime(1))
	end if

        call ingoing_direction(valveRad, valvePos, skimRad, skimPos, colRad, colPos, particleVector(1,:), particleStartPos(1,:))

!**************************************************
! TRANSVERSE TEMPERATURE BLURRING OF INGOING BEAM
!**************************************************
! NEW AND UNTESTED... There are some additional parameters needed as the subroutine is the sum of a Gaussian and Lorentian

	if(blurOn .eqv. .true.) then

	call transverse_temp(gMean, gStdDev, lGamma, lgFraction, colPos, (valvePos-colPos), particleTime(1), &
		particleSpeed(1), particleStartPos(1,:), particleVector(1,:))

	endif


!**************************************************
! OVERWRITE INGOING INPUTS IF DESIRED
! ingoing speeds, vectors and start positions
! are overwritten here if selected in inputs
! Note warnings are written to screen
!**************************************************

	if (owIngoingSp .eqv. .TRUE.) then
		Write(*,*) '*** Ingoing Speed Overwritten ***'
		particleSpeed(1) = owSpIn
	end if

	if (owIngoingTraj .eqv. .TRUE.) then
		Write(*,*) '*** Ingoing Trajectory Overwritten ***' 	
		particleVector(1,1) = owTrajxIn
		particleVector(1,2) = owTrajyIn
		particleVector(1,3) = owTrajzIn
	end if

	if (owIngoingPos .eqv. .TRUE.) then
		Write(*,*) '*** Ingoing Start Positions Overwritten ***'
		particleStartPos(1,1) = owPosxIn
		particleStartPos(1,2) = owPosyIn
	end if


!****************************************
! WRITE BEAM PROFILES TO FILE
! Need to double check this subroutine as it has been a while
!****************************************

	if(xyOutput .eqv. .true.) then

! If you ever want to write out a list of start positions and trajectories this is the place to do it
! I did it to debug but there is no need to write these files every time

! Write Ingoing Beam Profile at collimator position
		z = abs(particleStartPos(1,3) - colPos)  
		call getProfile(particleVector(1,:), particleStartPos(1,:), z, profile)
		Write(23,*) profile(1), profile(2)

! Write Ingoing Beam Profile at centre of probe beam 
        	z = abs(particleStartPos(1,3) - probeCentre)
		call getProfile(particleVector(1,:), particleStartPos(1,:), z, profile)
		Write(24,*) profile(1), profile(2)

! Write Scattered Beam Profile at centre of probe beam
        	z = abs(particleStartPos(2,3) - probeCentre)
        	call getProfile(particleVector(2,:), particleStartPos(2,:), z, profile)
        	Write(25,*) profile(1), profile(2)

! Write Beam profile at wheel (note as scattered trajectories start from where the ingoing beam hits the wheel this is the same for both)
        	z = abs(particleStartPos(1,3) - 0.0)
		call getProfile(particleVector(1,:), particleStartPos(1,:), z, profile)
  		Write(26,*) profile(1), profile(2)
	end if

! *****************************************
! ANGLE OF INCIDENCE
! changes the angle of incidence and starting point of the particle using a rotation matrix
! This is not needed for FM setup but will keep in anyway as long as angle of incidence is set to zero 
! or these next two lines are commented out angle of incidence will be zero
! *****************************************

! 	 call rotation(particleVector(1,:), incidenceAngle, particleVector(1,:))
!        call rotation(particleStartPos(1,:), incidenceAngle, particleStartPos(1,:))


! *****************************************
! TIME TAKEN TO TRAVEL TO WHEEL
! (NOT time of origin for scattered particle)
! NOTE: as the measurement is in the lab frame, distance travelled to wheel depends only on z-positions
! *****************************************

        tWheel = abs(particleStartPos(1,3) / (particleSpeed(1)*particleVector(1,3)))
        

! ******************************************
! SCATTERED PARTICLE PARAMETERS
! Creates scattered particle parameters based on ingoing beam particle
! The time a particle is scattered depends on time taken to reach the wheel and creation time of incoming particle
! Scattered start position is where the incoming particle hit the wheel (calculated from its initial traj and pos)
! ******************************************

        particleTime(2) = particleTime(1) + tWheel
        particleStartPos(2,1) = particleStartPos(1,1) + (particleVector(1,1)*tWheel*particleSpeed(1))
        particleStartPos(2,2) = particleStartPos(1,2) + (particleVector(1,2)*tWheel*particleSpeed(1))
        particleStartPos(2,3) = 0

        if (scattering .eqv. .TRUE.) then

!**************************** This bit is Adam's new IS / TD scattering section ******************************************
! IT IS COPIED EXACTLY FROM ADAM'S CODE AND HAS NOT UNDER GONE TESTING YET

		call random_number(rand1)

		! CASE 1: TD SCATTERING
		if(rand1 .gt. scatterFraction) then	
			call MB_speed(maxSpeed, temp, mass, mostLikelyProbability, particleSpeed(2))	! Generates speed from a Maxwell-Boltzmann distribution 
			call cosine_distribution(cosinePowerTD, particleVector(2,:))			! Generates particle distribution according to cos^n

		! CASE 2: IS SCATTERING
		else
			correctDirection = .false.
			
			! rejects particles not scattering in positive z direction from surface
			do while (correctDirection .eqv. .false.)
				! sets impulsive scattering direction based on some cosine direction in IS subroutine
				call cosine_distribution(cosinePowerIS, particleVector(2,:))
				! rotates scattered vector about the y-axis (this may not represent scattered distribution properly)
				call rotation(particleVector(2,:), exitAngle, particleVector(2,:))

				if(particleVector(2,3) .gt. 0) then
					correctDirection = .true.
				end if
			end do

			! Sets IS speed based on scattered direction using soft sphere model
			call deflection_angle(particleVector(1,:), particleVector(2,:), deflectionAngle)
			call soft_sphere_speed(massMol, energyTrans, surfaceMass, particleSpeed(1), deflectionAngle, particleSpeed(2))
		end if

		
!*********************************************	
! OVERWRITE SCATTERED PARAMETERS IF DESIRED
!*********************************************

	if (owOutgoingSp .eqv. .TRUE.) then
		Write(*,*) '*** Scattered Speed Overwritten ***'
		 particleSpeed(2) = owSpOut
	end if

	if (owOutgoingTraj .eqv. .TRUE.) then
		Write(*,*) '*** Scattered Trajectory Vector Overwritten ***'
		particleVector(2,1) = owTrajxOut
		particleVector(2,2) = owTrajyOut
		particleVector(2,3) = owTRajzOut
	endif

	if (owOutgoingPos .eqv. .TRUE.) then
		Write(*,*) '*** Scattered Trajectory Start Positions Overwritten ***'
		particleStartPos(2,1) = owPosxOut
		particleStartPos(2,2) = owPosyOut
	endif


        end if 											! Ends if scattering .eqv. true

!*************************************************************************************************
! USING THE TRAJECTORIES
! By this point all the work has been done to create trajectories, from here on is the probing 
!*************************************************************************************************

! Loops through ingoing trajectories (j=1) then scattered trajectories (j=2)

        do j = startParticle, vectorsPerParticle						! Do from startParticle to number of vectors per particle 


! **************************************************
! MULTI PASS SETUP
! **************************************************
             do k = 1, npass									! Do for each pass in the Herriott cell

	    probeCentreVector(1) = 0.0								! Probe centre position is stored as a vector
	    probeCentreVector(2) = yOffset(k)
	    probeCentreVector(3) = probeCentre + zOffset(k)

 
!**************************************************
! EULER TRANSFORM
! Rotates the vectors into the frame of the probe beam
!**************************************************

            call Euler(particleVector(j,:), th1(k), th2(k), EulerVector(j,:))
	    call Euler(particleStartPos(j,:), th1(k), th2(k), EulerStartPos(j,:))
	    call Euler(probeCentreVector(:), th1(k), th2(k), EulerCellCentre(:))	        
	

!***********************************************
! CHECKS BEAM INTERSECTION
! Finds out if trajectory passes through laser beam
!***********************************************

! Check if there is an intersection between trajectory and laser beam

           call checkIntersect(EulerVector(j,:), EulerStartPos(j,:), rLaser, EulerCellCentre(:), Intersect)


!***********************************************
! HIT WHEEL 
!***********************************************

	    ! Up until this point we have assumed that the wheel is much larger than the spot size of the ingoing molecular beam
	    ! therefore every particle is assumed to hit the wheel and therefoe generate a scattered particle originating at that same position.
	    ! Here we assume a finite size wheel and ensure that the scattered particles that would be generated from an ingoing particle
	    ! which did not hit the wheel are not probed. This is done by overwriting the Intersect parameter which otherwise allows them
	    ! to be probed.  


!	   if(j .eq. 2) then							! This used to only be run for scattered particles but is needed to see which
										! ingoing trajectories miss the wheel...
	   	call hitWheel(particleStartPos(2,:), wheelCentrex, wheelCentrey, wheelRadius, bathTop, hit)
 	    	
		if(hit .eqv. .false.) then					! If we have a trajectory missing the wheel we only add to the counter on the
	     		intersect = .false.					! first laser pass (the trajectory is independent of the probe beam)

			if(k .eq. 1) then					! If we counted every rejection our stats would be a factor of npass to large
				missWheel = missWheel + 1			! Output positions of trajectories which miss the wheel

				if(xyOutput .eqv. .true.) then
					Write(22,*) particleStartPos(2,1), particleStartPos(2,2)
				end if
			end if						
		endif								

!	   endif

!*************************************************
! FIND INTERSECTION DETAILS
!*************************************************

         
            ! If the trajectory intersects the probe beam we calculate xyz entry and exit position and use detection routines
	    ! otherwise it is ignored (note that because the ingoing and scattered trajectories are handled independently
	    ! a scattered trajectory does not necessarily have to come from a probed ingoing trajectory	

            if(Intersect) then


!*************************************************
! WRITE PROBED BEAM PROFILES
!*************************************************
! These provide a 2-D projection of the probe volume onto the centre of the Herriott cell 

	if(xyOutput .eqv. .true.) then

! Write ingoing beam Profile at centre of Herriot Cell (only probed trajectories only)
        	z = abs(particleStartPos(1,3) - probeCentre)
        	call getProfile(particleVector(1,:), particleStartPos(1,:), z, profile)
        	Write(27,*) profile(1), profile(2)

! Write scattered beam profile at centre of Herriot Cell (only probed trajectories)
       		z = abs(particleStartPos(2,3) - probeCentre)
        	call getProfile(particleVector(2,:), particleStartPos(2,:), z, profile)
        	Write(28,*) profile(1), profile(2)

	endif

!*************************************************
! GET INTERSECTION POSITIONS
!*************************************************


! Get the intersection positions 

                call getProbeIntersection(EulerVector(j,:), EulerStartPos(j,:), rLaser, EulerCellCentre(:), Pos)

	
! Note that intersect positions (Pos) are in the frame of the laser beam and need to undergo 
! a reverse Euler transform back into the lab frame 
		
		Call RevEuler(Pos(1,:), th1(k), th2(k), LabPos(1,:))
		Call RevEuler(Pos(2,:), th1(k), th2(k), LabPos(2,:))


!********************************************                
! GET THE ENTRY AND EXIT TIMES FROM THE LASER
!********************************************
! Intersection times (time when particle enters and exits the laser beam are calculated in the lab frame

                call getIntersectionTime(LabPos, particleStartPos(j,:), particleVector(j,:), particleSpeed(j), particleTime(j), entryTime, exitTime) ! Get entry and exit times
		


! We want to bin our data onto a series of descrete points, this routine directs inputs to the correct bin

                call startEndTimePoints(NumberOfTimePoints, entryTime, exitTime, probeStart, probeEnd, tStep, startTimePoint, endTimePoint)

! Perform the probe tasks i.e. create the appearance profile and the transverse velocity profile for each bin in the apperance
! profile. Note that here we are working in two frames: the appearance profile is calcuated in the lab frame while the
! transverse velocity information is calculated in the frame of the laser beam

                call Probe(AppProfile, TransSpeed, NumberOfTimePoints, entryTime, exitTime, startTimePoint, endTimePoint, probeStart, tStep, particleSpeed(j), EulerVector(j,:), EulerStartPos(j,:), speedStart, speedStep)


!************************************************
! OUTPUT COUNTERS
!************************************************

! Count ingoing trajectories
		if(j .eq. 1) then
			acceptedIngoing = acceptedIngoing + 1
		endif

! Count scattered trajectories 
		if(j .eq. 2) then
			acceptedScattered = acceptedScattered + 1
		endif

! Count all probed trajectories
    
                acceptedCounter = acceptedCounter + 1       ! This is a count of the trajectories that are probed

		arrayCount(j,k) = arrayCount(j,k) + 1

            end if ! intersect

	    end do ! npass

        end do	! vectors per particle
    
    end do 	! ncyc 

!************************************************
! CLOSE OUPTUT FILES
!************************************************

      Close(22)
      Close(23)
      Close(24)
      Close(25)
      Close(26)
      Close(27)
      Close(28)


!************************************************
! WRITING DATA TO FILE
!************************************************
	
! Write out Appearance Profiles and Transverse Speed Distribution
	Open(unit=20,file='Outputs/AppProfile.txt')
	Do i = 1, numberOfTimePoints
         filetime = (((i-1)*tstep)+probeStart)/1e-6             ! Calculate time point for AppProfile and Speed filenames
	 Write(20,*) filetime, AppProfile(i)			! Write Appearance Profile

	 Write(speedFilename,'("Outputs/Speed/Speed_",I3.3,".txt")')filetime	! Create filenames for each timestep
         Open(unit=100+i,file=speedFilename)				
	 Do j = 1, numberOfSpeedPoints				! For each range of speed points write transverse speed
	  Write(100+i,*) (((j-1)*speedStep)+speedStart), TransSpeed(i,j)
	 End do
         Close(unit=100+i)

	End do
	Close(20)

! Write Array Count stats
	Open(unit=21, file='Outputs/Statistics.txt')
	Write(21,*) ncyc, "Total ingoing trajectories Created "
	Write(21,*) missWheel, "Trajectories missing the wheel "
	Write(21,*) 1.0*missWheel/ncyc, "Fraction of trajectories missing the wheel "
	Write(21,*) ncyc*npass, "Possible ingoing probed trajectories "
	Write(21,*) acceptedIngoing, "Ingoing trajectories probed "
	Write(21,*) 1.0*acceptedIngoing/(ncyc*npass), "Ingoing acceptance ratio "
	Write(21,*) (ncyc - missWheel)*npass, "Possible scattered trajectories "
	Write(21,*) acceptedScattered, "Scattered trajectories probed "
	Write(21,*) 1.0*acceptedScattered/((ncyc-missWheel)*npass), "Scattered trajectory acceptance ratio "

	Do k = 1, npass
	 Write(21,*) k, arrayCount(1,k), arrayCount(2,k)	! Write laser pass no, ingoing traj probed, scattered traj probed
	End do
	Close(21)


!***************************************
! PRINT INFORMATION TO SCREEN
!***************************************

! Calculate run time
    call cpu_time(endTime)

    runTime = endTime - startTime

! Calculate acceptance ratio
    totalTraj = ncyc*vectorsPerParticle
    acceptanceRatio = real(acceptedCounter)/((real(ncyc)*real(vectorsPerParticle)))


! Print info to screen

    print *,
    print *, "Finished in", runTime, "seconds"
    print *,
    print *, ncyc, "Total ingoing trajectories Created "
    print *, missWheel, "Trajectories missing the wheel "
    print '(a, F4.2, a)', "        ", 1.0*missWheel/ncyc, " Fraction of trajectories missing the wheel "
    print *, ncyc*npass, "Possible ingoing probed trajectories "
    print *, acceptedIngoing, "Ingoing trajectories probed "
    print '(a, F4.2, a)', "        ", 1.0*acceptedIngoing/(ncyc*npass), " Ingoing acceptance ratio "
    print *, (ncyc - missWheel)*npass, "Possible scattered trajectories "
    print *, acceptedScattered, "Scattered trajectories probed "
    print '(a, F4.2, a)', "        ", 1.0*acceptedScattered/((ncyc-missWheel)*npass), " Scattered trajectory acceptance ratio "
    print *,

end program MCScatteringFM
