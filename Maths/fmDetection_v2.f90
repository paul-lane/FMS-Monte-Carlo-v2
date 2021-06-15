Module fmDetection_v2
    use mathConstants
    
! Array shared by entire class - REMEMBER ALWAYS TO ALLOCATE BEFORE USE
    
    contains

! ****************************************************************************************************************************
! EULER SUBROUTINE
! Performs an Euler transformation on a vector using the YZY system
! i.e. first performs a rotation about the Y axis, then about the new Z axis and then nominally about the new Y axis
! this third rotation is not necessary so is hard coded to zero
! THis subroutine is used to rotate the co-ordinates or vectors of a trajectory from the lab frame into the frame
! of the laser beam
! ************************************ 
    
		Subroutine Euler(input, th1, th2, output)

			implicit none

! Define variables and constants

			Double Precision, intent(in) :: input(3) 				! vector be rotated
			Double Precision, intent(out) :: output(3)	   			! rotated vector
			Double Precision, intent(in) :: th1, th2 				! angles to rotate by
			Double Precision :: theta1, theta2, theta3
			Double Precision :: a11, a12, a13, a21, a22, a23, a31, a32, a33 	! matrix elements
			Double Precision :: c1, s1, c2, s2, c3, s3  				! sine and cosine of angles

! Convert angles to radians

			theta1 = (th1/180.0)*pi                            
			theta2 = (th2/180.0)*pi
			theta3 = 0.0
					
! Calculate sine and cosine

			c1 = cos(theta1)										
			s1 = sin(theta1)
			c2 = cos(theta2)
			s2 = sin(theta2)
			c3 = cos(theta3)
			s3 = sin(theta3)

! Matrix elements of Euler transform

			a11 = (c1*c2*c3)-(s1*s3) 							
			a12 = -(c1*s2)			 							
			a13 = (c3*s1)+(c1*c2*s3)
			a21 = c3*s2
			a22 = c2
			a23 = s2*s3
			a31 = -(c1*s3)-(c2*c3*s1)
			a32 = (s1*s2)
			a33 = (c1*c3)-(c2*s1*s3)


! Multiply Euler matrix by input vector

			output(1) = (a11*input(1))+(a12*input(2))+(a13*input(3)) 
			output(2) = (a21*input(1))+(a22*input(2))+(a23*input(3))
			output(3) = (a31*input(1))+(a32*input(2))+(a33*input(3))

        End Subroutine Euler 


!*********************************************************************************************************************
! REVERSE EULER TRANSFORM
! Performs a reverse Euler transform using the YZY system (the reverse of the EULER routine above) - for more details see EULER 
! This is based on the EULER subroutine but the matrix has been transposed to allow the reverse action to be performed
! This is used to convert parameters calculated in the frame of the laser to be converted back into the lab frame
!**************************************************

	Subroutine RevEuler(input, th1, th2, output)

                        implicit none

! Define variables and constants

                        Double Precision, intent(in) :: input(3)                          ! vector be rotated
                        Double Precision, intent(out) :: output(3)                        ! rotated vector
                        Double Precision, intent(in) :: th1, th2                          ! angles to rotate by
                        Double Precision :: theta1, theta2, theta3
                        Double Precision :: a11, a12, a13, a21, a22, a23, a31, a32, a33   ! matrix elements
                        Double Precision :: c1, s1, c2, s2, c3, s3                        ! sine and cosine of angles


! Convert angles to radians

                        theta1 = (th1/180.0)*pi
                        theta2 = (th2/180.0)*pi
                        theta3 = 0.0

! Calculate sine and cosine

                        c1 = cos(theta1)
                        s1 = sin(theta1)
                        c2 = cos(theta2)
                        s2 = sin(theta2)
                        c3 = cos(theta3)
                        s3 = sin(theta3)

! Matrix elements of Reversed Euler transform
! note the elements have been transposed what was a12 in the EULER is noW a21 and a13 is now a31 etc

                        a11 = (c1*c2*c3)-(s1*s3)
                        a21 = -(c1*s2)
                        a31 = (c3*s1)+(c1*c2*s3)
                        a12 = c3*s2
                        a22 = c2
                        a32 = s2*s3
                        a13 = -(c1*s3)-(c2*c3*s1)
                        a23 = (s1*s2)
                        a33 = (c1*c3)-(c2*s1*s3)

! Multiply Euler matrix by input vector

                        output(1) = (a11*input(1))+(a12*input(2))+(a13*input(3))
                        output(2) = (a21*input(1))+(a22*input(2))+(a23*input(3))
                        output(3) = (a31*input(1))+(a32*input(2))+(a33*input(3))

	End subroutine RevEuler

!***********************************************************************************************************************
! CHECK INTERSECT
! This is a subroutine to determine whether a line intersects a circle
! It is used to determine if the particle trajectory passes therough a laser beam
!******************************************

	Subroutine checkIntersect(Traj, startPos, rlaser, cellCentre, Intersect)        

			Implicit none

! Define variables and constants
        
			Double Precision, intent(in) :: Traj(3), startPos(3), cellCentre(3) 	! vector inputs of trajectory start pos and beam centre
			Double Precision, intent(in) :: rLaser 				! radius of laser beam(s)
			Double Precision :: a, b, c					! coefficients for quadratic equation
			logical, intent(out) :: Intersect				! True/False depending whether the particle vector intersects the beam				


! Use the quadratic equation to determine intersect
! a, b and c are the roots of the quadratic eqn

			a = (traj(2)**2) + (traj(3)**2)
			b = (2*traj(3)*startPos(3))-(2*traj(3)*cellCentre(3))+(2*traj(2)*startPos(2))-(2*traj(2)*cellcentre(2))
			c = (startPos(3)**2)-(2*cellCentre(3)*startPos(3))+(cellCentre(3)**2)+(startPos(2)**2)&
				-(2*cellCentre(2)*startPos(2))+(cellCentre(2)**2)-(rlaser**2)


! If this condition is met the trajectory intersects the laser beam

			if((b**2).gt. (4*a*c)) then		
				Intersect = .TRUE.
            		else
                		Intersect = .FALSE.
            		Endif

        End Subroutine checkIntersect
        
	
!*******************************************************************************************************************
! GET PROBE INTERSECTION
! This calculates the co-ordinates at which a line intersects a circle
! Used to calculate the position at which the trajectory enters and exits the laser beam
!**************************************

	Subroutine getProbeIntersection(Traj, startPos, rLaser, cellCentre, Pos)        

           	 Implicit none

			Double Precision, intent(in) :: Traj(3), startPos(3), cellCentre(3) 	! vector inputs of traj, start pos at cell pos
			Double Precision, intent(in) :: rLaser 					! radius of laser beam(s) 
			Double Precision, intent(out) :: Pos(2,3)				! Two vector intersect positions
			Double Precision :: a, b, c						! coefficients for quadratic equation			
			Double Precision :: root1, root2					! two solutions to quadratic

! Use the quadratic equation to determine intersect
! a, b and c are the roots of the quadratic eqn

                        a = (traj(2)**2) + (traj(3)**2)
                        b = (2*traj(3)*startPos(3))-(2*traj(3)*cellCentre(3))+(2*traj(2)*startPos(2))-(2*traj(2)*cellcentre(2))
                        c = (startPos(3)**2)-(2*cellCentre(3)*startPos(3))+(cellCentre(3)**2)+(startPos(2)**2)&
                                -(2*cellCentre(2)*startPos(2))+(cellCentre(2)**2)-(rlaser**2)

! Calculate the solutions to the quadratic equations

			root1 = (-b + sqrt((b**2)-(4*a*c)))/(2*a)
			root2 = (-b - sqrt((b**2)-(4*a*c)))/(2*a)
			
        
! Use the parameterized equations to convert into  x, y, z positions for each root

			Pos(1,1) = (traj(1)*root1)+startPos(1)
			Pos(2,1) = (traj(1)*root2)+startPos(1)
			Pos(1,2) = (traj(2)*root1)+startPos(2)
			Pos(2,2) = (traj(2)*root2)+startPos(2)
			Pos(1,3) = (traj(3)*root1)+startPos(3)
			Pos(2,3) = (traj(3)*root2)+startPos(3)

        End subroutine getProbeIntersection
            
!*****************************************************************************************************************
! HIT WHEEL
! Uses chamber geometry and position of impact at the wheel position to determine 
! if the particle hit the wheel or not
!******************************************************

	Subroutine hitWheel(particleStartPos, wheelCentrex, wheelCentrey, wheelRadius, bathTop, hit)
	
	implicit none
	
	Double Precision, dimension(3), intent(in):: particleStartPos
	Double Precision, intent(in):: wheelCentrex, wheelCentrey, wheelRadius, bathTop
	Double Precision:: radialDist
	Logical, intent(out):: hit

! Calculate the radial distance of the trajectory from the centre of the wheel

	radialDist = sqrt(((particleStartPos(1)-wheelCentrex)**2) + ((particleStartPos(2)-wheelCentrey)**2))
	
! The only particles that hit the wheel are those which have a radial distance less than or equal to the radius of the wheel AND
! are above the top of the bath (particles within the radius of the wheel but below the bath top hit the bath and not the wheel)

	if((radialDist .le. wheelRadius) .and. (particleStartPos(2) .gt. bathTop)) then
		hit = .true.
	else
		hit = .false.		
	endif

	End subroutine hitWheel
        
!*****************************************************************************************************************
! GET INTERSECTION TIME
! Uses particle speed, creation time and intersection co-ordinates to calculate the time at which a particle
! enters and exits the laser beam
!******************************************************
        
        Subroutine getIntersectionTime(Pos, particleStartPos, particleVector, particleSpeed, particleTime, entryTime, exitTime)
            implicit none

! Define variables and constants
            Double Precision, dimension(2,3), intent(in) :: Pos					! Pos are co-ords for entry and exit of laser beam
            Double Precision, dimension(3), intent(in) :: particleStartPos, particleVector	! Particle vector and start pos
            Double Precision, intent(in) :: particleSpeed, particleTime				! Particle speed and creation time
            Double Precision, dimension(2) :: intersectionTime					! Intersection times
            Double Precision, intent(out) :: entryTime, exitTime					! Entry and exit times
	    DOuble Precision :: magStart, magVector						! Magnitudes of start and vector
            integer :: i

! Initialize to zero

            intersectionTime = 0

! Looks at both intersection points and calculates the time

            do i = 1, 2

! Calculate intersection times

		 intersectionTime(i) = particleTIme + (abs(particleStartPos(3) - Pos(i,3))/abs(particleVector(3)*particleSpeed))
           End do
            
! Sorts entry and exit times

            entryTime = minval(intersectionTime, mask = intersectionTime .gt. 0)
            exitTime = maxval(intersectionTime)

            End Subroutine getIntersectionTime

 
!********************************************************************************************************************
! START END TIME POINTS
! Uses the entry and exit time to direct output to the correct time bin
!************************************************* 

           
	subroutine startEndTimePoints(NumberOfTimepoints, entryTime, exitTime, probeStart, probeEnd, tStep, startTimePoint, endTimePoint)
            ! Uses the entry time and exit time to direct the output to the correct time
	    implicit none

! Define variables and constants

            integer, intent(in) :: NumberOfTimePoints
            Double Precision, intent(in) :: entryTime, exitTime, probeStart, probeEnd, tStep
            integer, intent(out) :: startTimePoint, endTimePoint

            
! If entry time is less than probe start time, then that particle is only recorded from the first time bin
! (if input parameters are setup correctly this should not be the case...)

            if (entryTime .lt. probeStart) then
                startTimePoint = 1

! Otherwise finds the time corresponding to the start of the bin that the particle finds itself in
! (if input parameters are setup correctly this should always be the case...)
            else

		startTimePoint = floor((entryTime - probeStart) /tStep) + 1

            end if

! If exit time is greater than end of probe, then end timepoint then particle is imaged all the way til the end of the probe time
! (Again if input parameters are setup correctly this shouldn't be the case)

            if (exitTime .gt. probeEnd) then
                endTimePoint = NumberOfTimepoints

! Otherwise finds the bin corresponding to the last time point that the particle is in the laser beam
            else

                endTimePoint = floor((exitTime - probeStart) / tStep) + 1

            end if

        end subroutine startEndTimePoints            
            
        
                
        
!*************************************************************************************************************************
! PROBE
! FM probe routine calculates the apperance profile and the transverse velocities
!**************************************************
            
        subroutine Probe(AppProfile, TransSpeed, NumberOfTimePoints, entryTime, exitTime, startTimePoint, endTimePoint, probeStart, tStep, particleSpeed, particleVector, particleStartPos, startSpeed, speedStep)

            implicit none

! Define variables and constants
            Double Precision, intent(inout), dimension(:) :: AppProfile
            Double Precision, intent(inout), dimension(:,:) :: TransSpeed
            integer, intent(in) :: NumberOfTimePoints, startTimePoint, endTimePoint
            Double Precision, intent(in) :: probeStart, tStep, particleSpeed, startSpeed, speedStep !, t0
            Double Precision, dimension(3), intent(in) :: particleVector, particleStartPos
            Double Precision, intent(in) :: entryTime, exitTime
	    integer :: t, s
            Double Precision :: currentTime, upLimit, lowLimit, weight, xSpeed


! Loops from entry timepoint to exit timepoint to avoid wasting cycles when particle is not within beam

            do t = startTimePoint, endTimePoint
                
! Iterate time from probe start time in steps of defined size
! currentTime refers to the start of this time window

                currentTime = probeStart + (t-1)*tStep 

! The weighting of a point is based on how long it spends in the beam in a particular time window
! To do this we use the entry and exit times of the vector and the start and end of each window
! We use lowlimit and uplimit for each time step if the particle spends the whole of the timestep in the
! the beam the lowlimit is the start of this time window and the upper limit the end of this window.
! If the particle enters the beam midway through the time window but stays until the end then the 
! low limit is its entry time and uplimit the end of the window.
! If the particle was already in the window at the start but leaves before the end it's low limit will be
! the start of the window and uplimit it's exit time.
! If the particle enters and exits in the same time window then the lowlimit will be it's entry time and 
! uplimit it's exit time.
! This ensures that the correct weighting is given to each particle in each time window.

                if(entryTime .GE. currentTime) then
                    lowLimit = entryTime

                else
                    lowLimit = currentTime

                endif
                
! currentTime + tStep is the end time of this window (also start of next)

                if(exitTime .GE. (currentTime + tStep)) then
                    upLimit = currentTime + tStep

                else
                    upLimit = exitTime

                endif
                
! weight in this window is therefore the difference between the upper and lower limits

                weight = upLimit - lowLimit
!                weight = 1.0				! Use if counting particle numbers rather than time in beam

! Add appropriate weight to array

                AppProfile(t) = AppProfile(t) + weight
                
! Transverse speed is x component of vector multiplied by speed

		xSpeed = abs(particleVector(1)*particleSpeed)		

! Transverse velocity bin is found 

                s = floor((xSpeed - startSpeed) / speedStep) + 1
                
! Transverse velocity array is updated

                TransSpeed(t,s) = TransSpeed(t,s) + weight
		
            End do
            
            End subroutine Probe


!************************************************************
! GET BEAM PROFILE AT GIVEN Z POSITION
!************************************************************

        Subroutine getProfile(particleVector, particleStartPos, z, profile)

                implicit none
                Double Precision, dimension(3), intent(in) :: particleVector, particleStartPos
                Double Precision, intent(in) :: z
                Double Precision, dimension(2), intent(out) :: profile

                profile(1) = (particleVector(1)*z) + particleStartPos(1)
                profile(2) = (particleVector(2)*z) + particleStartPos(2)
        End Subroutine



!********************************************************************************************************************
! SPEED TO ABSORPTION LINE
!*********************************************************


	    Subroutine speedToAbsorptionLine(TransSpeed, nSpeeds, nPad, nAbs, dSpeed, integral)

	    Implicit none
	
	    Double Precision, dimension(:), intent(in) :: TransSpeed
	    Double Precision, allocatable, intent(inout) :: integral(:)
	    Double Precision, intent(in) :: dSpeed
	    Double Precision :: sum 
	    Double Precision, allocatable :: array1(:), array2(:)
	    integer, intent(in) :: nSpeeds, nPad, nAbs
	    integer ::  i, n, centre 

	    centre = nAbs / 2						! Line centre is at midpoint of array
	    
	    allocate(array1(centre-nPad), array2(centre-nPad))

	    integral = 0.0						! Initialize array values
	    sum = 0.0
	    array1 = 0.0
 	    array2 = 0.0 	


	    Do i = 1, (centre - nPad)		
	     array1(i) =  (TransSpeed(i)+TransSpeed(i-1)) * (dspeed/2)
	    End do	

	    Do i = 1, (centre - nPad)
		n = (centre-nPad)-i+1
		array2(i) = array1(n) + sum
		sum = sum + array1(n)
	    End do


	    Do i = 1, nAbs
		if(i .LE. nPad) then
		 integral(i) = 0.0
		elseif((i .GT. nPad) .AND. (i .LE. centre)) then
		 integral(i) = array2(i-npad)
		else
		 integral(i) = integral(centre-(i-centre-1)) 		!!ADDED MINUS ONE
		endif                         
	    End do	         


	    End Subroutine speedToAbsorptionLine

!************************************************************************
! INTEGRAL ABSORPTION TO FM ABSORPRION
!***************************************

	    Subroutine intAbsToFM(FMabs, nFMpoints, modStep, integralAbs)
	    
	    Implicit none

	    Double Precision :: J0, J1, J2, J3
	    Double Precision :: J0J1, J1J2, J2J3 
	    Double Precision, intent(in) :: modStep
	    integer, intent(in):: nFMpoints
	    Double Precision, dimension(:), intent(inout) :: FMabs
	    Double Precision, dimension(:), intent(in) :: integralAbs
	    integer :: iDm1, iDp1, iDm2, iDp2, iDm3, iDp3, i
	    Double Precision :: Order1, Order2, Order3

	    J0 = 0.870363		! Bessel functions for modulation index of 0.9
	    J1 = 0.276393
	    J2 = 0.0509452
	    J3 = 0.00663612

	    J0J1 = J0 * J1
	    J1J2 = J1 * J2
	    J2J3 = J2 * J3

	    Do i = 1, nFMpoints
		iDm1 = ((3*modStep) + i) - modStep		! index of delta -1
		iDp1 = ((3*modStep) + i) + modStep		! index of delta +1
		iDm2 = ((3*modStep) + i) - (2 *	modStep)	! index of delta -2
		iDp2 = ((3*modStep) + i) + (2 * modStep)	! index of delta +2
		iDm3 = ((3*modStep) + i) - (3 * modStep)	! index of delta -3
 		iDp3 = ((3*modStep) + i) + (3 * modStep)	! index of delta +3
	    
		Order1 = (J0J1+J1J2) * (integralAbs(iDm1)-integralAbs(iDp1))
		Order2 = (J1J2+J2J3) * (integralAbs(iDm2)-integralAbs(iDp2))
		Order3 = J2J3 * (integralAbs(iDm3) - integralAbs(iDp3))
	
		FMabs(i) = Order1 + Order2 + Order3 

	    End do

	    End Subroutine intAbstoFM



	End Module fmDetection_v2         
