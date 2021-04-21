module speeds_v2
    use mathConstants

    contains

        ! Calculates speed of ingoing particle based on cumulative integral function
        subroutine ingoing_speed(x0, aMax, aMin, h, s, dist, pulseLength, speed, t0)

            ! variables relating to cumulative integral function of arrival times
            double precision, intent(in) :: x0, aMax, aMin, h, s
            double precision, intent(in) :: dist, pulseLength
            double precision, intent(out) :: speed, t0
            double precision :: t, x, arrivalTime

            ! Calculate random time of creation
            call random_number(t)
            t0 = (t*pulseLength) - (pulseLength/2.0)

            ! CaLculate TOF based on cumulative integral function from real data anf fit by Origin.
            ! Function in Origin is called Logistics5.
            call random_number(x)
            arrivalTime = x0/(((aMax-aMin)/(x-aMin))**(1.0D0/s)-1.0D0)**(1.0D0/h)
            speed = dist/(arrivalTime*1D-6)
        end subroutine ingoing_speed


	! Calculates speed of an ingoing particle based on cumulative integral of 3 Gaussians
	subroutine ingoing_speed_from_Gauss(w1, w2, w3, m1, m2, m3, std1, std2, std3, dist, pulseLength, speed, t0)

	    double precision, intent(in) :: w1, w2, w3			
	    double precision, intent(in) :: dist, pulseLength
	    double precision, intent(out) :: speed, t0
	    double precision :: t, x 	    
	    real*4 :: m1, m2, m3, std1, std2, std3
	    real*4 :: arrivalTime

	    ! Calculate random time of creation
	    call random_number(t)
	    t0 = (t*pulseLength) - (pulseLength/2.0)

	    ! Calculate TOF based on cumulative integral function from real data fitted to 3 Gaussians
	    call random_number(x)

	    if(x .le. w1) then
		call rnm(m1, std1, arrivalTime)
		
	    else if(x .le. (w1+w2)) then
		call rnm(m2, std2, arrivalTime)

	    else
		call rnm(m3, std3, arrivalTime)
	    end if

	    speed = dist/(arrivalTime*1.0D-6)
	
	end subroutine ingoing_speed_from_Gauss

	! subroutine to draw a random number from a Gaussian distribution
	subroutine rnm(rmean, sd,rand)

	double precision random
      	real*4 x,u,s,t,c0,c1,c2,d1,d2,d3,t2,t3,rand,rmean,sd

	call random_number(x)

	if(x.gt.0.999999)then
       		x=0.999999
      	endif

      	if(x.lt.1.0e-6)then
       		u=1.0e-06
      	else
       		u=x
      	endif

	s = 1.

      	if (u .gt. 0.5) then
      		u = 1. - u
      		s = -1.
      	end if

      	t = sqrt(- (2 * log(u)))
      	c0 = 2.515517
      	c1 = 0.802853
      	c2 = 0.010328
      	d1 = 1.432788
      	d2 = 0.189269
      	d3 = 0.001308
      	t2 = t * t
      	t3 = t * t2
      	x = s * (t - (((c0 + (c1 * t)) + (c2 * t2)) / (((1. + (d1 * t)) + (d2 * t2)) + (d3 * t3))))

      	rand = rmean + (x * sd)
      	return
      	end subroutine rnm




        ! Calculates speed based on Maxwell-Boltzmann Distribution of speeds
        subroutine MB_speed(maxSpeed, temp, mass, mostLikelyProbability, scatteredSpeed)

            logical :: hit
            double precision, intent(in) :: maxSpeed, temp, mass, mostLikelyProbability
            double precision, intent(inout) :: scatteredSpeed
            double precision :: rand1, rand2, probability, normalisedProbability

            hit = .FALSE.

            do while (hit .eqv. .FALSE.)
                call random_number(rand1)
                scatteredSpeed = rand1*maxSpeed

                probability = MB_probability(temp, scatteredSpeed, mass)

                ! Calculates the probability of the speed with respect to the most probable speed equalling 1.
                ! The Maxwell-Boltzmann distribution is already normalised to 1, meaning that the sum of all
                ! probabilities from zero to infinity will equal 1.
                ! It is possible to avoid this step, however, it would take a very long time to
                ! achieve a hit due to the small value of probability.
                normalisedProbability = probability/mostLikelyProbability

                call random_number(rand2)

                if (normalisedProbability .gt. rand2) then
                    hit = .TRUE.
                end if
            end do
        end subroutine MB_speed

        ! Generates random Lorentizan distributed value. Gamma is the scale parameter equal to HWHM of desired function.
        subroutine lorentzian_distribution(gamma, speed)
            implicit none

            double precision :: speed, rand
            double precision, intent(in) :: gamma

            call random_number(rand)

            speed = gamma*tan(pi*(rand-0.5D0))

        end subroutine lorentzian_distribution

        ! This method apparently is very efficient, however it generates two Gaussian distributed numbers at a time
        ! Make sure this is incorporated somehow to avoid wasting cycles
        subroutine gaussian_distribution(mean, sigma, z1, z2)
            implicit none

            double precision :: rand1, rand2, v1, v2, rSquared, z1, z2, mean, sigma

            do
                call random_number(rand1)
                call random_number(rand2)

                v1 = (2.0*rand1) - 1.0
                v2 = (2.0*rand2) - 1.0

                rSquared = (v1**2.0) + (v2**2.0)

                if (rSquared .lt. 1) then
                    z1 = v1*SQRT((-2.0*log(rSquared))/rSquared)
                    z2 = v2*SQRT((-2.0*log(rSquared))/rSquared)

                    z1 = mean + sigma*z1
                    z2 = mean + sigma*z2

                    EXIT
                end if
            end do
        end subroutine

        ! Adds a transverse speed to molecule exiting final apperture 
        subroutine transverse_temp_test(mean, sigma, gamma, gaussLorFraction, zPos, travelDistance, &
             startTime, speed, startPoint, vector) 
            implicit none 
 
            double precision, dimension(3) :: startPoint 
            double precision, dimension(3), intent(inout) :: vector 
            double precision, intent(in) :: mean, sigma, gaussLorFraction, gamma 
            double precision :: zPos, travelDistance, startTime, speed, GSpeed1, GSpeed2, LSpeed 
 
            startTime = (startTime + abs(travelDistance/(vector(3)*speed))) 
            startPoint(1) = startPoint(1) + (startTime*speed*vector(1)) 
            startPoint(2) = startPoint(2) + (startTime*speed*vector(2)) 
            startPoint(3) = zPos 
     
            vector = vector*speed 
 
            ! Generates Lorentzian speed and TWO Gaussian speeds 
            call lorentzian_distribution(gamma, LSpeed) 
            call gaussian_distribution(mean, sigma, GSpeed1, GSpeed2) 
 
            ! Modifies x-direction of vector with new speed 
            ! Gaussian and Lorentzian fraction weights contribution of each speed 
            vector(1) = (gaussLorFraction*GSpeed1) + ((1D0 - gaussLorFraction)* LSpeed) 
 
            call lorentzian_distribution(gamma, LSpeed) 
 
            ! Same as for x-direction, this modifies y-direction 
            vector(2) = (gaussLorFraction*GSpeed2) + ((1D0 - gaussLorFraction)* LSpeed) 
 
            ! Normalises vector so that molecule is still travelling at original speed overall 
            vector = vector/norm2(vector) 
 
        end subroutine transverse_temp_test 

        subroutine transverse_temp(mean, sigma, gamma, l_g_fraction, zPos, travelDistance, startTime, speed, startPoint, vector)
            implicit none

            double precision, dimension(3) :: startPoint
            double precision, dimension(3), intent(inout) :: vector
            double precision, intent(in) :: mean, sigma, l_g_fraction, gamma
            double precision :: zPos, travelDistance, startTime, speed, transSpeed, rand, z2

            startTime = (startTime + abs(travelDistance/(vector(3)*speed)))
            startPoint(1) = startPoint(1) + (startTime*speed*vector(1))
            startPoint(2) = startPoint(2) + (startTime*speed*vector(2))
            startPoint(3) = zPos

            vector = vector*speed

            call random_number(rand)

            if (rand .gt. l_g_fraction) then
                call lorentzian_distribution(gamma, transSpeed)
            else
                call gaussian_distribution(mean, sigma, transSpeed, z2)
            end if

            vector(1) = transSpeed

            call random_number(rand)

            if (rand .gt. l_g_fraction) then
                call lorentzian_distribution(gamma, transSpeed)
            else
                call gaussian_distribution(mean, sigma, transSpeed, z2)
            end if

            vector(2) = transSpeed

            vector = vector/norm2(vector)

        end subroutine transverse_temp

        ! Finds probability of particle travelling at given speed
        function MB_probability(temp, speed, mass) result(probability)

            double precision :: part1, part2, part3, speed, temp, mass, probability

            !part 1, 2, 3 correspond to individual parts of the maxwell-boltzmann distribution
            ! formula for calculating probability of a given speed
            part1 = 4.0D0*pi*speed*speed
            part2 = (mass/(2*pi*boltzmannConstant*temp))**(3.0D0/2.0D0)
            part3 = DEXP((-mass*speed*speed)/(2.0D0*boltzmannConstant*temp))

            probability = part1*part2*part3

        end function MB_probability

        ! Finds the most probable speed and its probability to use in normalisation 
        function MB_most_likely (temp, mass) result(mostLikelyProbability)

            double precision :: temp, mass, mostProbableSpeed, mostLikelyProbability

            mostProbableSpeed = sqrt((2.0D0*boltzmannConstant*temp)/mass)
            mostLikelyProbability = MB_probability(temp, mostProbableSpeed, mass)

        end function MB_most_likely

        subroutine deflection_angle(ingoing, outgoing, deflectionAngle)
            implicit none

            double precision, intent(in), dimension(3) :: ingoing, outgoing
            double precision, intent(out) :: deflectionAngle

            ! since this dot product finds the angle between the two vectors, it necessarily finds the deflection angle
            ! this is because the vectors are assumed to begin at the same point, and this is not the case with
            ! the ingoing and outgoing vectors, so the step where the angle is subtracted from 180 is not necessary
            deflectionAngle = acos(dot_product(ingoing,outgoing) / (norm2(ingoing)*norm2(outgoing))) * (360.0D0/(2*pi))

        end subroutine deflection_angle

        subroutine soft_sphere_speed(mass, internalRatio, surfaceMass, initialSpeed, deflectionAngle, finalSpeed)
            implicit none

            double precision :: initialEnergy, finalEnergy, massRatio, surfaceMass &
            ,part1, part2, part3, part4, part5, internalRatio, energyDiff, mass
            double precision, intent(in) :: initialSpeed, deflectionAngle
            double precision, intent(out) :: finalSpeed
            
            massRatio = mass/surfaceMass*1000.0D0
            initialEnergy = 0.5D0 * mass * initialSpeed * initialSpeed

            part1 = (2.0D0*massRatio)/((1+massRatio)**2.0D0)

            part2 = 1 + (massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2.0))
        
            part3 = cos(deflectionAngle*(2*pi/360.0D0))

            part4 = SQRT(1 - (massRatio*massRatio*(sin(deflectionAngle*((2*pi)/360.0D0))**2)) - internalRatio*(massRatio + 1))
        
            part5 = internalRatio*((massRatio + 1.0)/(2.0*massRatio))

            energyDiff = part1*(part2 - (part3 * part4) + part5) * initialEnergy

            finalEnergy = initialEnergy - energyDiff

            finalSpeed = SQRT(2*finalEnergy/mass)

        end subroutine soft_sphere_speed

end module speeds_v2
