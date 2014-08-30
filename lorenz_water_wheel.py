#/usr/bin/env python
"""
Author: Ambrose Granizo-Mackenzie
Breif: Animates and saves data for a Newtonian model of a Lorenz water wheel.
Notes: Highly based on the java code by Christopher Wellons (mosquitopsu@gmail.com)
        from http://git.nullprogram.com/ChaosWheel.git
"""


import os, pygame, pygame.locals
import random, math



def main(numBuckets):
    # Initialize pygame elements
    pygame.init()
    screen = pygame.display.set_mode((SIZE, SIZE))
    pygame.display.set_caption('Lorenz Water Wheel Simulation')
    pygame.mouse.set_visible(1)

    # Create The Backgound
    background = pygame.Surface(screen.get_size())
    background = background.convert()
    background.fill((250, 250, 250))
    
    # Display The Background
    screen.blit(background, (0, 0))
    pygame.display.flip()

    #Prepare Game Objects
    clock = pygame.time.Clock()

    # Math
    
    '''
    Chaotic 
    theta =  3.36172078085
    thetadot = -0.30030302512
    '''
    
    thetaI = 2.0 * math.pi * random.random()
    theta = thetaI
    thetadotI = random.random() - 0.5
    thetadot = thetadotI
    buckets = []
    for i in range(numBuckets):
        buckets.append(0.0)
        
    # main Loop
    counter = 0
    while counter <= MAX_COUNT:
        clock.tick(DELAY)
        # Handle input events
        for event in pygame.event.get():
            if event.type == pygame.locals.QUIT:
                counter = MAX_COUNT
            elif event.type == pygame.locals.KEYDOWN and event.key == pygame.locals.K_ESCAPE:
                counter = MAX_COUNT
            #elif event.type == pygame.locals.MOUSEBUTTONDOWN:
                #print "MBD"
            #elif event.type is pygame.locals.MOUSEBUTTONUP:
                #print "MDU"
        # Draw buckets
        diff = math.pi * 2.0 / len(buckets)
        size = min(screen.get_size())
        bucketSize = size / int((len(buckets) / 1.25))
        drawRadius = int(size / 2 - bucketSize)
        centerx = int(size/2)
        centery = int(size/2)
        for i in range(len(buckets)):
            angle = i * diff + theta - math.pi / 2
            x = centerx + int(math.cos(angle)*drawRadius)
            y = centery + int(math.sin(angle)*drawRadius)
            bucket = pygame.Rect(x - bucketSize / 2, y - bucketSize /2, bucketSize, bucketSize)
            pygame.draw.rect(background, (0,0,0), bucket, 2)
            height = int(bucketSize * buckets[i] / bucketFull)
            fluid = pygame.Rect(x - bucketSize / 2, y - bucketSize / 2 + (bucketSize - height), bucketSize, height)
            pygame.draw.rect(background, (0,0,255), fluid)
        
        # Refresh background
        screen.blit(background, (0, 0))
        background.fill((250, 250, 250))
        pygame.display.flip()
        theta, thetadot, buckets = updateState(DELAY / 1500.0, theta, thetadot, buckets)
        logState(theta, thetadot, buckets)

        counter+=1
    i = 0

    pygame.quit()
    
    # param = '_numBuckets-'+str(numBuckets)+'_radius-'+str(radius)+'_inertia-'+str(wheelInertia)+'_damping-'+str(damping)+'_gravity-'+str(gravity)+'_bucketFull'+str(bucketFull)+'_drainRate-'+str(drainRate)+'_fillRate-'+str(fillRate)+'_theteI-'+str(thetaI)+'_thetedotI-'+str(thetadotI)
    #param = '2'
    #'water_wheel'+param+'.txt'
    
    f = open(OUTPUTFILE, 'w')
    for i in range(len(rlRatio)):
        f.write(str(thetadots[i]) +'\t'+ str(rlRatio[i]) +'\t' + str(tbRatio[i])+'\n')
    f.close()
    
    
        
def updateState(tdot, theta, thetadot, buckets):
    
    # Store the original state
    thetaOrig = theta
    thetadotOrig = thetadot
    bucketsOrig = buckets
    
    # Variables for RK4
    dt = 0.0
    rateWeight = 1.0
    
    # Time deriatives of state
    ddtTheta = 0.0
    ddtThetadot = 0.0
    ddtBuckets = map(lambda x: 0.0, buckets)
    
    # Total derivatives approximations
    ddtThetaTotal = 0.0
    ddtThetadotTotal = 0.0
    ddtBucketsTotal = map(lambda x: 0, buckets)
    
    # 4th order Runge–Kutta Integration
    for r4idx in range(1,5):
        if (r4idx > 1) and (r4idx < 4):
            rateWeight = 2.0
            dt = tdot / 2.0
        elif (r4idx == 4):
            rateWeight = 1.0
            dt = tdot
        # System states to be used in RK4 step
        theta = thetaOrig + dt * dt * ddtTheta
        while theta < 0:
            theta += math.pi * 2.0
        while theta > math.pi * 2.0:
            theta-=math.pi * 2.0
        thetadot = thetadotOrig + dt * ddtThetadot
        
        for i in range(len(buckets)):
            val = float(bucketsOrig[i]) + dt * ddtBuckets[i]
            buckets[i] = min(bucketFull, max(0, val))
        
        # Differential equation for ddTheta (kinematics)
        ddtTheta = thetadot
        
        # Calculate inertia
        inertia = wheelInertia
        for i in range(len(buckets)):
            inertia += buckets[i] * radius * radius
        
        # Calculate torque
        torque = -1 * (damping * thetadot)
        diff = math.pi * 2.0 / float(len(buckets))
        for i in range(len(buckets)):
            torque += buckets[i] * radius * gravity * math.sin(theta + diff * i)
        
        # Differential equation for ddtThetadot (physics)
        ddtThetadot = torque / inertia
        
        # Differential equation for ddtBuckets (drain rate equation)
        for i in range(len(buckets)):
            ddtBuckets[i] = buckets[i] * -drainRate + inflow(theta + diff *i, buckets)
        
        # Log derivative approximations
        ddtThetaTotal += ddtTheta * rateWeight
        ddtThetadotTotal += ddtThetadot * rateWeight
        for i in range(len(ddtBucketsTotal)):
            ddtBuckets[i] = float(ddtBucketsTotal[i]) + ddtBuckets[i] * float(rateWeight)
        
        # End of RK4 loop
    
    theta = thetaOrig + 1.0 / 6.0 * ddtThetaTotal * tdot
    thetadot = thetadotOrig + 1.0 / 6.0 * ddtThetadotTotal * tdot
    for i in range(len(ddtBucketsTotal)):
        val = float(buckets[i]) + 1.0 / 6.0 * float(ddtBucketsTotal[i]) * float(tdot)
        buckets[i] = min(bucketFull, max(0, val))
        
    return (theta, thetadot, buckets)
            
def inflow(angle, buckets):
    lim = abs(math.cos(math.pi * 2.0 / len(buckets)))
    if (math.cos(angle) > lim): 
        return fillRate / 2.0 * (math.cos(len(buckets) * math.atan2(math.tan(angle), 1)) / 2.0 + 1)
    else:
        return 0
    
def logState(theta, thetadot, buckets):
    left = 0
    right = 0
    top = 0
    bottom = 0
    diff = math.pi * 2.0 / len(buckets)
    for i in range(len(buckets)):
        angle = theta + diff * i
        if math.cos(angle) > 0: right += buckets[i]
        left +=buckets[i]
        if math.sin(angle) > 0: top += buckets[i]
        bottom += buckets[i]
    rl = left / right
    tb = top / bottom 
    rlRatio.append(rl)
    tbRatio.append(tb)
    thetadots.append(thetadot)
    return (rlRatio, tbRatio)


if __name__ == '__main__':
    # Output file path
    OUTPUTFILE = "water_wheel.txt"
    
    # Simulation constants.
    SIZE = 300 # display size in pixels
    DELAY = 30 # milliseconds
    DEFAULT_BUCKETS = 10
    MIN_BUCKETS = 5

    # Time of the simulation
    TIME = 1200
    MAX_COUNT = 30 * TIME / DELAY

    # Simulation parameters.
    radius = 1 # meters
    wheelInertia = 0.1 # kg * m^2
    damping = 3.5 # m * kg / radians / sec
    gravity = 32.2 # m / s ^ 2
    bucketFull = 1.0 # kg
    drainRate = .2 # slug / sec / slug
    fillRate = .33 # kg / sec    
    
    thetadots = []
    rlRatio = []
    tbRatio = []
    
    main(DEFAULT_BUCKETS)
