import numpy as np
import matplotlib.pyplot as plt

####################
#
# Project 2
# PHYS 3926
#
# Author: Aidan Sirbu
#
####################



####################
# projectile() function models the projective motion of a baseball
#   after being struck with the option of accounting for air resistance
#
# Inputs:
# vel - int/float velocity in m/s
# phi - int/float launch angle in degrees
# tau - int/float time step in seconds
# air - 0 for no air resistance; 1 for air resistance
# method - e for Euler's method; ec for Euler-Cromer method; m for midpoint method;
#   all returns a plot of all methods using initial conditions, a range will not be returned
# plot - 0 for no plotting; 1 for plotting
#
# Output:
# - distance - float horizontal range of ball in meters
# - Plot of projectile motion
####################
def projectile(vel, phi, tau, air, method, plot):
    redundant = 1
    # CONSTANTS
    G = 9.81
    Cd = 0.35
    RHO = 1.2
    D = 0.074
    A = np.pi*(D/2)**2
    M = 0.145

    # Initial Conditions
    rx = 0
    ry = 1
    r = [rx, ry]
    phi = phi*np.pi/180
    vx = vel*np.cos(phi)
    vy = vel*np.cos(phi)
    v = [vx, vy]
    RX = 0
    RY = 1
    VX = vel*np.cos(phi)
    VY = vel*np.cos(phi)

    # Value Memory
    velList = [v]
    posList = [r]

    # Value Memory in the case of all methods
    velListE = [v]
    posListE = [r]
    velListEC = [v]
    posListEC = [r]
    velListM = [v]
    posListM = [r]

    # First alpha instance
    if air == 1:
        alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
    else:
        alpha = 0

    ##
    # Euler loop
    ##
    if method == "e":
        while redundant == 1:
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

    ##
    # Euler-Cromer loop
    ##
    elif method == "ec":
        while redundant == 1:
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

    ##
    # Midpoint loop
    ##
    elif method == "m":
        while redundant == 1:
            # Used for averaging velocities
            vxPrev = vx
            vyPrev = vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*(vx+vxPrev)/2
            ry = ry+tau*(vy+vyPrev)/2
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

    ##
    # All method loop
    ##
    if method == "all":
        ##
        # Euler loop
        ##
        while redundant == 1:
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velListE.append(v) # Update memory
            posListE.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

        # Reinitialize initial conditions
        vx = VX
        vy = VY
        rx = RX
        ry = RY
        v = [vx, vy]
        alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)

        ##
        # Euler-Cromer loop
        ##
        while redundant == 1:
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velListEC.append(v) # Update memory
            posListEC.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

        # Reinitialize initial conditions
        vx = VX
        vy = VY
        rx = RX
        ry = RY
        v = [vx, vy]
        alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)

        ##
        # Midpoint loop
        ##
        while redundant == 1:
            # Used for averaging velocities
            vxPrev = vx
            vyPrev = vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*(vx+vxPrev)/2
            ry = ry+tau*(vy+vyPrev)/2
            # Termination criterion
            if ry <= 0:
                break
            r = [rx, ry] # Update position vector
            velListM.append(v) # Update memory
            posListM.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0


    # Finding correct range
    if method != "all": # If all methods are used, multiple ranges will be reached, only a plot is needed
        lastIndex = len(posList) - 1
        xi = posList[lastIndex - 1]
        xf = posList[lastIndex]
        slope = (xf[1]-xi[1])/(xf[0]-xi[0])
        distance = xf[0]-xf[1]/slope

    # Plotting
    if method != "all":
        # Declare method used in plot title
        if method == "e":
            name = "Euler's"
        if method == "ec":
            name = "Euler-Cromer"
        if method == "m":
            name = "Midpoint"
        if plot == 1:
            xcomp = []
            ycomp = []
            i = 0
            while i < len(posList):
                xcomp.append(posList[i][0])
                ycomp.append(posList[i][1])
                i += 1
            plt.plot(xcomp, ycomp)
            plt.xlabel("Horizontal distance (m)")
            plt.ylabel("Vertical distance (m)")
            plt.title("Projectile motion of baseball using " + name + " method")
            plt.show()
    else:
        xcompE = []
        ycompE = []
        xcompEC = []
        ycompEC = []
        xcompM = []
        ycompM = []
        i = 0
        j = 0
        k = 0
        while i < len(posListE):
            xcompE.append(posListE[i][0])
            ycompE.append(posListE[i][1])
            i += 1
        while j < len(posListEC):
            xcompEC.append(posListEC[j][0])
            ycompEC.append(posListEC[j][1])
            j += 1
        while k < len(posListM):
            xcompM.append(posListM[k][0])
            ycompM.append(posListM[k][1])
            k += 1
        plt.plot(xcompE, ycompE, color="green")
        plt.plot(xcompEC, ycompEC, color="blue")
        plt.plot(xcompM, ycompM, color="red")
        plt.xlabel("Horizontal distance (m)")
        plt.ylabel("Vertical distance (m)")
        plt.title("Projectile motion of baseball using various ODE numerical methods")
        plt.legend(["Euler's Method", "Euler-Cromer Method", "Midpoint Method"], loc="upper left")
        plt.show()


    if method == "all":
        return "All methods were used, only a plot is returned"
    else:
        return distance


####################
# projectileFence() function models the projective motion of a baseball
#   after being struck with the option of accounting for air resistance
#
# Inputs:
# vel - int/float velocity in m/s
# phi - int/float launch angle in degrees
# tau - int/float time step in seconds
# air - 0 for no air resistance; 1 for air resistance
# method - string e for Euler's method; ec for Euler-Cromer method; m for midpoint method
# plot - 0 for no plotting; 1 for plotting
#
# Output:
# - height - float theoretical height of the ball at 400 meters (can be negative)
# - Plot of projectile motion
####################
def projectileFence(vel, phi, tau, air, method, plot):
    redundant = 1
    # CONSTANTS
    G = 9.81
    Cd = 0.35
    RHO = 1.2
    D = 0.074
    A = np.pi*(D/2)**2
    M = 0.145

    # Initial Conditions
    rx = 0
    ry = 1
    r = [rx, ry]
    phi = phi*np.pi/180
    vx = vel*np.cos(phi)
    vy = vel*np.sin(phi)
    v = [vx, vy]

    # Value Memory
    velList = [v]
    posList = [r]

    # First alpha instance
    if air == 1:
        alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
    else:
        alpha = 0

    ##
    # Euler loop
    ##
    if method == "e":
        while redundant == 1:
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Termination criterion
            if rx >= 121.92 or ry < 0: # 400 feet in meters or below 0 feet
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

    ##
    # Euler-Cromer loop
    ##
    elif method == "ec":
        while redundant == 1:
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*vx
            ry = ry+tau*vy
            # Termination criterion
            if rx >= 121.92 or ry < 0: # 400 feet in meters or below 0 feet
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0
    ##
    # Midpoint loop
    ##
    if method == "m":
        while redundant == 1:
            # Used for averaging velocities
            vxPrev = vx
            vyPrev = vy
            # Velocity step
            vx = vx*(1-tau*alpha)
            vy = vy*(1-tau*alpha)-tau*G
            v = [vx, vy] # Update velocity vector
            # Position step
            rx = rx+tau*(vx+vxPrev)/2
            ry = ry+tau*(vy+vyPrev)/2
            # Termination criterion
            if rx >= 121.92 or ry < 0: # 400 feet in meters or below 0 feet
                break
            r = [rx, ry] # Update position vector
            velList.append(v) # Update memory
            posList.append(r) # Update memory
            # Update alpha
            if air == 1:
                alpha = (Cd*RHO*A*np.linalg.norm(v))/(2*M)
            else:
                alpha = 0

    # Height at 400 feet
    lastIndex = len(posList) - 1
    height = posList[lastIndex][1]

    # Plotting
    if plot == 1:
        xcomp = []
        ycomp = []
        i = 0
        while i < len(posList):
            xcomp.append(posList[i][0])
            ycomp.append(posList[i][1])
            i += 1
        plt.plot(xcomp, ycomp)
        plt.show()


    return height


####################
# abhr() function utilizes projectile() function in order to calculate
#   the ab/hr ratio depending on a standard distribution of launch angle
#   and velocity
#
# Input:
# atBat - int number of simulated at bats to be run
# air - int 0 for no air resistance; 1 for air resistance
# method - string e for Euler's method; ec for Euler-Cromer method; m for midpoint method
# fence - float height of fence, if 0, function will use projectile() and determine
#   homeruns based on if the ball has passed 400 meters before hitting the ground;
#   if fence != 0, function will use projectileFence() and will deal a homerun if
#   the ball height at 400 feet exceeds the height of the fence
#
# Output:
# ab/hr - float at bat homerun ratio
####################
def abhr(atBats, air=1, method="e", fence=0.0):
    # CONSTANTS
    MUV = 44.704 # 100mph in m/s
    MUPHI = 45 # degrees
    SIGMAV = 6.706 # 15mph in m/s
    SIGMAPHI = 10 # degrees

    # Ratio
    ab = 0
    hr = 0

    for i in range(atBats):
        ab += 1
        v0 = SIGMAV * np.random.randn() + MUV
        phi0 = SIGMAPHI * np.random.randn() + MUPHI
        if fence == 0:
            distance = projectile(v0,phi0,0.1,air,method,0)
            if distance >= 121.92: # 400 feet in m
                hr += 1
        else:
            height = projectileFence(v0,phi0,0.1,air,method,0)
            if height > fence:
                hr += 1

    # Zero division exception handling
    if hr == 0:
        return "No homeruns"

    return ab/hr

##
# Part 1 Plotting (Figure 2.3)
##
print(projectile(50,45,0.1,1,"all",1))


##
# Part 2 AB/HR ratio
##
print("The AB/HR ratio of the RDH is: ", abhr(10000))

##
# Part 3 Fence
##

# How AB/HR changes with fence height
height = 0.5
heightList = []
abhrList = []
while height < 15.5:
    abhrList.append(abhr(1000, 1, "e", height))
    heightList.append(height)
    height += 0.5
plt.plot(heightList, abhrList)
plt.xlabel("Height of fence (m)")
plt.ylabel("AB/HR ratio")
plt.title("AB\HR ratio vs fence height")
plt.show()

# Fence height for AB/HR > 10
height = 0
highABHR = False
heightList = []
# Average 100 heights that yield an AB/HR greater than 10
for _ in range(100):
    # Loop until fence height is high enough to achieve AB/HR ratio greater than 10
    while not highABHR:
        if abhr(1000, 1, "e", height) > 10 :
            highABHR = True
        else:
            height += 1 # Increase height by 1
    heightList.append(height)
print("The fence height required for the AB\HR ratio of the RDH to be above 10 is: ", sum(heightList)/len(heightList))
