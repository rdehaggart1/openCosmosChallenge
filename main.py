#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 08:37:34 2021

@author: rory_haggart
"""

#TODO:
    # create a dashboard type thing that allows for variable parameters
    # print velocty? altitude?
    # add support for solar angles in RHP
    # add support for re-drawing when a parameter changes
    
# for plotting purposes
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib.patches import Arrow
from matplotlib import animation
import matplotlib

# for manipulation of data
import math
import numpy as np

# for system exit control
import sys

# GLOBAL CONSTANTS
# If we wanted a different central body, we could change these or create
    # a set of classes for different bodies
COLOR = "#0066aa"
NAME = "EARTH"

RADIUS = 6371       # radius of parent [km] (R_Earth = 6371)
MASS = 5.972e24     # mass of Earth [kg] (M_Earth = 5.972e24)
GRAV = 6.67408e-11  # universal gravitational constant [m3 kg-1 s-2]  
G = 1371            # incident solar radiation [W/m^2] (G_Earth = 1371)

# standard gravitational parameter of Earth [km3 s−2]
MU = GRAV * MASS / (1000*1000*1000)  


def main():    
    # semi-major (x-axis), semi-minor (y-axis) are used to define the orbit
    semiMajorAxis = RADIUS + 5000 # [km]
    semiMinorAxis = RADIUS + 5000 # [km]
    
    panelArea = 1               # area of the solar panel [m2]       
    panelAbsorptivity = 0.5     # absorptivity of the solar panel []
    panelEfficiency = 0.15      # efficiency of the solar panel
    
    # angle of the light rays from the sun. 0 -> | 180 <- | 90 ^
    solarAngle = 50 # [degrees]
    
    # change the backend to interactive (qt) so we can see animations
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='qt')
    except:
        pass 
    
    # create axes
    fig, ax = plt.subplots(1, 2)
    orbitAx = ax[0]
    powerAx = ax[1]
    
    FS = 15 # font size for plots
    
    orbitAx.set_facecolor((0, 0, 0)) # black background
    
    # grid
    orbitAx.grid()
    powerAx.grid()
    
    #fig.set_dpi(100)
    fig.set_size_inches(15, 8)
    fig.tight_layout(pad=7.0)
    
    # make sure out central body is actually shown as a circle
    orbitAx.set_aspect('equal')
    
    # plot and animate the satellite orbit around earth
    solarPanel = panel(fig, ax, 
                       panelArea, panelAbsorptivity, panelEfficiency,
                       semiMajorAxis, semiMinorAxis)
    
    # titles
    orbitAx.set_title("Perigee: {:.1f}km | Apogee: {:.1f}km".format(
                       solarPanel.perigee, solarPanel.apogee),
                       fontsize=FS)
    powerAx.set_title("Power", fontsize=FS)
    
    # labels
    orbitAx.set_xlabel("[km]", fontsize = FS)
    orbitAx.set_ylabel("[km]", fontsize = FS)
    orbitAx.tick_params(axis='both', which='major', labelsize=FS)
    powerAx.set_xlabel("Time [h]", fontsize = FS)
    powerAx.set_ylabel("Power Output [W]", fontsize = FS)
    powerAx.tick_params(axis='both', which='major', labelsize=FS)
        
    # make sure everything is visibile
    orbitAx.set_xlim([-1.1*(solarPanel.focus + solarPanel.semiMajor), 
                       1.1*(solarPanel.focus + solarPanel.semiMajor)])
    orbitAx.set_ylim([-1.1*solarPanel.semiMinor, 
                       1.1*solarPanel.semiMinor])
    powerAx.set_xlim(0, max(solarPanel.time))
    
    # plot arrows to represent the solar rays
    plotSun(fig, orbitAx, solarAngle)
    
    # show our plot
    plt.show()

    # calculate and plot the power output from the solar panel
    solarPanel.powerOutput(solarAngle)
    
    # animate!
    solarPanel.animate()

"""
   @brief  class structure for the solar panel. plot orbit & power out, and animate
   @param  figure handle
   @param  axis handle
   @param  absoption coefficient of the solar panel
   @param  efficiency of the solar panel
   @param  the semi-major axis of the orbit [km]
   @param  the semi-minor axis of the orbit [km]
"""
class panel:
    # initialise the scenario
    def __init__(self,
                 fig=None, 
                 ax=None,
                 area=1,
                 absorptivity=1,
                 efficiency = 0.15,
                 semiMajor=2*RADIUS, 
                 semiMinor=2*RADIUS):
        
        # figures
        self.fig = fig
        self.orbitAx = ax[0]
        self.powerAx = ax[1]
        
        # panel parameters
        self.area = area
        self.absorptivity = absorptivity
        self.efficiency = efficiency
        
        # orbit parameters
        self.semiMajor = semiMajor
        self.semiMinor = semiMinor
        
        # get the orbital data (position, velocity, angle)
        self.getOrbit()
        
        # central body
        central = Ellipse([0,0], RADIUS*2, RADIUS*2, 
                          linewidth=1, fill=1, color=COLOR)
        
        # satellite orbit, with the focus located at centre of central body
        orbit = Ellipse([self.focus,0], self.semiMajor*2, self.semiMinor*2, 
                        linewidth=1, fill=0, color='g')
        
        # size of the satellite and panel
        self.satWidth = 2000
        self.satHeight = 500
        self.panelWidth = 1.5*self.satWidth
        self.panelHeight = self.satHeight/2
        
        # satellite representation
        self.satellite = Rectangle((0, 0), self.satWidth, self.satHeight, 
                                   fc='w')
        # solar panel representation for clearer visualisation
        self.panel = Rectangle((0, 0), self.panelWidth, self.panelHeight, 
                               fc='#EEA533')
    
        # a marker for the power curve to trace data
        self.marker = Ellipse([0,0], 0.2, 0.2, 
                              linewidth=1, fill=1, color='k')
        
        # draw the central body and the satellite orbit
        self.orbitAx.add_artist(orbit)
        self.orbitAx.add_artist(central)
    
    # use the orbit definition to get the coords, angles, and speeds through orbit
    def getOrbit(self):
        self.X = [] # x position
        self.Y = [] # y position
        self.A = [] # angle from central body
        self.V = [] # velocity
        self.speedScaling = []  # a list for the speed up/slow down of orbit
        
        # distance from the centre of the orbit to the focus
        self.focus = np.sqrt(self.semiMajor**2 - self.semiMinor**2)
        # eccentricity of orbit
        self.eccentricity = self.focus / self.semiMajor
    
        # perigee and apogee (orbital altitude, not radius)
        self.perigee = self.semiMajor - self.focus - RADIUS # [km]
        self.apogee = self.semiMajor + self.focus - RADIUS  # [km]
    
        # might be concerning if satellite is crashing into central body...
        if self.perigee < 0:
            plt.close()
            sys.exit("Negative perigee. Exiting")
            
        # point of minimum radius (rMin) = point of maximum velocity (vMax) 
        rMin = np.sqrt((self.semiMajor - self.focus)**2)
        # v = sqrt(mu * (2/r - 1/a)) for an eliptical orbit
        self.vMax = np.sqrt(MU * (2/rMin - 1/self.semiMajor))
        
        # orbital period = 2*pi*sqrt(a^3 / mu)
        self.orbitalPeriod = 2 * np.pi * np.sqrt(self.semiMajor**3 / MU)
        
        # for 360 data points around the ellipse
        for i in range(360):
            # elipse is parameterised with x = semiMajor * sin (t), y = semiMinor * cos(t)
            # and x is offset by the focal length
            x = self.focus + self.semiMajor * np.sin(np.radians(i))
            y = self.semiMinor * np.cos(np.radians(i))
            
            # distance between satellite centre and central body at this point
            r = np.sqrt(x**2 + y**2)   
           
            # angle between satellite and central body at this point
                # 0 radians corresponds to x-axis in this representation
            theta = math.atan2(y,x)
            
            # calculate the tangential speed of the satellite at this point
            vTangential = np.sqrt(MU * (2/r - 1/self.semiMajor))
            
            # normalise the speed against the point of maximum velocity
            speed = vTangential/self.vMax
            
            # append the data for this point in space
            self.X.append(x)
            self.Y.append(y)
            self.A.append(theta)           
            self.V.append(speed)
            # for animation, we change the frame interval based on relative
                # speed of satellite. this list holds that data for each point
            self.speedScaling.append(10/self.V[i])
            
               
        # the speed scaling list works by changing the animation speed. If
            # we used a linear timescale, the marker would advance through
            # time at a variable rate. by using the speed scaling array
            # to build up the timeline, we account for this and allow the 
            # marker to advance at a constant rate
        self.time = np.cumsum([0] + self.speedScaling[:-1])
        
        # the speedscaling time doesn't match orbital period, so scale all values
        scaleFactor = (self.orbitalPeriod/self.time[-1])/3600
        self.time = [a*scaleFactor for a in self.time]
    
    # use the orbital data to calculate the corresponding power output
    def powerOutput(self, solarAngle):
        # ready to store the power output through the orbit
        self.powerOutput = []
        
        # for each of the 360 data points round the ellipse
        for i in range(360):
            # convert from -180:180 to a 0:360 representation
            angle = -np.degrees(self.A[i])+180
            
            # find the angular difference between the satellite and solar rays            
            diff = angle + solarAngle 
            
            # convert back to -180:180 in terms of difference for trig functions
            if diff > 180:
                diff = -1*(360 - diff)
            
            # if panel facing sun to some extent
            if(abs(diff) < 90):
                # radiation absorbed by a plate:
                    # incident radiation * area * absorptivity * cos(angle)
                absorbedPower =  G * self.area * self.absorptivity * np.cos(np.radians(diff))
            else:
                # one the solar panel is > 90 rotated relative to the 
                    # incoming radiation, nothing can be absorbed
                absorbedPower = 0
            
            # then the power output is the input * efficiency
            powerOutput = (absorbedPower) * self.efficiency
            self.powerOutput.append(powerOutput)
        
        # plot the power output of the array against the orbit timesteps
        self.powerAx.plot(self.time, self.powerOutput, color='#666666')
        
        # change the size of the marker ellipse to be circular w.r.t the
        # axis scaling
        self.marker.width = max(self.time)/50
        self.marker.height = max(self.powerOutput)/50   
    
    # initialise the animation by adding the patches to the plot
    def init(self):
        # add the satellite and panel to the figure
        self.orbitAx.add_patch(self.satellite)
        self.orbitAx.add_patch(self.panel)
        self.powerAx.add_patch(self.marker)
        return self.satellite,self.panel,self.marker,
      
    
    # animation function for the satellite
    def orbitAnimation(self, i):    
        # change the speed of the animation to reflect the velocity changes
            # x/speed means x is the shortest interval (i.e. lower x -> faster)
        self.anim.event_source.interval = self.speedScaling[i]
        
        # move the satellite to the new point on the ellipse trajectory
        self.satellite.set_xy([self.X[i] - self.satWidth/2, 
                               self.Y[i] - self.satHeight/2])
        # position the panel representation to the 'outside' of the satellite
        self.panel.set_xy([self.X[i] - self.panelWidth/2, 
                           self.Y[i] + self.satHeight/2])
        
        # move the marker along the power curve to show where we are
        self.marker.set_center([self.time[i], self.powerOutput[i]])
        
        # get the axis transformation data and use this to transform to display coords
        ts = self.orbitAx.transData
        coords = ts.transform([self.X[i], self.Y[i]])
        # perform a rotation relative to the axes of theta radians (calculated above)
            # and then a further 90 degrees so the 'face' of the panels is away from Earth
        tr = matplotlib.transforms.Affine2D().rotate_around(coords[0],coords[1], self.A[i] - np.pi/2)
        t = ts + tr
        
        # perform the rotative transformation
        self.satellite.set_transform(t)
        self.panel.set_transform(t)
        
        return self.satellite,self.panel,self.marker,
    
    def animate(self):
        # execute the animation
        self.anim = animation.FuncAnimation(self.fig, self.orbitAnimation, 
                                            init_func=self.init, 
                                            frames=360, 
                                            blit=True)
        #self.anim.save(r'animation.gif', fps=10)
    
    


"""
   @brief  draw a set of arrows that represent the incoming solar rays
   @param  figure handle
   @param  axis handle
   @param  the angle of the incoming solar rays [degrees] where 0: ->
"""
def plotSun(fig, ax, angle):
    # get the axis limits so we know what we're working with
    xLim = ax.get_xlim()
    yLim = ax.get_ylim()
    
    # the axis rectangle has its corners on the circumference of this circle
    boundingCircle = np.sqrt(xLim[1]**2 + yLim[1]**2)
    
    # list for the origins of the light rays
    rayOrigins = [None] * 3

    # so we can safely start the sun arrow on this circle
    rayOrigins[0] = (boundingCircle * np.sin(np.radians(270 - angle)), boundingCircle * np.cos(np.radians(270 - angle)))
    
    raySep = 10000
    
    if angle > 0:
        # a gradient value to maintain separation between each ray
        gradient = 1+(rayOrigins[0][1] / rayOrigins[0][0])
            
        # calculate separation in x,y
        raySepY = raySep / gradient 
        raySepX = raySep - raySepY
    
        # find the origins of the other light rays
        rayOrigins[1] = (rayOrigins[0][0] - raySepX, rayOrigins[0][1] + raySepY) 
        rayOrigins[2] = (rayOrigins[0][0] + raySepX, rayOrigins[0][1] - raySepY) 
    else:
        # a gradient value to maintain separation between each ray
        gradient = 1-(rayOrigins[0][1] / rayOrigins[0][0])
            
        # calculate separation in x,y
        raySepY = raySep / gradient 
        raySepX = raySep - raySepY
    
        # find the origins of the other light rays
        rayOrigins[1] = (rayOrigins[0][0] - raySepX, rayOrigins[0][1] - raySepY) 
        rayOrigins[2] = (rayOrigins[0][0] + raySepX, rayOrigins[0][1] + raySepY) 
    
    # plot the light rays
    for origin in rayOrigins:
        arrow = Arrow(origin[0], origin[1], 
                      -rayOrigins[0][0], -rayOrigins[0][1], 
                      width = 3000, color='y', alpha = 0.5, zorder=10)
        ax.add_artist(arrow)

if __name__ == "__main__":
    main()