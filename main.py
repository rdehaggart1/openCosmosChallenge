#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 08:37:34 2021

@author: rory_haggart
"""

#TODO:
    # perigee apogee stuff
    # create a dashboard type thing that allows for variable parameters
    # print orbital parameters that aren't tunable (eccentricity, ap/per, etc)
    # print velocty? altitude?
    # create a 'orbital parent' class? allows us to easily switch between
        # bodies. e.g. change distance from sun, radius, gravity, etc
    # add the body name into parent somehow
    
# for plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib import animation
import matplotlib

import numpy as np

    
# GLOBAL CONSTANTS
RADIUS = 6371       # radius of parent [km] (R_Earth = 6371)
MASS = 5.972e24     # mass of Earth [kg] (M_Earth = 5.972e24)
GRAV = 6.67408e-11  # universal gravitational constant [m3 kg-1 s-2]
MU = GRAV * MASS    # standard gravitational parameter of Earth [m3 s−2]

COLOR = "#0088ff"
NAME = "EARTH"

def main():
    # change the backend to interactive (qt) so we can see animations
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='qt')
    except:
        pass 
    
    semiMajorAxis = 25000 # km
    semiMinorAxis = 18000 # km

    # calculate the distance from the centre of the orbit to the focus
    focus = np.sqrt(semiMajorAxis**2 - semiMinorAxis**2)
    eccentricity = focus / semiMajorAxis
    print(eccentricity)
    
    # plot the satellite orbit around earth
    plotOrbit(semiMajorAxis,semiMinorAxis, focus)


"""
   @brief  draw the earth and the orbit of satellite as defined by params
   @param  the semi-major axis of the orbit
   @param  the semi-minor axis of the orbit
   @param  the (positive-x) focal point of the orbit
"""
def plotOrbit(semiMajor, semiMinor, focus):
    # todo: animate the orbit - preferably keep the whole orbit visible 
        # and have the point moving around the path
    
    # define earth
    earth = Ellipse([0,0], RADIUS*2, RADIUS*2, angle=0, linewidth=1, fill=1, color=COLOR)
    
    # define satellite orbit, with the focus located at centre of earth
        # todo: include an angle?
    orbit = Ellipse([focus,0], semiMajor*2, semiMinor*2, angle=0, linewidth=1, fill=0)
    
    satWidth = 2000
    satHeight = 500
    
    # define the satellite representation (a circle)
    satellite = Rectangle((0, 0), satWidth, satHeight, fc='k')
    
    # create axes
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    
    # plot the earth and the satellite orbit
    ax.add_artist(orbit)
    ax.add_artist(earth)
    
    labelX = -1 * RADIUS * np.sin(np.pi/4)
    labelY = -1 * RADIUS * np.cos(np.pi/4)

    # make sure everything is visibile
    ax.set_xlim([-1.1*(focus + semiMajor), 1.1*(focus + semiMajor)])
    ax.set_ylim([-1.1*semiMinor, 1.1*semiMinor])
    ax.grid()
    #fig.set_dpi(100)
    #fig.set_size_inches(7, 6.5)
    
    rMin = np.sqrt((semiMajor - focus)**2)
    vMax = np.sqrt(MU * (2/rMin - 1/semiMajor))
    
    # initialise the satellite animation
    def init():
        # add the satellite to the figure
        ax.add_patch(satellite)
        return satellite,
    
    # animation function for the satellite
    def animate(i):
        # elipse is parameterised with x = semiMajor * sin (t), y = semiMinor * cos(t)
        # and x is offset by the focal length
        x = focus + semiMajor * np.sin(np.radians(i))
        y = semiMinor * np.cos(np.radians(i))
        
        # distance between satellite centre and Earth centre at this point
        r = np.sqrt(x**2 + y**2)   
       
        # angle between satellite and Earth centre at this point
        theta = np.arctan(y/x)
        
        # calculate the tangential speed of the satellite at this point
        vTangential = np.sqrt(MU * (2/r - 1/semiMajor))
        
        # normalise the speed against the point of maximum velocity
        speed = vTangential/vMax
    
        # change the speed of the animation to reflect the velocity changes
            # x/speed means x is the shortest interval (i.e. lower x -> faster)
        anim.event_source.interval = 10/speed
        
        # move the satellite to the new point on the ellipse trajectory
        satellite.set_xy([x - satWidth/2, y - satHeight/2])
    
        ts = ax.transData
        coords = ts.transform([x,y])
        tr = matplotlib.transforms.Affine2D().rotate_around(coords[0],coords[1], theta - np.pi/2)
        t= ts + tr
    
        satellite.set_transform(t)
        
        
        return satellite,
        
    # execute the animation
    anim = animation.FuncAnimation(fig, animate, 
                                   init_func=init, 
                                   frames=360, 
                                   blit=True)
    
    plt.show()

if __name__ == "__main__":
    main()