#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 08:37:34 2021

@author: rory_haggart
"""

# for plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import animation

import numpy as np

def main():
    # change the backend to interactive (qt) so we can see animations
    try:
        import IPython
        shell = IPython.get_ipython()
        shell.enable_matplotlib(gui='qt')
    except:
        pass 
    
    semiMajorAxis = 20000 # km
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
    
    # CONSTANTS
    R_Earth = 6371          # radius of Earth [km]
    M_Earth = 5.972e24      # mass of Earth [kg]
    G = 6.67408e-11         # universal gravitational constant [m3 kg-1 s-2]
    mu_Earth = G * M_Earth  # standard gravitational parameter of Earth [m3 s−2]
    
    
    # define earth
    earth = Ellipse([0,0], R_Earth*2, R_Earth*2, angle=0, linewidth=1, fill=1, color="#0088ff")
    
    # define satellite orbit, with the focus located at centre of earth
        # todo: include an angle?
    orbit = Ellipse([focus,0], semiMajor*2, semiMinor*2, angle=0, linewidth=1, fill=0)
    
    # define the satellite representation (a circle)
    satellite = plt.Circle((1000, -1000), 1000, fc='k')
    
    # create axes
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    
    # plot the earth and the satellite orbit
    ax.add_artist(orbit)
    ax.add_artist(earth)
    
    # make sure everything is visibile
    ax.set_xlim([-1.1*(focus + semiMajor), 1.1*(focus + semiMajor)])
    ax.set_ylim([-1.1*semiMinor, 1.1*semiMinor])
    ax.grid()
    #fig.set_dpi(100)
    #fig.set_size_inches(7, 6.5)
    
    # todo:
        # calculate maximum velocity to normalise all velocities
        # work out how to vary speed of animation
        # do some perigee/apogee stuff
    
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
        
        # distance between satellite centre and Earth centre
        r = np.sqrt(x**2 + y**2)    
        
        vTangential = np.sqrt(mu_Earth * (2/r - 1/semiMajor))
        
        # move the satellite to the new point on the ellipse trajectory
        satellite.center = (x, y)
        return satellite,
    
    # execute the animation
    anim = animation.FuncAnimation(fig, animate, 
                                   init_func=init, 
                                   frames=360, 
                                   interval=50,
                                   blit=True)
    
    plt.show()

if __name__ == "__main__":
    main()