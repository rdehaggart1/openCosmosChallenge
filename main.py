#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 29 08:37:34 2021

@author: rory_haggart
"""

# for plotting
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np

def main():    
    semiMajorAxis = 20582 # km
    semiMinorAxis = 17000 # km

    # calculate the distance from the centre of the orbit to the focus
    focus = np.sqrt(semiMajorAxis**2 - semiMinorAxis**2)
    eccentricity = focus / semiMajorAxis
    print(eccentricity)
    
    # plot the satellite orbit around earth
    plotOrbit(semiMajorAxis,semiMinorAxis, focus)

# plot
def plotOrbit(semiMajor, semiMinor, focus):
    rEarth = 6371 # radius of the earth km
    
    # define earth
    earth = Ellipse([0,0], rEarth*2, rEarth*2, angle=0, linewidth=1, fill=1, color="#0088ff")
    
    # define satellite orbit, with the focus located at centre of earth
        # todo: include an angle?
    orbit = Ellipse([focus,0], semiMajor*2, semiMinor*2, angle=0, linewidth=1, fill=0)
    
    # create axes
    fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
    
    # plot the earth and the satellite orbit
    ax.add_artist(orbit)
    ax.add_artist(earth)
    
    # make sure everything is visibile
    ax.set_xlim([-1.1*(focus + semiMajor), 1.1*(focus + semiMajor)])
    ax.set_ylim([-1.1*semiMinor, 1.1*semiMinor])
    plt.show()

if __name__ == "__main__":
    main()