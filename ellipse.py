"""
Created on Mon Dec 12 19:07:12 2016
@author: lenne
"""

import numpy as np
import pylab as mp

def plot_ellipse(semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e4,ax=None,plot_kwargs=None,fill=False,fill_kwargs=None):
        '''
                - create an ellipse in polar coordinates then transform to cartesian
                - if given an axes, plot an ellipse with plot_kwargs
                - if not given an axes, create a basic figure and axes and plot
                major keywords are:
                semimaj : length of semimajor axis
                semimin : length of semiminor axis
                phi : angle in radians semimajor axis is above positive x axis
                x_cent : x center
                y_cent : y center
                theta_num : number of points to sample from 0 to 2pi
        '''
        # Generate data for ellipse structure
        theta = np.linspace(0,2*np.pi,theta_num)
        r = semimaj*semimin / np.sqrt((semimin*np.cos(theta))**2 + (semimaj*np.sin(theta))**2)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        data = np.array([x,y])
        R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
        data = np.dot(R,data)
        data[0] += x_cent
        data[1] += y_cent
      

        # Plot!
        if ax == None:
                fig,ax = mp.subplots()

        if plot_kwargs == None:
                ax.plot(data[0],data[1],color='b',linestyle='-')        
        else:
                ax.plot(data[0],data[1],**plot_kwargs)


        if fill == True:
                ax.fill(data[0],data[1],**fill_kwargs)
        
        return data
