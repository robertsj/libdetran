# This provides utilities for plotting things on a 
# 1D or 2D mesh or a slice of a 3D mesh.

import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt

# LaTeX labels are nicer
global __detranuselatex__ 
__detranuselatex__ = False
#try :
#  import os
#  print "Checking for LaTeX for nicer plotting labels..."
#  if (os.system("which latex > NULL")==0) :
#    from matplotlib import rc
#    rc('text', usetex=True)
#    rc('font', family='serif')
#    __detranuselatex__ = True
#except ImportError :
#  pass

def plot_quadrature(quad) :
    """ Plots a quadrature.
    """

    try :
        D = quad.dimension()
    except :
        print "Error getting quadrature dimension... maybe not a quadrature object?"
        return

    if D == 1 :

        # Get the abscissa and weights
        mu = np.asarray(quad.cosines(0))
        wt = np.asarray(quad.weights())
        # Plot
        plt.plot(mu, wt, 'bo')
        if __detranuselatex__ :
            plt.xlabel('$\mu$')
        else :
            plt.xlabel('mu')
        plt.ylabel('weight')
        plt.title(quad.name())
        plt.show()

    else :

        # Get the abscissa and weights
        mu  = np.asarray(quad.cosines(0))
        eta = np.asarray(quad.cosines(1))
        xi  = np.asarray(quad.cosines(2))
        wt  = np.asarray(quad.weights())
        # Plot.  Note, using my (old) version of matplotlib, the colors
        # are not translated to the scatter plot.  The sizes are, but
        # it's not really enough.  What's a good solution?
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        ax.view_init(30, 60) 
        myplot = ax.scatter(mu,eta,xi,c=wt, s=2000*wt, marker='^')
        labels = ['mu','eta','xi']
        if __detranuselatex__ :  
            labels = ['$\mu$', '$\eta$', '$\\xi$']
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_zlabel(labels[2])
        fig.colorbar(myplot)
        plt.title(quad.name())
        plt.show()