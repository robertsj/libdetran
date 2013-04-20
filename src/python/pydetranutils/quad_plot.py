# This provides utilities for plotting things on a 
# 1D or 2D mesh or a slice of a 3D mesh.
import numpy as np
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt

try :
    from matplotlib import rc
    rc('text', usetex=True)
    rc('font', family='serif')
except ImportError :
    print "Warning: LaTeX labels being skipped"

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
     plt.xlabel('$\mu$')
     plt.ylabel('weight')
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
     myplot = ax.scatter(mu,eta,xi,c=wt, s=100000*wt**2, marker='^')
     ax.set_xlabel('$\mu$')
     ax.set_ylabel('$\eta$')
     ax.set_zlabel('$\\xi$')
     fig.colorbar(myplot)
     plt.show()
