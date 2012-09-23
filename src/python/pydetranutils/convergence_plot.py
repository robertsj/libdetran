# This provides utilities for plotting convergence 
# information for solvers

import numpy as np
import matplotlib.pyplot as plt
        
def plot_eigen_convergence(solver) :
    """ Plots the eigenvalue and fission source error per iteration
    """
    # return vec_dbl of eigenvalues
    keffs  = np.asarray(solver.keff_history())
    # return vec_dbl of residuals
    resids = solver.residual_history()
    iterations = range(0, len(keffs))
    assert(len(keffs) == len(iterat))

def plot_convergence(solver) :
    """ Plots the residuals from a linear solve
    """

    # get the residuals as a vec_dbl
    resids = solver.residual_history()
    iters = range(0, len(resids))
    # plot
    plt.plot(resids, iters)
