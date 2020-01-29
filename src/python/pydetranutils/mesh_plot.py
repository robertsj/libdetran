# This provides utilities for plotting things on a 
# 1D or 2D mesh or a slice of a 3D mesh.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib    

def plot_mesh_function(mesh, f, title="", colormap = "hot", edges = False, mybounds = [], myticks = []) :
    """ Plot a mesh function.
    """
    if mesh.dimension() == 1 :
        # get the mesh points
        x = mesh_axes(mesh)
        # plot the map
        plt.plot(x, f)
        
    elif mesh.dimension() == 2 :

        # Get the mesh axes and then make a grid of them for plotting.
        x, y = mesh_axes(mesh)
        X, Y = np.meshgrid(x, y)
        # Reshape the function
        f = f.reshape(mesh.number_cells_x(), mesh.number_cells_y())
        if edges :
            plt.pcolor(X, Y, f, cmap=colormap, edgecolors='k')
        else :
            plt.pcolor(X, Y, f, cmap=colormap)
        plt.axis("scaled") 
        plt.xlabel("x [cm]")
        plt.ylabel("y [cm]")
        if len(myticks) :
          cbar = plt.colorbar(boundaries=mybounds,ticks=myticks)
        else : 
          cbar = plt.colorbar()
    else :
        print "not ready for 3d"
        return
    plt.title(title)
    # show the plot
    plt.show()

def plot_mesh_map(mesh, key, edges = False) :
    """ Plot a mesh map, optionally with edges explicitly displayed.
    """
    # no 3D for now
    if mesh.dimension() == 3 :
        print "not ready for 3d"
        return

    # Check that the map exists and return if missing
    if not mesh.mesh_map_exists(key) :
        print "Mesh map", key, " does not exist"
        return
    
    # Get the map
    map = np.asarray(mesh.mesh_map(key))
    if (mesh.dimension() == 2) :
        # reshape the map
        map = map.reshape(mesh.number_cells_x(), mesh.number_cells_y())
        
    # Choose a random color map for 2D plots
    unique_elements = np.unique(map)
    #num_unique_elements = len(unique_elements)
    num_unique_elements = max(unique_elements)+1
    colormap = matplotlib.colors.ListedColormap(np.random.rand(num_unique_elements, 3))
    bounds = np.linspace(-0.5, num_unique_elements - 0.5, num_unique_elements+1)
    ticks  = bounds[0:num_unique_elements]+0.5
    # Plot
    plot_mesh_function(mesh, map, key, colormap, edges, bounds, ticks)
        
        
def plot_multigroup_flux(mesh, state, edges = False) :
    """ Plot the multigroup fluxes.  
    
        For 1D, they are superimposed on one plot.  In 2D, they
        are split into subfigures for the number of groups.  Obviously,
        this can get cumbersome for many groups, so we kill it at 5+.
    """
    if mesh.dimension() == 1 :
        # get the mesh points
        x = mesh_axes(mesh)
        # plot the map
        plt.plot(x, f)
        
    elif mesh.dimension() == 2 :

        # Get the mesh axes and then make a grid of them for plotting.
        x, y = mesh_axes(mesh)
        X, Y = np.meshgrid(x, y)
        edgec = 'none'
        if edges :
            edgec = 'k'
        plt.pcolor(X, Y, f, cmap=colormap, edgecolors=edgec)
        
    else :
        print "not ready for 3d"
        return
    # show the plot
    plt.show()
    
def mesh_axes(mesh) :
    """ Get the fine mesh points for plotting.
    """
  
    if (mesh.dimension() == 1) :
        # for 1D, we take the cell center points
        x = np.zeros(mesh.number_cells_x())
        x[0] = mesh.dx(0) * 0.5
        for i in range(0, mesh.number_cells_x()-1) :
            x[i + 1] = x[i] + 0.5*(mesh.dx(i) + mesh.dx(i+1))
        return x    
        
    else :
        # for 2D, we take the mesh edges
        x = np.zeros(mesh.number_cells_x()+1)
        y = np.zeros(mesh.number_cells_y()+1)
        for i in range(0, mesh.number_cells_x()) :
            x[i + 1] = x[i] + mesh.dx(i)
        for j in range(0, mesh.number_cells_y()) :
            y[j + 1] = y[j] + mesh.dy(j)
        return (x, y)
    
