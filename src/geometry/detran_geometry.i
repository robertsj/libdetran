//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_geometry.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Mesh.hh"
%include "Mesh1D.hh"
%include "Mesh2D.hh"
%include "Mesh3D.hh"
%include "PinCell.hh"
%include "Assembly.hh"
%include "Core.hh"

%template(MeshSP)     detran::SP<detran::Mesh>;
%template(Mesh1DSP)   detran::SP<detran::Mesh1D>;
%template(Mesh2DSP)   detran::SP<detran::Mesh2D>;
%template(Mesh3DSP)   detran::SP<detran::Mesh3D>;
%template(PinCellSP)  detran::SP<detran::PinCell>;
%template(AssemblySP) detran::SP<detran::Assembly>;
%template(CoreSP)     detran::SP<detran::Core>;

namespace detran
{
// Save old versions for later post process (plot) use
//----------------------------------------------------------------------------//
//class Mesh1D : public Mesh
//{
//public:
//  Mesh1D(std::vector<int>    xfm, 
//         std::vector<double> xcme, 
//         std::vector<int>    mat_map);
//  Mesh1D(std::vector<double> xfme, 
//         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<int>    xfm, 
//                         std::vector<double> xcme, 
//                         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<double> xfme, 
//                         std::vector<int>    mat_map); 
//  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  %pythoncode 
//  %{
//      
//  def check(self) :
//    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
//    """  
//    try :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      self.PLOT = True
//    except ImportError :
//      self.PLOT = False
//      
//  def plot_mesh_map(self, map_key) :
//    """ Plot the mesh map.
//    """
//    self.check()  
//    if self.PLOT :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      # Get the map.
//      m_map = np.asarray(self.mesh_map(map_key))
//      x = self.mesh_axes()
//      plt.plot(x, m_map)
//      plt.title('mesh map')
//      plt.show()
//    else :
//      print "Warning: matplotlib or numpy not found; skipping plot."
//      
//  def plot_flux(self, f) :
//    """ Plot a flux or flux-shaped vector on the mesh."
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = self.mesh_axes()
//    plt.plot(x, f)
//    plt.title('flux map')
//    plt.show()
//  
//  def mesh_axes(self) :
//    """ Get the fine mesh esges defining cells for plotting.
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = np.zeros(self.number_cells_x())
//    x[0] = 0.5*self.dx(0)
//    for i in range(1, self.number_cells_x()) :
//      x[i] = x[i-1] + self.dx(i)
//    return x
// 
//  %}  
//  
//};
//
////----------------------------------------------------------------------------//
//class Mesh2D : public Mesh
//{
//public:
//  Mesh2D(std::vector<int>    xfm, 
//         std::vector<int>    yfm, 
//         std::vector<double> xcme, 
//         std::vector<double> ycme, 
//         std::vector<int>    mat_map);
//  Mesh2D(std::vector<double> xfme, 
//         std::vector<double> yfme, 
//         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<int>    xfm, 
//                         std::vector<int>    yfm, 
//                         std::vector<double> xcme, 
//                         std::vector<double> ycme, 
//                         std::vector<int>    mat_map);
//  static SP<Mesh> Create(std::vector<double> xfme, 
//                         std::vector<double> yfme, 
//                         std::vector<int>    mat_map); 
//  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
//  
//  %pythoncode 
//  %{
//    
//  def check(self) :
//    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
//    """  
//    try :
//      import matplotlib.pyplot as plt
//      import numpy as np
//      self.PLOT = True
//    except ImportError :
//      self.PLOT = False
//      
//  def plot_mesh_map(self, map_key, edges=False) :
//    """ Plot the mesh map.
//    """
//    self.check()  
//    if self.PLOT :
//    
//      # Imports
//      import matplotlib
//      import matplotlib.pyplot as plt
//      import numpy as np
//      
//      # Get the map.
//      m_map = np.asarray(self.mesh_map(map_key))
//      if len(m_map) == 0 :
//        print "Oops! Map has size zero!"
//        return
//            
//      # Reshape it.
//      m_map2 = m_map.reshape(self.number_cells_x(), self.number_cells_y())
//      
//      # Define random color map.
//      colormap = matplotlib.colors.ListedColormap ( np.random.rand ( 10000,3))
//      
//      # Get the mesh axes and then make a grid of them for plotting.
//      x, y = self.mesh_axes()
//      X, Y = np.meshgrid(x, y)
//      
//      if edges == False :
//        edges = None
//      else :
//        edges = "black"
//
//      # Create the plot a
//      plt.pcolor(X, Y, m_map2, cmap=colormap)
//      plt.show()
//      
//    else :
//      print "Warning: matplotlib or numpy not found; skipping plot."
//      
//  def plot_flux(self, f) :
//    """ Plot a flux or flux-shaped vector on the mesh."
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    # Reshape flux.
//    f = f.reshape(self.number_cells_x(), self.number_cells_y())
//    #print f
//    # Plot me. 
//    plt.imshow(f, interpolation='nearest', cmap=plt.cm.hot)
//    plt.title('flux map')
//    plt.colorbar()
//    plt.show()
//  
//  def mesh_axes(self) :
//    """ Get the fine mesh esges defining cells for plotting.
//    """
//    import matplotlib.pyplot as plt
//    import numpy as np
//    x = np.zeros(self.number_cells_x()+1)
//    y = np.zeros(self.number_cells_y()+1)
//    for i in range(0, self.number_cells_x()) :
//      x[i + 1] = x[i] + self.dx(i)
//    for j in range(0, self.number_cells_y()) :
//      y[j + 1] = y[j] + self.dy(j)
//    print y
//    return (x, y)
//         
//  %}
//};
  

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of detran_geometry.i
//---------------------------------------------------------------------------//





