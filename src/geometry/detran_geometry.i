//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_geometry.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

%module detran_geometry
%{
#include "Definitions.hh"
#include "SP.hh"
#include "Mesh.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"  
#include "PinCell.hh"
#include "Assembly.hh"
#include "Core.hh"
%}

// Load the standard library interfaces.
%include std_vector.i
%include std_string.i

// Load the vector maps.  Note, I am a bit unhappy with
// how it's all used.  They work *if* I declare the class
// interface below.  Otherwise, just including e.g.
// Mesh2D doesn't allow the maps, since I'm using 
// typedefs on the input arguments.  There should be an
// easy way around this, but I'm not a SWIG pro.
%include std_vec_typemap.i
%apply (std::vector<int>    INPUTVECTOR, 
        std::vector<int>    INPUTVECTOR, 
        std::vector<double> INPUTVECTOR, 
        std::vector<double> INPUTVECTOR, 
        std::vector<int>    INPUTVECTOR)
      {(std::vector<int>    xfm, 
        std::vector<int>    yfm, 
        std::vector<double> xcme, 
        std::vector<double> ycme, 
        std::vector<int>    mat_map)}

%include "Definitions.hh"
%include "SP.hh"
%include "Mesh.hh"

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

%template(MeshSP)     detran::SP<detran::Mesh>;
%template(Mesh1DSP)   detran::SP<detran::Mesh1D>;
%template(Mesh2DSP)   detran::SP<detran::Mesh2D>;
%template(Mesh3DSP)   detran::SP<detran::Mesh3D>;
%template(PinCellSP)  detran::SP<detran::PinCell>;
%template(AssemblySP) detran::SP<detran::Assembly>;
%template(CoreSP)     detran::SP<detran::Core>;

namespace detran
{
// Simple redefinition of the basic interfaces along with
// some Python-specific extensions for plotting, etc.

//----------------------------------------------------------------------------//
class Mesh1D : public Mesh
{
public:
  Mesh1D(std::vector<int>    xfm, 
         std::vector<double> xcme, 
         std::vector<int>    mat_map);
  Mesh1D(std::vector<double> xfme, 
         std::vector<int>    mat_map);
  static SP<Mesh> Create(std::vector<int>    xfm, 
                         std::vector<double> xcme, 
                         std::vector<int>    mat_map);
  static SP<Mesh> Create(std::vector<double> xfme, 
                         std::vector<int>    mat_map); 
  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
  %pythoncode 
  %{
      
  def check(self) :
    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
    """  
    try :
      import matplotlib.pyplot as plt
      import numpy as np
      self.PLOT = True
    except ImportError :
      self.PLOT = False
      
  def plot_mesh_map(self, map_key) :
    """ Plot the mesh map.
    """
    self.check()  
    if self.PLOT :
      import matplotlib.pyplot as plt
      import numpy as np
      # Get the map.
      m_map = np.asarray(self.mesh_map(map_key))
      x = self.mesh_axes()
      plt.plot(x, m_map)
      plt.title('mesh map')
      plt.show()
    else :
      print "Warning: matplotlib or numpy not found; skipping plot."
      
  def plot_flux(self, f) :
    """ Plot a flux or flux-shaped vector on the mesh."
    """
    import matplotlib.pyplot as plt
    import numpy as np
    x = self.mesh_axes()
    plt.plot(x, f)
    plt.title('flux map')
    plt.show()
  
  def mesh_axes(self) :
    """ Get the fine mesh esges defining cells for plotting.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.zeros(self.number_cells_x())
    x[0] = 0.5*self.dx(0)
    for i in range(1, self.number_cells_x()) :
      x[i] = x[i-1] + self.dx(i)
    return x
 
  %}  
  
};

//----------------------------------------------------------------------------//
class Mesh2D : public Mesh
{
public:
  Mesh2D(std::vector<int>    xfm, 
         std::vector<int>    yfm, 
         std::vector<double> xcme, 
         std::vector<double> ycme, 
         std::vector<int>    mat_map);
  Mesh2D(std::vector<double> xfme, 
         std::vector<double> yfme, 
         std::vector<int>    mat_map);
  static SP<Mesh> Create(std::vector<int>    xfm, 
                         std::vector<int>    yfm, 
                         std::vector<double> xcme, 
                         std::vector<double> ycme, 
                         std::vector<int>    mat_map);
  static SP<Mesh> Create(std::vector<double> xfme, 
                         std::vector<double> yfme, 
                         std::vector<int>    mat_map); 
  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
  
  %pythoncode 
  %{
    
  def check(self) :
    """ Check for matplotlib and numpy.  Call at top of all plotting methods.
    """  
    try :
      import matplotlib.pyplot as plt
      import numpy as np
      self.PLOT = True
    except ImportError :
      self.PLOT = False
      
  def plot_mesh_map(self, map_key, edges=False) :
    """ Plot the mesh map.
    """
    self.check()  
    if self.PLOT :
    
      # Imports
      import matplotlib
      import matplotlib.pyplot as plt
      import numpy as np
      
      # Get the map.
      m_map = np.asarray(self.mesh_map(map_key))
      if len(m_map) == 0 :
        print "Oops! Map has size zero!"
        return
            
      # Reshape it.
      m_map2 = m_map.reshape(self.number_cells_x(), self.number_cells_y())
      
      # Define random color map.
      colormap = matplotlib.colors.ListedColormap ( np.random.rand ( 10000,3))
      
      # Get the mesh axes and then make a grid of them for plotting.
      x, y = self.mesh_axes()
      X, Y = np.meshgrid(x, y)
      
      if edges == False :
        edges = None
      else :
        edges = "black"

      # Create the plot a
      plt.pcolor(X, Y, m_map2, cmap=colormap)
      plt.show()
      
    else :
      print "Warning: matplotlib or numpy not found; skipping plot."
      
  def plot_flux(self, f) :
    """ Plot a flux or flux-shaped vector on the mesh."
    """
    import matplotlib.pyplot as plt
    import numpy as np
    # Reshape flux.
    f = f.reshape(self.number_cells_x(), self.number_cells_y())
    #print f
    # Plot me. 
    plt.imshow(f, interpolation='nearest', cmap=plt.cm.hot)
    plt.title('flux map')
    plt.colorbar()
    plt.show()
  
  def mesh_axes(self) :
    """ Get the fine mesh esges defining cells for plotting.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    x = np.zeros(self.number_cells_x()+1)
    y = np.zeros(self.number_cells_y()+1)
    for i in range(0, self.number_cells_x()) :
      x[i + 1] = x[i] + self.dx(i)
    for j in range(0, self.number_cells_y()) :
      y[j + 1] = y[j] + self.dy(j)
    print y
    return (x, y)
         
  %}
};
  
//----------------------------------------------------------------------------//
class Mesh3D : public Mesh
{
public:
  Mesh3D(std::vector<int> xfm,
         std::vector<int> yfm,
         std::vector<int> zfm,
         std::vector<double> xcme,
         std::vector<double> ycme,
         std::vector<double> zcme,
         std::vector<int> mat_map);
  Mesh3D(std::vector<double> xfme,
         std::vector<double> yfme,
         std::vector<double> zfme,
         std::vector<int> mat_map);
  static SP<Mesh> Create(std::vector<int> xfm,
                         std::vector<int> yfm,
                         std::vector<int> zfm,
                         std::vector<double> xcme,
                         std::vector<double> ycme,
                         std::vector<double> zcme,
                         std::vector<int> mat_map);
  static SP<Mesh> Create(std::vector<double> xfme,
                         std::vector<double> yfme,
                         std::vector<double> zfme,
                         std::vector<int> mat_map);
  void add_coarse_mesh_map(std::string map_key, std::vector<int> mesh_map);
  void add_mesh_map(std::string map_key, std::vector<int> mesh_map);
};

//----------------------------------------------------------------------------//
class PinCell
{
public:
  PinCell(double pitch, 
          std::vector<double> radii, 
          std::vector<int> mat_map);
  static SP<PinCell> Create(double pitch,
                            std::vector<double> radii,
                            std::vector<int> mat_map);
  SP<Mesh> mesh();
  const Mesh2D& mesh_ref() const;
  void meshify(int number_meshes, bool flag);
}; // end class PinCell

class Assembly
{
public:
  explicit Assembly(int dimension);
  static SP<Assembly> Create(int dimension);
  void add_pincell(SP<detran::PinCell>);
  void finalize(std::vector<int> pincell_map);
  SP<Mesh> mesh();
  const Mesh2D& mesh_ref() const;
}; // end class Assembly

class Core
{
public:
  explicit Core(int dimension);
  static SP<Core> Create(int dimension);
  void add_assembly(SP<detran::Assembly>);
  void finalize(std::vector<int> assembly_map);
  SP<Mesh> mesh();
  const Mesh2D& mesh_ref() const;
}; // end class Core

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of detran_geometry.i
//---------------------------------------------------------------------------//





