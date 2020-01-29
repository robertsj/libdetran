//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  RegionFactory.cc
 *  @brief RegionFactory
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "RegionFactory.hh"
#include "utilities/SoftEquivalence.hh"
#include <cmath>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
RegionFactory::vec_region
RegionFactory::CreatePinCell(SP_pincell pin)
{
  Require(pin);

  // set the pitch, using the x pitch if y not set
  Point pitch = pin->pitch();
  if (detran_utilities::soft_equiv(pitch.y(), 0.0))
    pitch = pitch * Point(1, 0, 1) + Point(0, pitch.x(), 0);

  size_t div = pin->division();

  //--------------------------------------------------------------------------//
  // SURFACES
  //--------------------------------------------------------------------------//

  // cell bounding planes
  vec_surface bounds(6);
  bounds[0] = QF::CreatePlaneX(0.0);        // west
  bounds[1] = QF::CreatePlaneX(pitch.x());  // east
  bounds[2] = QF::CreatePlaneY(0.0);        // south
  bounds[3] = QF::CreatePlaneY(pitch.y());  // north
  if (pitch.z() > 0.0)
  {
    bounds[4] = QF::CreatePlaneZ(0.0);        // bottom
    bounds[5] = QF::CreatePlaneZ(pitch.z());  // top
  }

  // partitioning planes
  vec_surface partitions;
  vec_int left_sense;
  bool hv   = div == PinCell::DIVISION_HV   || div == PinCell::DIVISION_HV_DIAG;
  bool diag = div == PinCell::DIVISION_DIAG || div == PinCell::DIVISION_HV_DIAG;
  if (hv)
  {
    partitions.push_back(QF::CreatePlaneY(0.5*pitch.y()-pin->pin_center().y()));
    left_sense.push_back(1);
  }
  if (diag)
  {
    partitions.push_back(QF::CreatePlane(pitch.y()/pitch.x(), -1, 0, 0));
    left_sense.push_back(0);
  }
  if (hv)
  {
    partitions.push_back(QF::CreatePlaneX(0.5*pitch.x()-pin->pin_center().x()));
    left_sense.push_back(0);
  }
  if (diag)
  {
    partitions.push_back(QF::CreatePlane(pitch.y()/pitch.x(), 1, 0, pitch.y()));
    left_sense.push_back(0);
  }

  // cylinders
  vec_surface cylinders;
  Point center = pin->pin_center() + pitch * 0.5;
  const vec_dbl &radii = pin->radii();
  for (size_t i = 0; i < radii.size(); ++i)
    cylinders.push_back(QF::CreateCylinderZ(center.x(), center.y(), radii[i]));

  //--------------------------------------------------------------------------//
  // REGIONS
  //--------------------------------------------------------------------------//

  vec_region regions;

  // bounding box for all regions is the same and is just the pin cell + eps
  double eps = 1e-8;
  Point bbox_min = Point(0) - eps;
  Point bbox_max = pitch + eps;

  size_t number_divisions = std::max(1, (int)(2 * partitions.size()));
  size_t number_regions   = (1 + radii.size()) * number_divisions;

  // regions from the center out
  for (size_t r = 0; r < radii.size() + 1; ++r)
  {
    // counter-clockwise around from quadrant 0
    for (size_t a = 0; a < number_divisions; ++a)
    {
      Assert(r < pin->mat_map().size());
      SP_region region = Region::Create(pin->mat_map()[r], bbox_min, bbox_max);

      // cylinder bounds
      if (r < radii.size())
        region->append(cylinders[r], false);  // inside this cylinder
      if (r > 0 && cylinders.size() > 0)
        region->append(cylinders[r-1], true); // outside this cylinder

      // outer bounding surfaces
      if (r == radii.size())
      {
        // \todo not all surfaces needed
        for (size_t b = 0; b < 6; ++b)
          if(bounds[b]) region->append(bounds[b], (bool)((b + 1) % 2));
      }

      if (partitions.size())
      {
        int s0  = (a    ) % partitions.size();
        int s1  = (a + 1) % partitions.size(); // 2/2=1 3/2=1 4/2=2
        bool b0 = (a    ) / partitions.size() ? ! left_sense[s0] :   left_sense[s0];
        bool b1 = ((a + 1) / partitions.size()) % 2 ?   left_sense[s1] : ! left_sense[s1];
//        std::cout << " a = " << a
//                  << " sense0=" << left_sense[s0]
//                  << " sense1=" << left_sense[s1]
//                  << " b0=" << b0
//                  << " b1=" << b1 <<  " " << left_sense[s1]  << " " << (a + 1) /partitions.size()
//                  << std::endl;
        region->append(partitions[s0], b0);
        region->append(partitions[s1], b1);
      }

      regions.push_back(region);
    }
  }

  Ensure(regions.size() == number_regions);
  return regions;
}

//----------------------------------------------------------------------------//
RegionFactory::vec_region
RegionFactory::CreateAssembly(SP_assembly assembly)
{
  Require(assembly);
  vec_region regions;

  // Get all the unique pin regions
  std::vector<vec_region> pin_regions;
  Assembly::vec_pincell pin_cells = assembly->pin_cells();
  for (size_t i = 0; i < pin_cells.size(); ++i)
  {
    pin_regions.push_back(CreatePinCell(pin_cells[i]));
  }
  Point pitch = pin_cells[0]->pitch();

  for (size_t i = 0; i < assembly->dimension(0); ++i)
  {
    for (size_t j = 0; j < assembly->dimension(1); ++j)
    {
      Point origin(pitch.x()*i, pitch.y()*j);
      size_t p = assembly->pincell_map()[i + j*assembly->dimension(0)];
      for (size_t r = 0; r < pin_regions[p].size(); ++r)
      {
        regions.push_back(CreateTranslatedRegion(pin_regions[p][r], origin));
      }
    }
  }

  return regions;
}

//----------------------------------------------------------------------------//
RegionFactory::vec_region
RegionFactory::CreateCore(SP_core core)
{
  Require(core);
  vec_region regions;

//  // Get all the unique regions
//  std::vector<vec_region> assembly_regions;
//  Core::vec_assembly assemblies = core->assemblies();
//  for (size_t i = 0; i < assemblies.size(); ++i)
//  {
//    assembly_regions.push_back(CreatePinCell(assemblies[i]));
//  }
//  Point pitch = assemblies[0]->pitch();
//
//  for (size_t i = 0; i < core->dimension(0); ++i)
//  {
//    for (size_t j = 0; j < core->dimension(1); ++j)
//    {
//      Point O(pitch.x()*i, pitch.y()*j);
//      size_t p = core->assembly_map()[i + j*core->dimension(0)];
//      for (size_t r = 0; r < assembly_regions[p].size(); ++r)
//      {
//        regions.push_back(CreateTranslatedRegion(assembly_regions[p][r], O));
//      }
//    }
//  }

  return regions;
}

//----------------------------------------------------------------------------//
RegionFactory::SP_region
RegionFactory::CreateTranslatedRegion(SP_region r, const Point &origin)
{
  Require(r);
  size_t m = r->attribute("MATERIAL");
  SP_region new_r(new Region(m, r->bound_min()-origin, r->bound_max()-origin));
  Region::SP_node new_n(new CSG_Translation(r->top_node(), origin));
  new_r->append(new_n, INTERSECTION);
  return new_r;
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of RegionFactory.cc
//----------------------------------------------------------------------------//
