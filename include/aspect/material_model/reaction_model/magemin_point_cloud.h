/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#ifndef _aspect_material_model_reaction_model_magemin_point_cloud_h
#define _aspect_material_model_reaction_model_magemin_point_cloud_h

#include <cstdint>
#include <vector>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * Input coordinates indexed by nanoflann. Thermodynamic results remain
       * in the hash table; this structure only records how to find them again.
       */
      struct MAGEMinLookupPoint
      {
        std::vector<double> coordinates;
        std::uint64_t assemblage_mask;
        bool uses_fine_cache;
      };


      /**
       * The minimal dataset interface required by nanoflann.
       */
      struct MAGEMinPointCloud
      {
        std::vector<MAGEMinLookupPoint> points;

        std::size_t
        kdtree_get_point_count() const
        {
          return points.size();
        }

        double
        kdtree_get_pt(const std::size_t point_index,
                      const std::size_t coordinate) const
        {
          return points[point_index].coordinates[coordinate];
        }

        template <class BoundingBox>
        bool
        kdtree_get_bbox(BoundingBox &) const
        {
          return false;
        }
      };
    }
  }
}

#endif // !_aspect_material_model_reaction_model_magemin_point_cloud_h
