/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#ifndef _aspect_material_model_reaction_model_unstructured_lookup_h
#define _aspect_material_model_reaction_model_unstructured_lookup_h

#include <nanoflann.hpp>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      /**
       * An unstructured lookup that checks a linearly probed hash table first
       * and optionally searches the hash-owned points with a KD-tree.
       */
      class UnstructuredLookup
      {
        public:
          enum class Source
          {
            miss,
            hash,
            kd_tree
          };

          enum class HashFunction
          {
            ptX,
            automatic
          };

          enum class KDTreeValidation
          {
            accept_closest,
            physical_checks
          };

          /**
           * MAGEMin information needed to validate KD-tree neighbours. The
           * thermodynamic result can remain in values; these fields contain
           * only the information used by the physical checks.
           */
          struct MAGEMinMetadata
          {
            std::uint64_t assemblage_mask = 0;
            double melt_fraction = 0.0;
            std::vector<double> mineral_proportions;
            std::vector<std::vector<double>> mineral_oxides;
          };

          struct Result
          {
            Source source = Source::miss;
            std::vector<double> values;
            const MAGEMinMetadata *magemin_metadata = nullptr;
          };

          UnstructuredLookup(const std::vector<double> &inverse_bin_widths,
                             const bool use_kd_tree = true,
                             const std::size_t number_of_neighbors = 4,
                             const HashFunction hash_function = HashFunction::ptX);

          UnstructuredLookup(const UnstructuredLookup &) = delete;
          UnstructuredLookup &operator=(const UnstructuredLookup &) = delete;

          /** Clear the lookup and set a new power-of-two hash-table size. */
          void
          resize_hash_table(const std::size_t new_size);

          void
          set_kdtree_validation(const KDTreeValidation validation,
                                const std::size_t minimum_neighbors = 3,
                                const double maximum_melt_range = 0.001,
                                const double maximum_mineral_range = 0.10,
                                const double maximum_oxide_range = 0.01);

          void
          insert(const std::vector<double> &inputs,
                 const std::vector<double> &outputs);

          void
          insert(const std::vector<double> &inputs,
                 const std::vector<double> &outputs,
                 const MAGEMinMetadata &magemin_metadata);

          Result
          lookup(const std::vector<double> &inputs) const;

          std::size_t
          size() const;

        private:
          struct Entry
          {
            std::vector<std::int64_t> key;
            std::vector<double> values;
            MAGEMinMetadata magemin_metadata;
            std::uint64_t cached_hash = 0;
            bool occupied = false;
            bool has_magemin_metadata = false;
          };

          /** Non-owning nanoflann wrapper around occupied hash-table entries. */
          struct PointCloud
          {
            PointCloud(const std::vector<Entry> &table,
                       const std::vector<std::size_t> &occupied_slots);

            std::size_t
            kdtree_get_point_count() const;

            double
            kdtree_get_pt(const std::size_t point_index,
                          const std::size_t coordinate) const;

            template <class BoundingBox>
            bool
            kdtree_get_bbox(BoundingBox &) const
            {
              return false;
            }

            const std::vector<Entry> &table;
            const std::vector<std::size_t> &occupied_slots;
          };

          using Distance = nanoflann::L2_Simple_Adaptor<double, PointCloud>;
          using KDTree = nanoflann::KDTreeSingleIndexIncrementalAdaptor<
                         Distance, PointCloud, -1, std::size_t>;

          std::vector<std::int64_t>
          make_key(const std::vector<double> &inputs) const;

          std::uint64_t
          hash(const std::vector<std::int64_t> &key) const;

          std::uint64_t
          ptX_hash(const std::vector<std::int64_t> &key) const;

          std::uint64_t
          murmur_hash(const std::vector<std::int64_t> &key) const;

          std::size_t
          find_entry(const std::vector<std::int64_t> &key,
                     const std::uint64_t full_hash) const;

          void
          insert_internal(const std::vector<double> &inputs,
                          const std::vector<double> &outputs,
                          const MAGEMinMetadata *magemin_metadata);

          bool
          kdtree_validation(const std::vector<std::size_t> &point_indices,
                            const std::size_t number_of_points) const;

          const MAGEMinMetadata *
          metadata(const Entry &entry) const;

          static constexpr std::size_t invalid_index =
            std::numeric_limits<std::size_t>::max();
          static constexpr std::size_t initial_hash_table_size = 1024;

          const std::vector<double> inverse_bin_widths;
          const bool use_kd_tree;
          const std::size_t number_of_neighbors;
          const HashFunction hash_function;

          KDTreeValidation validation;
          std::size_t minimum_validation_neighbors;
          double maximum_melt_range;
          double maximum_mineral_range;
          double maximum_oxide_range;

          std::size_t output_dimension;
          std::size_t occupancy;
          std::size_t shift;
          std::vector<Entry> table;

          // Dense KD point number -> occupied slot in the linearly probed table.
          std::vector<std::size_t> occupied_slots;
          PointCloud point_cloud;
          std::unique_ptr<KDTree> kd_tree;
      };
    }
  }
}

#endif
