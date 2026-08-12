/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.
*/

#include <aspect/material_model/reaction_model/unstructured_lookup.h>

#include <deal.II/base/exceptions.h>

#include <algorithm>
#include <cmath>
#include <utility>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      UnstructuredLookup::PointCloud::PointCloud(
        const std::vector<Entry> &table,
        const std::vector<std::size_t> &occupied_slots)
        : table(table),
          occupied_slots(occupied_slots)
      {}



      std::size_t
      UnstructuredLookup::PointCloud::kdtree_get_point_count() const
      {
        return occupied_slots.size();
      }



      double
      UnstructuredLookup::PointCloud::kdtree_get_pt(
        const std::size_t point_index,
        const std::size_t coordinate) const
      {
        return static_cast<double>(
                 table[occupied_slots[point_index]].key[coordinate]);
      }



      UnstructuredLookup::UnstructuredLookup(
        const std::vector<double> &inverse_bin_widths,
        const bool use_kd_tree,
        const std::size_t number_of_neighbors,
        const HashFunction hash_function)
        : inverse_bin_widths(inverse_bin_widths),
          use_kd_tree(use_kd_tree),
          number_of_neighbors(number_of_neighbors),
          hash_function(hash_function),
          validation(KDTreeValidation::accept_closest),
          minimum_validation_neighbors(3),
          maximum_melt_range(0.001),
          maximum_mineral_range(0.10),
          maximum_oxide_range(0.01),
          output_dimension(0),
          occupancy(0),
          shift(64 - 10),
          table(initial_hash_table_size),
          occupied_slots(),
          point_cloud(table, occupied_slots),
          kd_tree(std::make_unique<KDTree>(inverse_bin_widths.size(),
                                           point_cloud))
      {
        AssertThrow(!inverse_bin_widths.empty(),
                    dealii::ExcMessage("An unstructured lookup needs at least one input dimension."));
        AssertThrow(hash_function != HashFunction::ptX ||
                    inverse_bin_widths.size() >= 2,
                    dealii::ExcMessage("A ptX lookup requires pressure and temperature input dimensions."));
        AssertThrow(number_of_neighbors > 0,
                    dealii::ExcMessage("The number of KD-tree neighbors must be positive."));
        for (const double scale : inverse_bin_widths)
          AssertThrow(std::isfinite(scale) && scale > 0.0,
                      dealii::ExcMessage("All unstructured lookup inverse bin widths must be finite and positive."));
        occupied_slots.reserve(initial_hash_table_size);
      }



      void
      UnstructuredLookup::resize_hash_table(const std::size_t new_size)
      {
        AssertThrow(new_size > 0 && (new_size & (new_size - 1)) == 0,
                    dealii::ExcMessage("Unstructured lookup hash-table size must be a power of two."));

        table.clear();
        table.resize(new_size);
        occupied_slots.clear();
        occupied_slots.reserve(new_size);
        occupancy = 0;
        output_dimension = 0;
        shift = 64 - static_cast<std::size_t>(std::log2(new_size));
        kd_tree = std::make_unique<KDTree>(inverse_bin_widths.size(),
                                           point_cloud);
      }



      void
      UnstructuredLookup::set_kdtree_validation(
        const KDTreeValidation new_validation,
        const std::size_t minimum_neighbors,
        const double new_maximum_melt_range,
        const double new_maximum_mineral_range,
        const double new_maximum_oxide_range)
      {
        AssertThrow(minimum_neighbors > 0,
                    dealii::ExcMessage("KD-tree physical validation needs at least one neighbor."));
        AssertThrow(std::isfinite(new_maximum_melt_range) &&
                    new_maximum_melt_range >= 0.0 &&
                    std::isfinite(new_maximum_mineral_range) &&
                    new_maximum_mineral_range >= 0.0 &&
                    std::isfinite(new_maximum_oxide_range) &&
                    new_maximum_oxide_range >= 0.0,
                    dealii::ExcMessage("KD-tree physical validation ranges must be finite and nonnegative."));

        validation = new_validation;
        minimum_validation_neighbors = minimum_neighbors;
        maximum_melt_range = new_maximum_melt_range;
        maximum_mineral_range = new_maximum_mineral_range;
        maximum_oxide_range = new_maximum_oxide_range;
      }



      std::vector<std::int64_t>
      UnstructuredLookup::make_key(const std::vector<double> &inputs) const
      {
        AssertThrow(inputs.size() == inverse_bin_widths.size(),
                    dealii::ExcMessage("Unstructured lookup input dimension does not match its configuration."));

        std::vector<std::int64_t> key(inputs.size());
        for (std::size_t i = 0; i < inputs.size(); ++i)
          {
            AssertThrow(std::isfinite(inputs[i]),
                        dealii::ExcMessage("Unstructured lookup inputs must be finite."));
            key[i] = static_cast<std::int64_t>(
                       std::llround(std::max(0.0, inputs[i]) * inverse_bin_widths[i]));
          }
        return key;
      }



      std::uint64_t
      UnstructuredLookup::ptX_hash(
        const std::vector<std::int64_t> &key) const
      {
        const std::uint64_t pressure_bin = static_cast<std::uint64_t>(key[0]);
        const std::uint64_t temperature_bin = static_cast<std::uint64_t>(key[1]);
        std::uint64_t result =
          (temperature_bin << 32) | (pressure_bin & 0xffffffffULL);

        for (std::size_t i = 2; i < key.size(); ++i)
          {
            result ^= static_cast<std::uint64_t>(key[i]);
            result *= 0xc6a4a7935bd1e995ULL;
            result ^= result >> 47;
          }

        return result * 11400714819323198485ULL;
      }



      std::uint64_t
      UnstructuredLookup::murmur_hash(
        const std::vector<std::int64_t> &key) const
      {
        std::uint64_t result = 0;
        for (const std::int64_t value : key)
          {
            result ^= static_cast<std::uint64_t>(value);
            result *= 0xc6a4a7935bd1e995ULL;
            result ^= result >> 47;
          }

        return result * 11400714819323198485ULL;
      }



      std::uint64_t
      UnstructuredLookup::hash(const std::vector<std::int64_t> &key) const
      {
        if (hash_function == HashFunction::ptX)
          return ptX_hash(key);

        if (hash_function == HashFunction::automatic)
          return murmur_hash(key);

        AssertThrow(false, dealii::ExcNotImplemented());
        return 0;
      }



      std::size_t
      UnstructuredLookup::find_entry(
        const std::vector<std::int64_t> &key,
        const std::uint64_t full_hash) const
      {
        const std::size_t initial_slot =
          shift == 64 ? 0 : static_cast<std::size_t>(full_hash >> shift);
        for (std::size_t probe = 0; probe < table.size(); ++probe)
          {
            const std::size_t slot =
              (initial_slot + probe) & (table.size() - 1);
            if (!table[slot].occupied)
              return invalid_index;
            if (table[slot].cached_hash == full_hash &&
                table[slot].key == key)
              return slot;
          }
        return invalid_index;
      }



      void
      UnstructuredLookup::insert(const std::vector<double> &inputs,
                                 const std::vector<double> &outputs)
      {
        insert_internal(inputs, outputs, nullptr);
      }



      void
      UnstructuredLookup::insert(
        const std::vector<double> &inputs,
        const std::vector<double> &outputs,
        const MAGEMinMetadata &magemin_metadata)
      {
        insert_internal(inputs, outputs, &magemin_metadata);
      }



      void
      UnstructuredLookup::insert_internal(
        const std::vector<double> &inputs,
        const std::vector<double> &outputs,
        const MAGEMinMetadata *magemin_metadata)
      {
        AssertThrow(!outputs.empty(),
                    dealii::ExcMessage("An unstructured lookup entry needs at least one output value."));
        if (occupancy == 0)
          output_dimension = outputs.size();
        AssertThrow(outputs.size() == output_dimension,
                    dealii::ExcMessage("Unstructured lookup output dimension changed between entries."));
        for (const double value : outputs)
          AssertThrow(std::isfinite(value),
                      dealii::ExcMessage("Unstructured lookup outputs must be finite."));

        if (magemin_metadata != nullptr)
          {
            AssertThrow(std::isfinite(magemin_metadata->melt_fraction),
                        dealii::ExcMessage("MAGEMin melt fraction must be finite."));
            for (const double value : magemin_metadata->mineral_proportions)
              AssertThrow(std::isfinite(value),
                          dealii::ExcMessage("MAGEMin mineral proportions must be finite."));
            AssertDimension(magemin_metadata->mineral_oxides.size(),
                            magemin_metadata->mineral_proportions.size());
            for (const std::vector<double> &phase : magemin_metadata->mineral_oxides)
              for (const double value : phase)
                AssertThrow(std::isfinite(value),
                            dealii::ExcMessage("MAGEMin phase compositions must be finite."));
          }

        std::vector<std::int64_t> key = make_key(inputs);
        const std::uint64_t full_hash = hash(key);
        const std::size_t existing_entry = find_entry(key, full_hash);
        if (existing_entry != invalid_index)
          {
            table[existing_entry].values = outputs;
            if (magemin_metadata != nullptr)
              {
                table[existing_entry].magemin_metadata = *magemin_metadata;
                table[existing_entry].has_magemin_metadata = true;
              }
            return;
          }

        const std::size_t initial_slot =
          shift == 64 ? 0 : static_cast<std::size_t>(full_hash >> shift);
        for (std::size_t probe = 0; probe < table.size(); ++probe)
          {
            const std::size_t slot =
              (initial_slot + probe) & (table.size() - 1);
            if (!table[slot].occupied)
              {
                table[slot].key = std::move(key);
                table[slot].values = outputs;
                table[slot].cached_hash = full_hash;
                table[slot].occupied = true;
                if (magemin_metadata != nullptr)
                  {
                    table[slot].magemin_metadata = *magemin_metadata;
                    table[slot].has_magemin_metadata = true;
                  }

                ++occupancy;
                occupied_slots.push_back(slot);
                kd_tree->addPoint(occupied_slots.size() - 1);
                return;
              }
          }

        AssertThrow(false,
                    dealii::ExcMessage("Unstructured lookup hash table is full."));
      }



      const UnstructuredLookup::MAGEMinMetadata *
      UnstructuredLookup::metadata(const Entry &entry) const
      {
        return entry.has_magemin_metadata ? &entry.magemin_metadata : nullptr;
      }



      bool
      UnstructuredLookup::kdtree_validation(
        const std::vector<std::size_t> &point_indices,
        const std::size_t number_of_points) const
      {
        if (number_of_points == 0)
          return false;
        if (validation == KDTreeValidation::accept_closest)
          return true;
        if (number_of_points < minimum_validation_neighbors)
          return false;

        const Entry &reference = table[occupied_slots[point_indices[0]]];
        if (!reference.has_magemin_metadata)
          return false;
        const MAGEMinMetadata &reference_metadata = reference.magemin_metadata;

        double minimum_melt = reference_metadata.melt_fraction;
        double maximum_melt = reference_metadata.melt_fraction;
        for (std::size_t neighbor = 0; neighbor < number_of_points; ++neighbor)
          {
            const Entry &entry = table[occupied_slots[point_indices[neighbor]]];
            if (!entry.has_magemin_metadata ||
                entry.magemin_metadata.assemblage_mask !=
                reference_metadata.assemblage_mask ||
                entry.magemin_metadata.mineral_proportions.size() !=
                reference_metadata.mineral_proportions.size() ||
                entry.magemin_metadata.mineral_oxides.size() !=
                reference_metadata.mineral_oxides.size())
              return false;

            minimum_melt = std::min(minimum_melt,
                                    entry.magemin_metadata.melt_fraction);
            maximum_melt = std::max(maximum_melt,
                                    entry.magemin_metadata.melt_fraction);
          }
        if (maximum_melt - minimum_melt > maximum_melt_range)
          return false;

        for (std::size_t mineral = 0;
             mineral < reference_metadata.mineral_proportions.size();
             ++mineral)
          {
            double minimum = reference_metadata.mineral_proportions[mineral];
            double maximum = minimum;
            for (std::size_t neighbor = 1; neighbor < number_of_points; ++neighbor)
              {
                const MAGEMinMetadata &neighbor_metadata =
                  table[occupied_slots[point_indices[neighbor]]].magemin_metadata;
                minimum = std::min(minimum,
                                   neighbor_metadata.mineral_proportions[mineral]);
                maximum = std::max(maximum,
                                   neighbor_metadata.mineral_proportions[mineral]);
              }
            if (maximum - minimum > maximum_mineral_range)
              return false;

            if (mineral >= 64 ||
                (reference_metadata.assemblage_mask &
                 (std::uint64_t(1) << mineral)) == 0)
              continue;

            const std::size_t number_of_oxides =
              reference_metadata.mineral_oxides[mineral].size();
            for (std::size_t neighbor = 1; neighbor < number_of_points; ++neighbor)
              if (table[occupied_slots[point_indices[neighbor]]]
                  .magemin_metadata.mineral_oxides[mineral].size() !=
                  number_of_oxides)
                return false;

            for (std::size_t oxide = 0; oxide < number_of_oxides; ++oxide)
              {
                double minimum_oxide =
                  reference_metadata.mineral_oxides[mineral][oxide];
                double maximum_oxide = minimum_oxide;
                for (std::size_t neighbor = 1;
                     neighbor < number_of_points;
                     ++neighbor)
                  {
                    const double value =
                      table[occupied_slots[point_indices[neighbor]]]
                      .magemin_metadata.mineral_oxides[mineral][oxide];
                    minimum_oxide = std::min(minimum_oxide, value);
                    maximum_oxide = std::max(maximum_oxide, value);
                  }
                if (maximum_oxide - minimum_oxide > maximum_oxide_range)
                  return false;
              }
          }

        return true;
      }



      UnstructuredLookup::Result
      UnstructuredLookup::lookup(const std::vector<double> &inputs) const
      {
        const std::vector<std::int64_t> key = make_key(inputs);
        const std::uint64_t full_hash = hash(key);
        const std::size_t exact_entry = find_entry(key, full_hash);
        if (exact_entry != invalid_index)
          return {Source::hash,
                  table[exact_entry].values,
                  metadata(table[exact_entry])
                 };

        if (!use_kd_tree || occupancy == 0)
          return {};

        std::vector<double> query(inputs.size());
        for (std::size_t i = 0; i < inputs.size(); ++i)
          query[i] = std::max(0.0, inputs[i]) * inverse_bin_widths[i];

        const std::size_t requested_neighbors =
          std::min(number_of_neighbors, occupancy);
        std::vector<std::size_t> point_indices(requested_neighbors);
        std::vector<double> squared_distances(requested_neighbors);
        const std::size_t neighbors_found =
          kd_tree->knnSearch(query.data(), requested_neighbors,
                             point_indices.data(), squared_distances.data());

        if (!kdtree_validation(point_indices, neighbors_found))
          return {};

        const Entry &closest_entry =
          table[occupied_slots[point_indices[0]]];
        return {Source::kd_tree,
                closest_entry.values,
                metadata(closest_entry)
               };
      }



      std::size_t
      UnstructuredLookup::size() const
      {
        return occupancy;
      }
    }
  }
}
