// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/AbstractVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

namespace Acts {

class Layer;
class Surface;
class SurfaceArray;
class TrackingVolume;
struct Polyhedron;

static ViewConfig s_viewSensitive = ViewConfig({0, 180, 240});
static ViewConfig s_viewPassive = ViewConfig({240, 280, 0});
static ViewConfig s_viewVolume = ViewConfig({220, 220, 0});
static ViewConfig s_viewGrid = ViewConfig({220, 0, 0});
static ViewConfig s_viewLine = ViewConfig({0, 0, 220});

struct GeometryView3D {
  /// Helper method to draw Polyhedron objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param polyhedron The surface to be drawn
  /// @param viewConfig The drawing configuration
  static void drawPolyhedron(IVisualization3D& helper,
                             const Polyhedron& polyhedron,
                             const ViewConfig& viewConfig = s_viewVolume);

  /// Helper method to draw Surface objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param surface The surface to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param transform An option additional transform
  /// @param ViewConfig The drawing configuration
  static void drawSurface(IVisualization3D& helper, const Surface& surface,
                          const GeometryContext& gctx,
                          const Transform3& transform = Transform3::Identity(),
                          const ViewConfig& ViewConfig = s_viewSensitive);

  /// Helper method to draw SurfaceArray objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param surfaceArray The surface to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param transform An option additional transform
  /// @param sensitiveConfig The drawing configuration for sensitive surfaces
  /// @param passiveConfig The drawing configuration for passive surfaces
  /// @param gridConfig The drawing configuraiton for grid
  /// @param outputDir Directory to write to
  static void drawSurfaceArray(
      IVisualization3D& helper, const SurfaceArray& surfaceArray,
      const GeometryContext& gctx,
      const Transform3& transform = Transform3::Identity(),
      const ViewConfig& sensitiveConfig = s_viewSensitive,
      const ViewConfig& passiveConfig = s_viewPassive,
      const ViewConfig& gridConfig = s_viewGrid,
      const std::string& outputDir = ".");

  /// Helper method to draw AbstractVolume objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param volume The volume to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param transform An option additional transform
  /// @param viewConfig The drawing configuration for boundary surfaces
  static void drawVolume(IVisualization3D& helper, const AbstractVolume& volume,
                         const GeometryContext& gctx,
                         const Transform3& transform = Transform3::Identity(),
                         const ViewConfig& viewConfig = s_viewVolume);

  /// Helper method to draw AbstractVolume objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param layer The tracking layer to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param layerConfig The drawing configuration for passive surfaces
  /// @param sensitiveConfig The drawing configuration for sensitive surfaces
  /// @param gridConfig The drawing configuraiton for grid display
  /// @param outputDir Directory to write to
  static void drawLayer(IVisualization3D& helper, const Layer& layer,
                        const GeometryContext& gctx,
                        const ViewConfig& layerConfig = s_viewPassive,
                        const ViewConfig& sensitiveConfig = s_viewSensitive,
                        const ViewConfig& gridConfig = s_viewGrid,
                        const std::string& outputDir = ".");

  /// Helper method to draw AbstractVolume objects
  ///
  /// @param [in,out] helper The visualization helper
  /// @param tVolume The tracking volume to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param containerView The drawing configuration for a container volume
  /// @param volumeView The drawing configuration for the navigation level
  /// volume
  /// @param layerView The drawing configuration for passive surfaces
  /// @param sensitiveView The drawing configuration for sensitive surfaces
  /// @param gridView The drawing configuraiton for grid display
  /// @param writeIt The prescription to write it or not
  /// @param tag The (optional) additional output tag
  /// @param outputDir Directory to write to
  static void drawTrackingVolume(
      IVisualization3D& helper, const TrackingVolume& tVolume,
      const GeometryContext& gctx,
      const ViewConfig& containerView = s_viewVolume,
      const ViewConfig& volumeView = s_viewVolume,
      const ViewConfig& layerView = s_viewPassive,
      const ViewConfig& sensitiveView = s_viewSensitive,
      const ViewConfig& gridView = s_viewGrid, bool writeIt = true,
      const std::string& tag = "", const std::string& outputDir = ".");

  /// Helper method to draw lines - base for all lines
  ///
  /// @param [in,out] helper The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param arrows [ -1 | 0 | 1 | 2 ] = [ start | none | end | both ]
  /// @param arrowLength wrt halflength
  /// @param arrowWidth wrt thickness
  /// @param viewConfig The drawing configuration for this segment
  static void drawSegmentBase(IVisualization3D& helper, const Vector3& start,
                              const Vector3& end, int arrows = 0,
                              double arrowLength = 0., double arrowWidth = 0.,
                              const ViewConfig& viewConfig = s_viewLine);

  /// Convenience function : line
  ///
  /// @param [in,out] helper The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param viewConfig The drawing configuration for this segment
  static void drawSegment(IVisualization3D& helper, const Vector3& start,
                          const Vector3& end,
                          const ViewConfig& viewConfig = s_viewLine);

  /// Convenience function : arrow pointing back
  ///
  /// @param [in,out] helper The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param arrowLength wrt thickness
  /// @param arrowWidth wrt thickness
  /// @param viewConfig The drawing configuration for this segment
  static void drawArrowBackward(IVisualization3D& helper, const Vector3& start,
                                const Vector3& end, double arrowLength,
                                double arrowWidth,
                                const ViewConfig& viewConfig = s_viewLine);

  /// Convenience function : arrow pointing forwad
  ///
  /// @param [in,out] helper The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param arrowLength wrt thickness
  /// @param arrowWidth wrt thickness
  /// @param viewConfig The drawing configuration for this segment
  static void drawArrowForward(IVisualization3D& helper, const Vector3& start,
                               const Vector3& end, double arrowLength,
                               double arrowWidth,
                               const ViewConfig& viewConfig = s_viewLine);

  /// Convenience function : arrow pointing both directions
  ///
  /// @param [in,out] helper The visualization helper
  /// @param start The start point
  /// @param end The end point
  /// @param arrowLength wrt thickness
  /// @param arrowWidth wrt thickness
  /// @param viewConfig The drawing configuration for this segment
  static void drawArrowsBoth(IVisualization3D& helper, const Vector3& start,
                             const Vector3& end, double arrowLength,
                             double arrowWidth,
                             const ViewConfig& viewConfig = s_viewLine);
};

}  // namespace Acts