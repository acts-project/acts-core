// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootMaterialTrackWriter.hpp"

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <ios>
#include <iostream>
#include <stdexcept>

#include <TFile.h>
#include <TTree.h>

using Acts::VectorHelpers::eta;
using Acts::VectorHelpers::perp;
using Acts::VectorHelpers::phi;

ActsExamples::RootMaterialTrackWriter::RootMaterialTrackWriter(
    const ActsExamples::RootMaterialTrackWriter::Config& config,
    Acts::Logging::Level level)
    : WriterT(config.collection, "RootMaterialTrackWriter", level),
      m_cfg(config) {
  // An input collection name and tree name must be specified
  if (m_cfg.collection.empty()) {
    throw std::invalid_argument("Missing input collection");
  } else if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  // Setup ROOT I/O
  m_outputFile = TFile::Open(m_cfg.filePath.c_str(), m_cfg.fileMode.c_str());
  if (m_outputFile == nullptr) {
    throw std::ios_base::failure("Could not open '" + m_cfg.filePath);
  }

  m_outputFile->cd();
  m_outputTree =
      new TTree(m_cfg.treeName.c_str(), "TTree from RootMaterialTrackWriter");
  if (m_outputTree == nullptr)
    throw std::bad_alloc();

  // Set the branches
  m_outputTree->Branch("event_id", &m_eventId);
  m_outputTree->Branch("v_x", &m_v_x);
  m_outputTree->Branch("v_y", &m_v_y);
  m_outputTree->Branch("v_z", &m_v_z);
  m_outputTree->Branch("v_px", &m_v_px);
  m_outputTree->Branch("v_py", &m_v_py);
  m_outputTree->Branch("v_pz", &m_v_pz);
  m_outputTree->Branch("v_phi", &m_v_phi);
  m_outputTree->Branch("v_eta", &m_v_eta);
  m_outputTree->Branch("t_X0", &m_tX0);
  m_outputTree->Branch("t_L0", &m_tL0);
  m_outputTree->Branch("mat_x", &m_step_x);
  m_outputTree->Branch("mat_y", &m_step_y);
  m_outputTree->Branch("mat_z", &m_step_z);
  m_outputTree->Branch("mat_dx", &m_step_dx);
  m_outputTree->Branch("mat_dy", &m_step_dy);
  m_outputTree->Branch("mat_dz", &m_step_dz);
  m_outputTree->Branch("mat_step_length", &m_step_length);
  m_outputTree->Branch("mat_X0", &m_step_X0);
  m_outputTree->Branch("mat_L0", &m_step_L0);
  m_outputTree->Branch("mat_A", &m_step_A);
  m_outputTree->Branch("mat_Z", &m_step_Z);
  m_outputTree->Branch("mat_rho", &m_step_rho);

  if (m_cfg.prePostStep) {
    m_outputTree->Branch("mat_sx", &m_step_sx);
    m_outputTree->Branch("mat_sy", &m_step_sy);
    m_outputTree->Branch("mat_sz", &m_step_sz);
    m_outputTree->Branch("mat_ex", &m_step_ex);
    m_outputTree->Branch("mat_ey", &m_step_ey);
    m_outputTree->Branch("mat_ez", &m_step_ez);
  }
  if (m_cfg.storeSurface) {
    m_outputTree->Branch("sur_id", &m_sur_id);
    m_outputTree->Branch("sur_type", &m_sur_type);
    m_outputTree->Branch("sur_x", &m_sur_x);
    m_outputTree->Branch("sur_y", &m_sur_y);
    m_outputTree->Branch("sur_z", &m_sur_z);
    m_outputTree->Branch("sur_range_min", &m_sur_range_min);
    m_outputTree->Branch("sur_range_max", &m_sur_range_max);
  }
  if (m_cfg.storeVolume) {
    m_outputTree->Branch("vol_id", &m_vol_id);
  }
}

ActsExamples::RootMaterialTrackWriter::~RootMaterialTrackWriter() {}

ActsExamples::ProcessCode ActsExamples::RootMaterialTrackWriter::endRun() {
  // write the tree and close the file
  ACTS_INFO("Writing ROOT output File : " << m_cfg.filePath);
  m_outputFile->cd();
  m_outputTree->Write();
  m_outputFile->Close();
  return ActsExamples::ProcessCode::SUCCESS;
}

ActsExamples::ProcessCode ActsExamples::RootMaterialTrackWriter::writeT(
    const AlgorithmContext& ctx,
    const std::vector<Acts::RecordedMaterialTrack>& materialTracks) {
  // Exclusive access to the tree while writing
  std::lock_guard<std::mutex> lock(m_writeMutex);

  m_eventId = ctx.eventNumber;
  // Loop over the material tracks and write them out
  for (auto& mtrack : materialTracks) {
    // Clearing the vector first
    m_step_sx.clear();
    m_step_sy.clear();
    m_step_sz.clear();
    m_step_x.clear();
    m_step_y.clear();
    m_step_z.clear();
    m_step_ex.clear();
    m_step_ey.clear();
    m_step_ez.clear();
    m_step_dx.clear();
    m_step_dy.clear();
    m_step_dz.clear();
    m_step_length.clear();
    m_step_X0.clear();
    m_step_L0.clear();
    m_step_A.clear();
    m_step_Z.clear();
    m_step_rho.clear();

    m_sur_id.clear();
    m_sur_type.clear();
    m_sur_x.clear();
    m_sur_y.clear();
    m_sur_z.clear();
    m_sur_range_min.clear();
    m_sur_range_max.clear();

    m_vol_id.clear();

    // Reserve the vector then
    size_t mints = mtrack.second.materialInteractions.size();
    m_step_sx.reserve(mints);
    m_step_sy.reserve(mints);
    m_step_sz.reserve(mints);
    m_step_x.reserve(mints);
    m_step_y.reserve(mints);
    m_step_z.reserve(mints);
    m_step_ex.reserve(mints);
    m_step_ey.reserve(mints);
    m_step_ez.reserve(mints);
    m_step_dx.reserve(mints);
    m_step_dy.reserve(mints);
    m_step_dz.reserve(mints);
    m_step_length.reserve(mints);
    m_step_X0.reserve(mints);
    m_step_L0.reserve(mints);
    m_step_A.reserve(mints);
    m_step_Z.reserve(mints);
    m_step_rho.reserve(mints);

    m_sur_id.reserve(mints);
    m_sur_type.reserve(mints);
    m_sur_x.reserve(mints);
    m_sur_y.reserve(mints);
    m_sur_z.reserve(mints);
    m_sur_range_min.reserve(mints);
    m_sur_range_max.reserve(mints);

    m_vol_id.reserve(mints);

    // reset the global counter
    if (m_cfg.recalculateTotals) {
      m_tX0 = 0.;
      m_tL0 = 0.;
    } else {
      m_tX0 = mtrack.second.materialInX0;
      m_tL0 = mtrack.second.materialInL0;
    }

    // set the track information at vertex
    m_v_x = mtrack.first.first.x();
    m_v_y = mtrack.first.first.y();
    m_v_z = mtrack.first.first.z();
    m_v_px = mtrack.first.second.x();
    m_v_py = mtrack.first.second.y();
    m_v_pz = mtrack.first.second.z();
    m_v_phi = phi(mtrack.first.second);
    m_v_eta = eta(mtrack.first.second);

    // an now loop over the material
    for (auto& mint : mtrack.second.materialInteractions) {
      auto direction = mint.direction.normalized();

      // The material step position information
      m_step_x.push_back(mint.position.x());
      m_step_y.push_back(mint.position.y());
      m_step_z.push_back(mint.position.z());
      m_step_dx.push_back(direction.x());
      m_step_dy.push_back(direction.y());
      m_step_dz.push_back(direction.z());

      if (m_cfg.prePostStep) {
        Acts::Vector3 prePos =
            mint.position - 0.5 * mint.pathCorrection * direction;
        Acts::Vector3 posPos =
            mint.position + 0.5 * mint.pathCorrection * direction;

        m_step_sx.push_back(prePos.x());
        m_step_sy.push_back(prePos.y());
        m_step_sz.push_back(prePos.z());
        m_step_ex.push_back(posPos.x());
        m_step_ey.push_back(posPos.y());
        m_step_ez.push_back(posPos.z());
      }

      // Store surface information
      if (m_cfg.storeSurface) {
        const Acts::Surface* surface = mint.surface;
        Acts::GeometryIdentifier slayerID;
        if (surface) {
          auto sfIntersection = surface->intersect(
              ctx.geoContext, mint.position, mint.direction, true);
          slayerID = surface->geometryId();
          m_sur_id.push_back(slayerID.value());
          m_sur_type.push_back(surface->type());
          m_sur_x.push_back(sfIntersection.intersection.position.x());
          m_sur_y.push_back(sfIntersection.intersection.position.y());
          m_sur_z.push_back(sfIntersection.intersection.position.z());

          const Acts::SurfaceBounds& surfaceBounds = surface->bounds();
          const Acts::RadialBounds* radialBounds =
              dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
          const Acts::CylinderBounds* cylinderBounds =
              dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);

          if (radialBounds) {
            m_sur_range_min.push_back(radialBounds->rMin());
            m_sur_range_max.push_back(radialBounds->rMax());
          } else if (cylinderBounds) {
            m_sur_range_min.push_back(
                -cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ));
            m_sur_range_max.push_back(
                cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ));
          } else {
            m_sur_range_min.push_back(0);
            m_sur_range_max.push_back(0);
          }
        } else {
          slayerID.setVolume(0);
          slayerID.setBoundary(0);
          slayerID.setLayer(0);
          slayerID.setApproach(0);
          slayerID.setSensitive(0);
          m_sur_id.push_back(slayerID.value());
          m_sur_type.push_back(-1);

          m_sur_x.push_back(0);
          m_sur_y.push_back(0);
          m_sur_z.push_back(0);
          m_sur_range_min.push_back(0);
          m_sur_range_max.push_back(0);
        }
      }

      // store volume information
      if (m_cfg.storeVolume) {
        const Acts::Volume* volume = mint.volume;
        Acts::GeometryIdentifier vlayerID;
        if (volume) {
          vlayerID = volume->geometryId();
          m_vol_id.push_back(vlayerID.value());
        } else {
          vlayerID.setVolume(0);
          vlayerID.setBoundary(0);
          vlayerID.setLayer(0);
          vlayerID.setApproach(0);
          vlayerID.setSensitive(0);
          m_vol_id.push_back(vlayerID.value());
        }
      }

      // the material information
      const auto& mprops = mint.materialSlab;
      m_step_length.push_back(mprops.thickness());
      m_step_X0.push_back(mprops.material().X0());
      m_step_L0.push_back(mprops.material().L0());
      m_step_A.push_back(mprops.material().Ar());
      m_step_Z.push_back(mprops.material().Z());
      m_step_rho.push_back(mprops.material().massDensity());
      // re-calculate if defined to do so
      if (m_cfg.recalculateTotals) {
        m_tX0 += mprops.thicknessInX0();
        m_tL0 += mprops.thicknessInL0();
      }
    }
    // write to
    m_outputTree->Fill();
  }

  // return success
  return ActsExamples::ProcessCode::SUCCESS;
}
