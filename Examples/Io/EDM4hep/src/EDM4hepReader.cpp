// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/EDM4hep/EDM4hepReader.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/EDM4hep/EDM4hepUtil.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <algorithm>
#include <iomanip>
#include <map>
#include <stdexcept>

#include <DD4hep/Detector.h>
#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <podio/Frame.h>
#include <podio/ObjectID.h>

namespace ActsExamples {

EDM4hepReader::EDM4hepReader(const Config& config, Acts::Logging::Level level)
    : m_cfg(config),
      m_logger(Acts::getDefaultLogger("EDM4hepParticleReader", level)) {
  if (m_cfg.outputParticlesGenerator.empty()) {
    throw std::invalid_argument(
        "Missing output collection generator particles");
  }
  if (m_cfg.outputParticlesSimulation.empty()) {
    throw std::invalid_argument(
        "Missing output collection simulated particles");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing output collection sim hits");
  }

  m_eventsRange = std::make_pair(0, reader().getEntries("events"));

  m_outputParticlesGenerator.initialize(m_cfg.outputParticlesGenerator);
  m_outputParticlesSimulation.initialize(m_cfg.outputParticlesSimulation);
  m_outputSimHits.initialize(m_cfg.outputSimHits);

  m_cfg.trackingGeometry->visitSurfaces([&](const auto* surface) {
    const auto* detElement = dynamic_cast<const Acts::DD4hepDetectorElement*>(
        surface->associatedDetectorElement());

    if (detElement == nullptr) {
      ACTS_ERROR("Surface has no associated detector element");
      return;
    }

    const auto translation = detElement->sourceElement()
                                 .nominal()
                                 .worldTransformation()
                                 .GetTranslation();
    Acts::Vector3 position;
    position << translation[0], translation[1], translation[2];
    position *= Acts::UnitConstants::cm;

    m_surfaceMap.insert({detElement->sourceElement().key(), surface});
  });
}

Acts::PodioUtil::ROOTReader& EDM4hepReader::reader() {
  bool exists = false;
  auto& reader = m_reader.local(exists);
  if (!exists) {
    reader.openFile(m_cfg.inputPath);
  }

  return reader;
}

std::string EDM4hepReader::name() const {
  return "EDM4hepReader";
}

std::pair<std::size_t, std::size_t> EDM4hepReader::availableEvents() const {
  return m_eventsRange;
}

namespace {
std::string vid(unsigned int vtx) {
  return "V" + std::to_string(vtx);
}

std::string pid(const SimParticle& particle) {
  return "P" + std::to_string(particle.particleId().value());
}

std::string plabel(const SimParticle& particle) {
  using namespace Acts::UnitLiterals;
  std::stringstream ss;
  ss << particle.pdg() << "\\n(" << particle.particleId() << ")\\n"
     << "p=" << std::setprecision(3) << particle.absoluteMomentum() / 1_GeV
     << " GeV";
  return ss.str();
}

}  // namespace

void EDM4hepReader::graphviz(std::ostream& os,
                             const std::vector<SimParticle>& particles,
                             const ParentRelationship& parents) const {
  os << "digraph Event {\n";

  std::set<unsigned int> primaryVertices;

  for (const auto& particle : particles) {
    if (particle.particleId().generation() == 0) {
      primaryVertices.insert(particle.particleId().vertexPrimary());

      os << vid(particle.particleId().vertexPrimary()) << " -> "
         << pid(particle) << ";\n";
    }

    os << pid(particle) << " [label=\"" << plabel(particle) << "\"];\n";
  }

  for (const auto [childIdx, parentIdx] : parents) {
    const auto& child = particles[childIdx];
    const auto& parent = particles[parentIdx];
    os << pid(parent) << " -> " << pid(child);

    if (parent.particleId().vertexSecondary() ==
        child.particleId().vertexSecondary()) {
      os << " [style=dashed]";
    }

    os << ";\n";
  }

  for (unsigned int vtx : primaryVertices) {
    os << vid(vtx) << " [label=\"PV" << vtx << "\"];\n";
  }

  os << "}";
}

ProcessCode EDM4hepReader::read(const AlgorithmContext& ctx) {
  podio::Frame frame = reader().readEntry("events", ctx.eventNumber);
  const auto& mcParticleCollection =
      frame.get<edm4hep::MCParticleCollection>(m_cfg.inputParticles);

  ACTS_DEBUG("Reading EDM4hep inputs");

  std::vector<SimParticle> unorderedParticlesInitial;

  // Read particles from the input file
  // Find particles without parents and group them by vtx position to find
  // primary vertices
  std::vector<std::pair<Acts::Vector3, std::vector<edm4hep::MCParticle>>>
      primaryVertices;
  for (const auto& particle : mcParticleCollection) {
    if (particle.parents_size() > 0) {
      // not a primary vertex
      continue;
    }
    const auto& vtx = particle.getVertex();
    Acts::Vector3 vtxPos = {vtx[0], vtx[1], vtx[2]};
    vtxPos /= Acts::UnitConstants::mm;

    // linear search for vector
    auto it = std::ranges::find_if(primaryVertices, [&vtxPos](const auto& v) {
      return v.first == vtxPos;
    });

    if (it == primaryVertices.end()) {
      ACTS_DEBUG("Found primary vertex at " << vtx.x << ", " << vtx.y << ", "
                                            << vtx.z);
      primaryVertices.push_back({vtxPos, {particle}});
    } else {
      it->second.push_back(particle);
    }
  }

  ACTS_DEBUG("Found " << primaryVertices.size() << " primary vertices");

  // key: child, value: parent
  ParentRelationship parentRelationship;

  // key: input particle index, value: index in the unordered particle
  // container
  std::unordered_map<int, std::size_t> edm4hepParticleMap;

  std::size_t nPrimaryVertices = 0;
  // Walk the particle tree
  for (const auto& [vtxPos, particles] : primaryVertices) {
    nPrimaryVertices += 1;
    ACTS_DEBUG("Walking particle tree for primary vertex at "
               << vtxPos.x() << ", " << vtxPos.y() << ", " << vtxPos.z());
    std::size_t nParticles = 0;
    std::size_t nSecondaryVertices = 0;
    std::size_t maxGen = 0;
    auto startSize = unorderedParticlesInitial.size();
    for (const auto& inParticle : particles) {
      nParticles += 1;
      SimParticle particle =
          EDM4hepUtil::readParticle(inParticle)
              .withParticleId(SimBarcode{}
                                  .setParticle(nParticles)
                                  .setVertexPrimary(nPrimaryVertices));
      ACTS_VERBOSE("+ add particle " << particle);
      ACTS_VERBOSE("  - at " << particle.position().transpose());
      ACTS_VERBOSE("  - createdInSim: " << inParticle.isCreatedInSimulation());
      ACTS_VERBOSE("  - vertexIsNotEndpointOfParent: "
                   << inParticle.vertexIsNotEndpointOfParent());
      ACTS_VERBOSE("  - isStopped: " << inParticle.isStopped());
      ACTS_VERBOSE("  - endpoint: " << inParticle.getEndpoint().x << ", "
                                    << inParticle.getEndpoint().y << ", "
                                    << inParticle.getEndpoint().z);
      const auto pid = particle.particleId();
      unorderedParticlesInitial.push_back(std::move(particle));
      edm4hepParticleMap[inParticle.getObjectID().index] =
          unorderedParticlesInitial.size() - 1;
      processChildren(inParticle, pid, unorderedParticlesInitial,
                      parentRelationship, edm4hepParticleMap,
                      nSecondaryVertices, maxGen);
    }
    ACTS_VERBOSE("Primary vertex complete, produced "
                 << (unorderedParticlesInitial.size() - startSize)
                 << " particles and " << nSecondaryVertices
                 << " secondary vertices in " << maxGen << " generations");
    setSubParticleIds(std::next(unorderedParticlesInitial.begin(), startSize),
                      unorderedParticlesInitial.end());
  }

  ACTS_DEBUG("Found " << unorderedParticlesInitial.size() << " particles");

  std::vector<SimParticle> particlesGeneratorUnordered;
  particlesGeneratorUnordered.reserve(mcParticleCollection.size());
  std::vector<SimParticle> particlesSimulatedUnordered;
  particlesSimulatedUnordered.reserve(mcParticleCollection.size());

  for (const auto& inParticle : mcParticleCollection) {
    auto particleIt = edm4hepParticleMap.find(inParticle.getObjectID().index);
    if (particleIt == edm4hepParticleMap.end()) {
      ACTS_ERROR("Particle " << inParticle.getObjectID().index
                             << " not found in particle map");
      continue;
    }
    const std::size_t index = particleIt->second;
    const auto& particleInitial = unorderedParticlesInitial.at(index);
    if (!inParticle.isCreatedInSimulation()) {
      particlesGeneratorUnordered.push_back(particleInitial);
    }
    SimParticle particleSimulated = particleInitial;

    float time = inParticle.getTime() * Acts::UnitConstants::ns;
    for (const auto& daughter : inParticle.getDaughters()) {
      if (!daughter.vertexIsNotEndpointOfParent()) {
        time = daughter.getTime() * Acts::UnitConstants::ns;
        break;
      }
    }

    particleSimulated.final().setPosition4(
        inParticle.getEndpoint()[0] * Acts::UnitConstants::mm,
        inParticle.getEndpoint()[1] * Acts::UnitConstants::mm,
        inParticle.getEndpoint()[2] * Acts::UnitConstants::mm, time);

    Acts::Vector3 momentumFinal = {inParticle.getMomentumAtEndpoint()[0],
                                   inParticle.getMomentumAtEndpoint()[1],
                                   inParticle.getMomentumAtEndpoint()[2]};
    particleSimulated.final().setDirection(momentumFinal.normalized());
    particleSimulated.final().setAbsoluteMomentum(momentumFinal.norm());

    ACTS_VERBOSE("- Updated particle initial -> final, position: "
                 << particleInitial.fourPosition().transpose() << " -> "
                 << particleSimulated.final().fourPosition().transpose());
    ACTS_VERBOSE("                                     momentum: "
                 << particleInitial.fourMomentum().transpose() << " -> "
                 << particleSimulated.final().fourMomentum().transpose());

    particlesSimulatedUnordered.push_back(particleSimulated);
  }

  std::ranges::sort(particlesGeneratorUnordered, detail::CompareParticleId{});
  std::ranges::sort(particlesSimulatedUnordered, detail::CompareParticleId{});

  SimParticleContainer particlesGenerator{particlesGeneratorUnordered.begin(),
                                          particlesGeneratorUnordered.end()};
  SimParticleContainer particlesSimulated{particlesSimulatedUnordered.begin(),
                                          particlesSimulatedUnordered.end()};

  if (!m_cfg.graphvizOutput.empty()) {
    std::string path = perEventFilepath(m_cfg.graphvizOutput, "particles.dot",
                                        ctx.eventNumber);
    std::ofstream dot(path);
    graphviz(dot, unorderedParticlesInitial, parentRelationship);
  }

  std::vector<SimHit> simHitsUnordered;

  ACTS_DEBUG("Reading sim hits from " << m_cfg.inputSimHits.size()
                                      << " sim hit collections");
  for (const auto& name : m_cfg.inputSimHits) {
    const auto& inputHits = frame.get<edm4hep::SimTrackerHitCollection>(name);

    simHitsUnordered.reserve(simHitsUnordered.size() + inputHits.size());

    for (const auto& hit : inputHits) {
      auto simHit = EDM4hepUtil::readSimHit(
          hit,
          [&](const auto& inParticle) {
            ACTS_VERBOSE(
                "SimHit has source particle: "
                << Acts::EDM4hepUtil::getParticle(hit).getObjectID().index);
            auto it = edm4hepParticleMap.find(inParticle.getObjectID().index);
            if (it == edm4hepParticleMap.end()) {
              ACTS_ERROR(
                  "SimHit has source particle that we did not see before");
              return SimBarcode{};
            }
            const auto& particle = unorderedParticlesInitial.at(it->second);
            ACTS_VERBOSE("- " << inParticle.getObjectID().index << " -> "
                              << particle.particleId());
            return particle.particleId();
          },
          [&](std::uint64_t cellId) {
            ACTS_VERBOSE("CellID: " << cellId);

            const auto& vm =
                m_cfg.dd4hepDetector->dd4hepDetector().volumeManager();

            const auto detElement = vm.lookupDetElement(cellId);

            ACTS_VERBOSE(" -> detElement: " << detElement.name());
            ACTS_VERBOSE("   -> id: " << detElement.id());
            ACTS_VERBOSE("   -> key: " << detElement.key());

            Acts::Vector3 position;
            position << detElement.nominal()
                            .worldTransformation()
                            .GetTranslation()[0],
                detElement.nominal().worldTransformation().GetTranslation()[1],
                detElement.nominal().worldTransformation().GetTranslation()[2];
            position *= Acts::UnitConstants::cm;

            ACTS_VERBOSE("   -> detElement position: " << position.transpose());

            auto it = m_surfaceMap.find(detElement.key());
            if (it == m_surfaceMap.end()) {
              ACTS_ERROR("Unable to find surface for detElement "
                         << detElement.name() << " with cellId " << cellId);
              throw std::runtime_error("Unable to find surface for detElement");
            }
            const auto* surface = it->second;
            if (surface == nullptr) {
              ACTS_ERROR("Unable to find surface for detElement "
                         << detElement.name() << " with cellId " << cellId);
              throw std::runtime_error("Unable to find surface for detElement");
            }
            ACTS_VERBOSE("   -> surface: " << surface->geometryId());
            return surface->geometryId();
          });

      simHitsUnordered.push_back(std::move(simHit));
    }
  }

  std::ranges::sort(simHitsUnordered, detail::CompareGeometryId{});

  SimHitContainer simHits{simHitsUnordered.begin(), simHitsUnordered.end()};

  if (m_cfg.sortSimHitsInTime) {
    ACTS_DEBUG("Sorting sim hits in time");
    std::multimap<ActsFatras::Barcode, std::size_t> hitsByParticle;

    for (std::size_t i = 0; i < simHits.size(); ++i) {
      hitsByParticle.insert({simHits.nth(i)->particleId(), i});
    }

    for (auto it = hitsByParticle.begin(), end = hitsByParticle.end();
         it != end; it = hitsByParticle.upper_bound(it->first)) {
      ACTS_DEBUG("Particle " << it->first << " has "
                             << hitsByParticle.count(it->first) << " hits");

      std::vector<std::size_t> hitIndices;
      hitIndices.reserve(hitsByParticle.count(it->first));
      for (auto hitIndex : makeRange(hitsByParticle.equal_range(it->first))) {
        hitIndices.push_back(hitIndex.second);
      }

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        ACTS_VERBOSE("Before sorting:");
        for (const auto& hitIdx : hitIndices) {
          ACTS_VERBOSE(" - " << hitIdx << " / " << simHits.nth(hitIdx)->index()
                             << " " << simHits.nth(hitIdx)->time());
        }
      }

      std::ranges::sort(hitIndices, {}, [&simHits](std::size_t h) {
        return simHits.nth(h)->time();
      });

      for (std::size_t i = 0; i < hitIndices.size(); ++i) {
        auto& hit = *simHits.nth(hitIndices[i]);
        SimHit updatedHit{hit.geometryId(),     hit.particleId(),
                          hit.fourPosition(),   hit.momentum4Before(),
                          hit.momentum4After(), static_cast<std::int32_t>(i)};
        hit = updatedHit;
      }

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        ACTS_VERBOSE("After sorting:");
        for (const auto& hitIdx : hitIndices) {
          ACTS_VERBOSE(" - " << hitIdx << " / " << simHits.nth(hitIdx)->index()
                             << " " << simHits.nth(hitIdx)->time());
        }
      }
    }
  }

  m_outputParticlesGenerator(ctx, std::move(particlesGenerator));
  m_outputParticlesSimulation(ctx, std::move(particlesSimulated));

  m_outputSimHits(ctx, std::move(simHits));

  return ProcessCode::SUCCESS;
}

void EDM4hepReader::processChildren(
    const edm4hep::MCParticle& inParticle, SimBarcode parentId,
    std::vector<SimParticle>& particles, ParentRelationship& parentRelationship,
    std::unordered_map<int, std::size_t>& particleMap,
    std::size_t& nSecondaryVertices, std::size_t& maxGen) const {
  constexpr auto indent = [&](std::size_t n) {
    std::string result;
    for (std::size_t i = 0; i < n; ++i) {
      result += "  ";
    }
    return result;
  };

  const std::size_t gen = parentId.generation();
  maxGen = std::max(maxGen, gen);

  ACTS_VERBOSE(indent(gen) << "  - processing daughters for input particle "
                           << inParticle.id());
  ACTS_VERBOSE(indent(gen) << "    -> found " << inParticle.daughters_size()
                           << " daughter(s)");

  bool parentDecayed =
      std::any_of(inParticle.daughters_begin(), inParticle.daughters_end(),
                  [](const edm4hep::MCParticle& daughter) {
                    return !daughter.vertexIsNotEndpointOfParent();
                  });
  std::size_t secondaryVertex = 0;
  if (parentDecayed) {
    ACTS_VERBOSE(indent(gen) << "    -> parent decays");
    secondaryVertex = ++nSecondaryVertices;
  }

  std::size_t parentIndex = particles.size() - 1;

  std::size_t nParticles = 0;
  for (const auto& daughter : inParticle.getDaughters()) {
    SimParticle particle = EDM4hepUtil::readParticle(daughter);

    auto pid = parentId.makeDescendant(nParticles++);
    if (daughter.vertexIsNotEndpointOfParent()) {
      // incoming particle survived, interaction via descendant
    } else {
      // incoming particle decayed
      pid = pid.setVertexSecondary(secondaryVertex);
    }
    particle.setParticleId(pid);

    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "+ add particle " << particle);
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - generation: " << particle.particleId().generation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "  - at " << particle.position().transpose());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - createdInSim: "
                 << daughter.isCreatedInSimulation());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - vertexIsNotEndpointOfParent: "
                 << daughter.vertexIsNotEndpointOfParent());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - isStopped: " << daughter.isStopped());
    ACTS_VERBOSE(indent(particle.particleId().generation())
                 << "     - endpoint: " << daughter.getEndpoint().x << ", "
                 << daughter.getEndpoint().y << ", "
                 << daughter.getEndpoint().z);

    particles.push_back(std::move(particle));
    particleMap[daughter.getObjectID().index] = particles.size() - 1;
    parentRelationship[particles.size() - 1] = parentIndex;
    processChildren(daughter, pid, particles, parentRelationship, particleMap,
                    nSecondaryVertices, maxGen);
  }
}

void EDM4hepReader::setSubParticleIds(std::vector<SimParticle>::iterator begin,
                                      std::vector<SimParticle>::iterator end) {
  std::vector<std::size_t> numByGeneration;
  numByGeneration.reserve(10);

  for (auto it = begin; it != end; ++it) {
    auto& particle = *it;
    const auto pid = particle.particleId();
    if (pid.generation() >= numByGeneration.size()) {
      numByGeneration.resize(pid.generation() + 1, 0);
    }
    unsigned int nextSubParticle = numByGeneration[pid.generation()]++;

    auto newPid = particle.particleId().setSubParticle(nextSubParticle);
    particle.setParticleId(newPid);
  }
}

}  // namespace ActsExamples
