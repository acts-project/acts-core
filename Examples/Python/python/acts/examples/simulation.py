from typing import Optional, Union, Any, List
from pathlib import Path
from collections import namedtuple
from collections.abc import Iterable

import acts
from acts.examples import (
    RandomNumbers,
    EventGenerator,
    FixedMultiplicityGenerator,
    CsvParticleWriter,
    ParticlesPrinter,
    RootParticleWriter,
    RootVertexWriter,
)

# Defaults (given as `None` here) use class defaults defined in
# Examples/Algorithms/Generators/ActsExamples/Generators/ParametricParticleGenerator.hpp
MomentumConfig = namedtuple(
    "MomentumConfig",
    ["min", "max", "transverse", "logUniform"],
    defaults=[None, None, None, None],
)
EtaConfig = namedtuple(
    "EtaConfig", ["min", "max", "uniform"], defaults=[None, None, None]
)
PhiConfig = namedtuple("PhiConfig", ["min", "max"], defaults=[None, None])
ParticleConfig = namedtuple(
    "ParticleConfig",
    ["num", "pdg", "randomizeCharge", "charge", "mass"],
    defaults=[None, None, None, None, None],
)
ParticleSelectorConfig = namedtuple(
    "ParticleSelectorConfig",
    [
        "rho",  # (min,max)
        "absZ",  # (min,max)
        "time",  # (min,max)
        "phi",  # (min,max)
        "eta",  # (min,max)
        "absEta",  # (min,max)
        "pt",  # (min,max)
        "m",  # (min,max)
        "hits",  # (min,max)
        "measurements",  # (min,max)
        "removeCharged",  # bool
        "removeNeutral",  # bool
        "removeSecondaries",  # bool
    ],
    defaults=[(None, None)] * 10 + [None] * 3,
)


def _getParticleSelectionKWargs(config: ParticleSelectorConfig) -> dict:
    return {
        "rhoMin": config.rho[0],
        "rhoMax": config.rho[1],
        "absZMin": config.absZ[0],
        "absZMax": config.absZ[1],
        "timeMin": config.time[0],
        "timeMax": config.time[1],
        "phiMin": config.phi[0],
        "phiMax": config.phi[1],
        "etaMin": config.eta[0],
        "etaMax": config.eta[1],
        "absEtaMin": config.absEta[0],
        "absEtaMax": config.absEta[1],
        "ptMin": config.pt[0],
        "ptMax": config.pt[1],
        "mMin": config.m[0],
        "mMax": config.m[1],
        "hitsMin": config.hits[0],
        "hitsMax": config.hits[1],
        "measurementsMin": config.measurements[0],
        "measurementsMax": config.measurements[1],
        "removeCharged": config.removeCharged,
        "removeNeutral": config.removeNeutral,
        "removeSecondaries": config.removeSecondaries,
    }


@acts.examples.NamedTypeArgs(
    momentumConfig=MomentumConfig,
    etaConfig=EtaConfig,
    phiConfig=PhiConfig,
    particleConfig=ParticleConfig,
)
def addParticleGun(
    s: acts.examples.Sequencer,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    momentumConfig: MomentumConfig = MomentumConfig(),
    etaConfig: EtaConfig = EtaConfig(),
    phiConfig: PhiConfig = PhiConfig(),
    particleConfig: ParticleConfig = ParticleConfig(),
    multiplicity: int = 1,
    vtxGen: Optional[EventGenerator.VertexGenerator] = None,
    printParticles: bool = False,
    rnd: Optional[RandomNumbers] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the particle generation using the particle gun

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    momentumConfig : MomentumConfig(min, max, transverse, logUniform)
        momentum configuration: minimum momentum, maximum momentum, transverse, log-uniform
    etaConfig : EtaConfig(min, max, uniform)
        pseudorapidity configuration: eta min, eta max, uniform
    phiConfig : PhiConfig(min, max)
        azimuthal angle configuration: phi min, phi max
    particleConfig : ParticleConfig(num, pdg, randomizeCharge, charge, mass)
        particle configuration: number of particles, particle type, charge flip
    multiplicity : int, 1
        number of generated vertices
    vtxGen : VertexGenerator, None
        vertex generator module
    printParticles : bool, False
        print generated particles
    rnd : RandomNumbers, None
        random number generator
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    rnd = rnd or RandomNumbers(seed=228)

    evGen = EventGenerator(
        level=customLogLevel(),
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=multiplicity),
                vertex=vtxGen
                or acts.examples.GaussianVertexGenerator(
                    mean=acts.Vector4(0, 0, 0, 0),
                    stddev=acts.Vector4(0, 0, 0, 0),
                ),
                particles=acts.examples.ParametricParticleGenerator(
                    **acts.examples.defaultKWArgs(
                        p=(momentumConfig.min, momentumConfig.max),
                        pTransverse=momentumConfig.transverse,
                        pLogUniform=momentumConfig.logUniform,
                        eta=(etaConfig.min, etaConfig.max),
                        phi=(phiConfig.min, phiConfig.max),
                        etaUniform=etaConfig.uniform,
                        numParticles=particleConfig.num,
                        pdg=particleConfig.pdg,
                        randomizeCharge=particleConfig.randomizeCharge,
                        charge=particleConfig.charge,
                        mass=particleConfig.mass,
                    )
                ),
            )
        ],
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        randomNumbers=rnd,
    )
    s.addReader(evGen)

    s.addWhiteboardAlias("particles", evGen.config.outputParticles)
    s.addWhiteboardAlias("vertices_truth", evGen.config.outputVertices)

    s.addWhiteboardAlias("particles_generated_selected", evGen.config.outputParticles)

    if printParticles:
        s.addAlgorithm(
            ParticlesPrinter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            CsvParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            RootParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "particles.root"),
            )
        )

        s.addWriter(
            RootVertexWriter(
                level=customLogLevel(),
                inputVertices=evGen.config.outputVertices,
                filePath=str(outputDirRoot / "vertices.root"),
            )
        )

    return s


def addPythia8(
    s: acts.examples.Sequencer,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    nhard: int = 1,
    npileup: int = 200,
    beam: Optional[
        Union[acts.PdgParticle, Iterable]
    ] = None,  # default: acts.PdgParticle.eProton
    cmsEnergy: Optional[float] = None,  # default: 14 * acts.UnitConstants.TeV
    hardProcess: Optional[Iterable] = None,  # default: ["HardQCD:all = on"]
    pileupProcess: Iterable = ["SoftQCD:all = on"],
    vtxGen: Optional[EventGenerator.VertexGenerator] = None,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    printParticles: bool = False,
    printPythiaEventListing: Optional[Union[None, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the particle generation using Pythia8

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the particle gun steps (returned from addParticleGun)
    rnd : RandomNumbers, None
        random number generator
    nhard, npileup : int, 1, 200
        Number of hard-scatter and pileup vertices
    beam : PdgParticle|[PdgParticle,PdgParticle], eProton
        beam particle(s)
    cmsEnergy : float, 14 TeV
        CMS energy
    hardProcess, pileupProcess : [str], ["HardQCD:all = on"], ["SoftQCD:all = on"]
        hard and pileup processes
    vtxGen : VertexGenerator, None
        vertex generator module
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    printParticles : bool, False
        print generated particles
    printPythiaEventListing
        None or "short" or "long"
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    rnd = rnd or acts.examples.RandomNumbers()

    vtxGen = vtxGen or acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 0, 0), mean=acts.Vector4(0, 0, 0, 0)
    )

    if not isinstance(beam, Iterable):
        beam = (beam, beam)

    if printPythiaEventListing is None:
        printShortEventListing = False
        printLongEventListing = False
    elif printPythiaEventListing == "short":
        printShortEventListing = True
        printLongEventListing = False
    elif printPythiaEventListing == "long":
        printShortEventListing = False
        printLongEventListing = True
    else:
        raise RuntimeError("Invalid pythia config")

    generators = []
    if nhard is not None and nhard > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=customLogLevel(),
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=hardProcess,
                        printLongEventListing=printLongEventListing,
                        printShortEventListing=printShortEventListing,
                    ),
                ),
            )
        )
    if npileup > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
                vertex=vtxGen,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=customLogLevel(),
                    **acts.examples.defaultKWArgs(
                        pdgBeam0=beam[0],
                        pdgBeam1=beam[1],
                        cmsEnergy=cmsEnergy,
                        settings=pileupProcess,
                    ),
                ),
            )
        )

    evGen = acts.examples.EventGenerator(
        level=customLogLevel(),
        generators=generators,
        outputParticles="particles_generated",
        outputVertices="vertices_generated",
        randomNumbers=rnd,
    )
    s.addReader(evGen)

    s.addWhiteboardAlias("particles", evGen.config.outputParticles)
    s.addWhiteboardAlias("vertices_truth", evGen.config.outputVertices)

    s.addWhiteboardAlias("particles_generated_selected", evGen.config.outputParticles)

    if printParticles:
        s.addAlgorithm(
            acts.examples.ParticlesPrinter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
            )
        )

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()

        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                outputDir=str(outputDirCsv),
                outputStem="particles",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()

        s.addWriter(
            acts.examples.RootParticleWriter(
                level=customLogLevel(),
                inputParticles=evGen.config.outputParticles,
                filePath=str(outputDirRoot / "particles.root"),
            )
        )

        s.addWriter(
            acts.examples.RootVertexWriter(
                level=customLogLevel(),
                inputVertices=evGen.config.outputVertices,
                filePath=str(outputDirRoot / "vertices.root"),
            )
        )

    return s


def addGenParticleSelection(
    s: acts.examples.Sequencer,
    config: ParticleSelectorConfig,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection after generation.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    config: ParticleSelectorConfig
        the particle selection configuration
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    selector = acts.examples.ParticleSelector(
        **acts.examples.defaultKWArgs(**_getParticleSelectionKWargs(config)),
        level=customLogLevel(),
        inputParticles="particles_generated",
        outputParticles="tmp_particles_generated_selected",
    )
    s.addAlgorithm(selector)

    s.addWhiteboardAlias("particles_selected", selector.config.outputParticles)
    s.addWhiteboardAlias(
        "particles_generated_selected", selector.config.outputParticles
    )


def addFatras(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    enableInteractions: bool = True,
    pMin: Optional[float] = None,
    inputParticles: str = "particles_generated_selected",
    outputParticles: str = "particles_simulated",
    outputSimHits: str = "simhits",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirObj: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """This function steers the detector simulation using Fatras

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Fatras steps (returned from addFatras)
    trackingGeometry : tracking geometry
    field : magnetic field
    rnd : RandomNumbers
        random number generator
    enableInteractions : Enable the particle interactions in the simulation
    pMin : Minimum monmentum of particles simulated by FATRAS
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    outputDirObj : Path|str, path, None
        the output folder for the Obj output, None triggers no output
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    alg = acts.examples.FatrasSimulation(
        **acts.examples.defaultKWArgs(
            level=customLogLevel(),
            inputParticles=inputParticles,
            outputParticles=outputParticles,
            outputSimHits=outputSimHits,
            randomNumbers=rnd,
            trackingGeometry=trackingGeometry,
            magneticField=field,
            generateHitsOnSensitive=True,
            emScattering=enableInteractions,
            emEnergyLossIonisation=enableInteractions,
            emEnergyLossRadiation=enableInteractions,
            emPhotonConversion=enableInteractions,
            pMin=pMin,
        )
    )
    s.addAlgorithm(alg)

    s.addWhiteboardAlias("particles", outputParticles)

    s.addWhiteboardAlias("particles_simulated_selected", outputParticles)

    addSimWriters(
        s,
        alg.config.outputSimHits,
        outputParticles,
        outputDirCsv,
        outputDirRoot,
        outputDirObj,
        logLevel,
    )

    return s


def addSimWriters(
    s: acts.examples.Sequencer,
    simHits: str = "simhits",
    particlesSimulated: str = "particles_simulated",
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirObj: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            acts.examples.CsvParticleWriter(
                level=customLogLevel(),
                outputDir=str(outputDirCsv),
                inputParticles=particlesSimulated,
                outputStem="particles_simulated",
            )
        )
        s.addWriter(
            acts.examples.CsvSimHitWriter(
                level=customLogLevel(),
                inputSimHits=simHits,
                outputDir=str(outputDirCsv),
                outputStem="hits",
            )
        )

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        s.addWriter(
            acts.examples.RootParticleWriter(
                level=customLogLevel(),
                inputParticles=particlesSimulated,
                filePath=str(outputDirRoot / "particles_simulation.root"),
            )
        )
        s.addWriter(
            acts.examples.RootSimHitWriter(
                level=customLogLevel(),
                inputSimHits=simHits,
                filePath=str(outputDirRoot / "hits.root"),
            )
        )

    if outputDirObj is not None:
        outputDirObj = Path(outputDirObj)
        if not outputDirObj.exists():
            outputDirObj.mkdir()
        s.addWriter(
            acts.examples.ObjSimHitWriter(
                level=customLogLevel(),
                inputSimHits=simHits,
                outputDir=str(outputDirObj),
                outputStem="hits",
            )
        )


# holds the Geant4Handle for potential reuse
__geant4Handle = None


def addGeant4(
    s: acts.examples.Sequencer,
    detector: Optional[Any],
    trackingGeometry: Union[acts.TrackingGeometry, acts.Detector],
    field: acts.MagneticFieldProvider,
    rnd: acts.examples.RandomNumbers,
    volumeMappings: List[str] = [],
    materialMappings: List[str] = ["Silicon"],
    inputParticles: str = "particles_generated_selected",
    outputParticles: str = "particles_simulated",
    outputSimHits: str = "simhits",
    recordHitsOfSecondaries=True,
    keepParticlesWithoutHits=True,
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    outputDirObj: Optional[Union[Path, str]] = None,
    logLevel: Optional[acts.logging.Level] = None,
    killVolume: Optional[acts.Volume] = None,
    killAfterTime: float = float("inf"),
    killSecondaries: bool = False,
    physicsList: str = "FTFP_BERT",
    regionList: List[Any] = [],
) -> None:
    """This function steers the detector simulation using Geant4

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Geant4 steps (returned from addGeant4)
    trackingGeometry : tracking geometry or detector
    field : magnetic field
    rnd : RandomNumbers, None
        random number generator
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    outputDirObj : Path|str, path, None
        the output folder for the Obj output, None triggers no output
    killVolume: acts.Volume, None
        if given, particles are killed when going outside this volume.
    killAfterTime: float
        if given, particle are killed after the global time since event creation exceeds the given value
    killSecondaries: bool
        if given, secondary particles are removed from simulation
    """

    from acts.examples.geant4 import Geant4Simulation, SensitiveSurfaceMapper

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    global __geant4Handle

    smmConfig = SensitiveSurfaceMapper.Config()
    smmConfig.volumeMappings = volumeMappings
    smmConfig.materialMappings = materialMappings
    sensitiveMapper = SensitiveSurfaceMapper.create(
        smmConfig, customLogLevel(), trackingGeometry
    )

    alg = Geant4Simulation(
        level=customLogLevel(),
        geant4Handle=__geant4Handle,
        detector=detector,
        randomNumbers=rnd,
        inputParticles=inputParticles,
        outputParticles=outputParticles,
        outputSimHits=outputSimHits,
        sensitiveSurfaceMapper=sensitiveMapper,
        magneticField=field,
        physicsList=physicsList,
        killVolume=killVolume,
        killAfterTime=killAfterTime,
        killSecondaries=killSecondaries,
        recordHitsOfCharged=True,
        recordHitsOfNeutrals=False,
        recordHitsOfPrimaries=True,
        recordHitsOfSecondaries=recordHitsOfSecondaries,
        recordPropagationSummaries=False,
        keepParticlesWithoutHits=keepParticlesWithoutHits,
    )
    __geant4Handle = alg.geant4Handle
    s.addAlgorithm(alg)

    s.addWhiteboardAlias("particles", outputParticles)

    s.addWhiteboardAlias("particles_simulated_selected", outputParticles)

    addSimWriters(
        s,
        alg.config.outputSimHits,
        outputParticles,
        outputDirCsv,
        outputDirRoot,
        outputDirObj,
        logLevel=logLevel,
    )

    return s


def addSimParticleSelection(
    s: acts.examples.Sequencer,
    config: ParticleSelectorConfig,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection after simulation.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    config: ParticleSelectorConfig
        the particle selection configuration
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    selector = acts.examples.ParticleSelector(
        **acts.examples.defaultKWArgs(**_getParticleSelectionKWargs(config)),
        level=customLogLevel(),
        inputParticles="particles_simulated",
        outputParticles="tmp_particles_simulated_selected",
    )
    s.addAlgorithm(selector)

    s.addWhiteboardAlias("particles_selected", selector.config.outputParticles)
    s.addWhiteboardAlias(
        "particles_simulated_selected", selector.config.outputParticles
    )


def addDigitization(
    s: acts.examples.Sequencer,
    trackingGeometry: Union[acts.TrackingGeometry, acts.Detector],
    field: acts.MagneticFieldProvider,
    digiConfigFile: Union[Path, str],
    outputDirCsv: Optional[Union[Path, str]] = None,
    outputDirRoot: Optional[Union[Path, str]] = None,
    rnd: Optional[acts.examples.RandomNumbers] = None,
    doMerge: Optional[bool] = None,
    minEnergyDeposit: Optional[float] = None,
    logLevel: Optional[acts.logging.Level] = None,
) -> acts.examples.Sequencer:
    """This function steers the digitization step

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the Digitization steps (returned from addDigitization)
    trackingGeometry : tracking geometry or detector
    field : magnetic field
    digiConfigFile : Path|str, path
        Configuration (.json) file for digitization or smearing description
    outputDirCsv : Path|str, path, None
        the output folder for the Csv output, None triggers no output
    outputDirRoot : Path|str, path, None
        the output folder for the Root output, None triggers no output
    rnd : RandomNumbers, None
        random number generator
    """

    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    rnd = rnd or acts.examples.RandomNumbers()

    digiCfg = acts.examples.DigitizationAlgorithm.Config(
        digitizationConfigs=acts.examples.readDigiConfigFromJson(
            str(digiConfigFile),
        ),
        surfaceByIdentifier=trackingGeometry.geoIdSurfaceMap(),
        randomNumbers=rnd,
        inputSimHits="simhits",
        outputMeasurements="measurements",
        outputMeasurementParticlesMap="measurement_particles_map",
        outputMeasurementSimHitsMap="measurement_simhits_map",
        outputParticleMeasurementsMap="particle_measurements_map",
        outputSimHitMeasurementsMap="simhit_measurements_map",
        **acts.examples.defaultKWArgs(
            doMerge=doMerge,
        ),
    )

    # Not sure how to do this in our style
    if minEnergyDeposit is not None:
        digiCfg.minEnergyDeposit = minEnergyDeposit

    digiAlg = acts.examples.DigitizationAlgorithm(digiCfg, customLogLevel())

    s.addAlgorithm(digiAlg)

    if outputDirRoot is not None:
        outputDirRoot = Path(outputDirRoot)
        if not outputDirRoot.exists():
            outputDirRoot.mkdir()
        rmwConfig = acts.examples.RootMeasurementWriter.Config(
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputSimHits=digiAlg.config.inputSimHits,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            filePath=str(outputDirRoot / f"{digiAlg.config.outputMeasurements}.root"),
            surfaceByIdentifier=trackingGeometry.geoIdSurfaceMap(),
        )
        s.addWriter(acts.examples.RootMeasurementWriter(rmwConfig, customLogLevel()))

    if outputDirCsv is not None:
        outputDirCsv = Path(outputDirCsv)
        if not outputDirCsv.exists():
            outputDirCsv.mkdir()
        s.addWriter(
            acts.examples.CsvMeasurementWriter(
                level=customLogLevel(),
                inputMeasurements=digiAlg.config.outputMeasurements,
                inputClusters=digiAlg.config.outputClusters,
                inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
                outputDir=str(outputDirCsv),
            )
        )

    return s


def addDigiParticleSelection(
    s: acts.examples.Sequencer,
    config: ParticleSelectorConfig,
    logLevel: Optional[acts.logging.Level] = None,
) -> None:
    """
    This function steers the particle selection after digitization.

    Parameters
    ----------
    s: Sequencer
        the sequencer module to which we add the ParticleSelector
    config: ParticleSelectorConfig
        the particle selection configuration
    """
    customLogLevel = acts.examples.defaultLogging(s, logLevel)

    selector = acts.examples.ParticleSelector(
        **acts.examples.defaultKWArgs(**_getParticleSelectionKWargs(config)),
        level=customLogLevel(),
        inputParticles="particles_simulated_selected",
        inputParticleMeasurementsMap="particle_measurements_map",
        outputParticles="tmp_particles_digitized_selected",
    )
    s.addAlgorithm(selector)

    s.addWhiteboardAlias("particles_selected", selector.config.outputParticles)
    s.addWhiteboardAlias(
        "particles_digitized_selected", selector.config.outputParticles
    )
