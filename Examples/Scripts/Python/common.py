import acts
import acts.examples


def getOpenDataDetector(mdecorator=None):
    import acts.examples.dd4hep

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    if mdecorator is None:
        mdecorator = acts.JsonMaterialDecorator(
            rConfig=config,
            jFileName="thirdparty/OpenDataDetector/config/odd-material-mapping-config.json",
            level=acts.logging.ERROR,
        )

    trackingGeometry, deco = detector.finalize(dd4hepConfig, mdecorator)

    return detector, trackingGeometry, deco


def addPythia8(
    sequencer: acts.examples.Sequencer,
    rnd: acts.examples.RandomNumbers,
    nhard=1,
    npileup=200,
    beam0=acts.PdgParticle.eProton,
    beam1=acts.PdgParticle.eProton,
    cmsEnergy=14 * acts.UnitConstants.TeV,
    vertexStddev: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
    vertexMean: acts.Vector4 = acts.Vector4(0, 0, 0, 0),
):
    vertexGenerator = acts.examples.GaussianVertexGenerator(
        stddev=vertexStddev, mean=vertexMean
    )

    generators = []
    if nhard > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=nhard),
                vertex=vertexGenerator,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=acts.logging.INFO,
                    pdgBeam0=beam0,
                    pdgBeam1=beam1,
                    cmsEnergy=cmsEnergy,
                    settings=["HardQCD:all = on"],
                ),
            )
        )
    if npileup > 0:
        generators.append(
            acts.examples.EventGenerator.Generator(
                multiplicity=acts.examples.FixedMultiplicityGenerator(n=npileup),
                vertex=vertexGenerator,
                particles=acts.examples.pythia8.Pythia8Generator(
                    level=acts.logging.INFO,
                    pdgBeam0=beam0,
                    pdgBeam1=beam1,
                    cmsEnergy=cmsEnergy,
                    settings=["SoftQCD:all = on"],
                ),
            )
        )

    # Input
    evGen = acts.examples.EventGenerator(
        level=acts.logging.INFO,
        generators=generators,
        outputParticles="particles_input",
        randomNumbers=rnd,
    )

    sequencer.addReader(evGen)

    return evGen
