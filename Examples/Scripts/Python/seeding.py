#!/usr/bin/env python3
from pathlib import Path
from enum import Enum
import argparse

import acts
import acts.examples
from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants


# Graciously taken from https://stackoverflow.com/a/60750535/4280680
class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums
    """

    def __init__(self, **kwargs):
        # Pop off the type value
        enum_type = kwargs.pop("enum", None)

        # Ensure an Enum subclass is provided
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, Enum):
            raise TypeError("type must be an Enum when using EnumAction")

        # Generate choices from the Enum
        kwargs.setdefault("choices", tuple(e.name for e in enum_type))

        super(EnumAction, self).__init__(**kwargs)

        self._enum = enum_type

    def __call__(self, parser, namespace, values, option_string=None):
        for e in self._enum:
            if e.name == values:
                setattr(namespace, self.dest, e)
                break
        else:
            raise ValueError("%s is not a validly enumerated algorithm." % values)


from acts.examples.reconstruction import SeedingAlgorithm


def runSeeding(
    trackingGeometry,
    field,
    outputDir,
    s=None,
    seedingAlgorithm=SeedingAlgorithm.Default,
):
    from acts.examples.simulation import (
        addParticleGun,
        EtaConfig,
        PhiConfig,
        ParticleConfig,
        addFatras,
        addDigitization,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )
    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    addParticleGun(
        s,
        EtaConfig(-2.0, 2.0),
        ParticleConfig(4, acts.PdgParticle.eMuon, True),
        PhiConfig(0.0, 360.0 * u.degree),
        multiplicity=2,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )

    srcdir = Path(__file__).resolve().parent.parent.parent.parent
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(1.0 * u.GeV, None),
            eta=(-2.5, 2.5),
            measurements=(9, None),
            removeNeutral=True,
        ),
    )

    from acts.examples.reconstruction import (
        addSeeding,
        SeedFinderConfigArg,
        SeedFinderOptionsArg,
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        SeedFinderConfigArg(
            r=(None, 200 * u.mm),  # rMin=default, 33mm
            deltaR=(1 * u.mm, 60 * u.mm),
            collisionRegion=(-250 * u.mm, 250 * u.mm),
            z=(-2000 * u.mm, 2000 * u.mm),
            maxSeedsPerSpM=1,
            sigmaScattering=50,
            radLengthPerSeed=0.1,
            minPt=500 * u.MeV,
            impactMax=3 * u.mm,
        ),
        SeedFinderOptionsArg(
            bFieldInZ=2 * u.T,
        ),
        acts.logging.VERBOSE,
        seedingAlgorithm=seedingAlgorithm,
        geoSelectionConfigFile=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        outputDirRoot=outputDir,
    )
    return s


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Example script to run seed finding",
    )

    p.add_argument(
        "--algorithm",
        action=EnumAction,
        enum=SeedingAlgorithm,
        default=SeedingAlgorithm.Default,
        help="Select the seeding algorithm to use",
    )

    args = p.parse_args()
    # detector = getOpenDataDetector()
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runSeeding(
        trackingGeometry, field, outputDir=Path.cwd(), seedingAlgorithm=args.algorithm
    ).run()
