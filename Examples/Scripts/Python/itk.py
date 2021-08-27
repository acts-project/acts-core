#!/usr/bin/env python3
import sys
from pathlib import Path

from acts.examples import (
    TGeoDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
    Interval,
)

import acts

from acts import MaterialMapJsonConverter, UnitConstants as u


def runITk(
    trackingGeometry,
    decorators,
    outputDir: Path,
    events=1,
    outputObj=True,
    outputCsv=False,
):

    for ievt in range(events):
        eventStore = WhiteBoard(name=f"EventStore#{ievt}", level=acts.logging.INFO)
        ialg = 0

        context = AlgorithmContext(ialg, ievt, eventStore)

        for cdr in decorators:
            r = cdr.decorate(context)
            if r != ProcessCode.SUCCESS:
                raise RuntimeError("Failed to decorate event context")

        if outputCsv:
            writer = CsvTrackingGeometryWriter(
                level=acts.logging.INFO,
                trackingGeometry=trackingGeometry,
                outputDir=str(outputDir / "csv"),
                writePerEvent=True,
            )
            writer.write(context)

        if outputObj:
            writer = ObjTrackingGeometryWriter(
                level=acts.logging.INFO,
                outputDir=str(outputDir / "obj"),
            )
            writer.write(context, trackingGeometry)


def buildITkGeometry(geo_dir: Path):
    Volume = TGeoDetector.Config.Volume
    LayerTriplet = TGeoDetector.Config.LayerTriplet

    matDeco = acts.IMaterialDecorator.fromFile(
        geo_dir / "atlas/itk-hgtd/material-maps-ITk-HGTD.json",
        level=acts.logging.VERBOSE,
    )

    return TGeoDetector.create(
        fileName=str(geo_dir / "atlas/itk-hgtd/ATLAS-ITk-HGTD.tgeo.root"),
        mdecorator=matDeco,
        buildBeamPipe=True,
        unitScalor=1.0,  # explicit units
        beamPipeRadius=29.0 * u.mm,
        beamPipeHalflengthZ=3000.0 * u.mm,
        beamPipeLayerThickness=0.8 * u.mm,
        surfaceLogLevel=acts.logging.WARNING,
        layerLogLevel=acts.logging.WARNING,
        volumeLogLevel=acts.logging.WARNING,
        volumes=[
            Volume(
                name="InnerPixels",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((0 * u.mm, 135 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -250 * u.mm),
                    central=(-250 * u.mm, 250 * u.mm),
                    positive=(250 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(negative=-1.0, central=5 * u.mm, positive=-1.0),
                splitTolZ=LayerTriplet(
                    negative=10 * u.mm, central=-1.0, positive=10 * u.mm
                ),
            ),
            Volume(
                name="OuterPixels",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("Pixel::Pixel"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["Pixel::siLog"]),
                sensitiveAxes=LayerTriplet("YZX"),
                rRange=LayerTriplet((135 * u.mm, 350 * u.mm)),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -377 * u.mm),
                    central=(-377 * u.mm, 377 * u.mm),
                    positive=(3000 * u.mm, 377 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=15 * u.mm, central=5 * u.mm, positive=15 * u.mm
                ),
                splitTolZ=LayerTriplet(
                    negative=20 * u.mm, central=-1.0, positive=20 * u.mm
                ),
            ),
            Volume(
                name="Strips",
                layers=LayerTriplet(True),
                subVolumeName=LayerTriplet("*"),
                binToleranceR=(5 * u.mm, 5 * u.mm),
                binToleranceZ=(5 * u.mm, 5 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(
                    negative=["SCT::ECSensor*"],
                    central=["SCT::BRLSensor*"],
                    positive=["SCT::ECSensor*"],
                ),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(-1.0, 1050 * u.mm),
                    central=(380 * u.mm, 1050 * u.mm),
                    positive=(-1.0, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-3000 * u.mm, -1400 * u.mm),
                    central=(-1400 * u.mm, 1400 * u.mm),
                    positive=(1400 * u.mm, 3000 * u.mm),
                ),
                splitTolR=LayerTriplet(
                    negative=-1.0,
                    central=35 * u.mm,
                    positive=-1.0,
                ),
                splitTolZ=LayerTriplet(
                    negative=35 * u.mm, central=-1.0, positive=35 * u.mm
                ),
            ),
            Volume(
                name="HGTD",
                layers=LayerTriplet(positive=True, central=False, negative=True),
                subVolumeName=LayerTriplet("HGTD::HGTD"),
                binToleranceR=(15 * u.mm, 15 * u.mm),
                binToleranceZ=(15 * u.mm, 15 * u.mm),
                binTolerancePhi=(0.025 * u.mm, 0.025 * u.mm),
                sensitiveNames=LayerTriplet(["HGTD::HGTDSiSensor*"]),
                sensitiveAxes=LayerTriplet("XYZ"),
                rRange=LayerTriplet(
                    negative=(0 * u.mm, 1050 * u.mm),
                    positive=(0 * u.mm, 1050 * u.mm),
                ),
                zRange=LayerTriplet(
                    negative=(-4000 * u.mm, -3000 * u.mm),
                    postive=(3000 * u.mm, 4000 * u.mm),
                ),
                splitTolR=LayerTriplet(-1.0),
                splitTolZ=LayerTriplet(negative=10 * u.mm, positive=10 * u.mm),
            ),
        ],
    )


if "__main__" == __name__:

    assert len(sys.argv == 1), "Please provide the detector example input directory"
    geo_example_dir = Path(sys.argv[0])
    assert geo_example_dir.exists(), "Detector example input directory missing"

    detector, trackingGeometry, decorators = buildITkGeometry(geo_example_dir)

    runITk(
        trackingGeometry=trackingGeometry,
        decorators=decorators,
        outputDir=Path.cwd(),
        outputCsv=False,
        outputObj=True,
    )
