import pytest

from helpers import dd4hepEnabled

import acts.examples


def count_surfaces(geo):
    __tracebackhide__ = True
    nSurfaces = 0

    def visit(srf):
        nonlocal nSurfaces
        nSurfaces += 1

    geo.visitSurfaces(visit)

    return nSurfaces


def test_generic_geometry():
    detector, geo, contextDecorators = acts.examples.GenericDetector.create()
    assert detector is not None
    assert geo is not None
    assert contextDecorators is not None

    assert count_surfaces(geo) == 18728


@pytest.mark.skipif(not dd4hepEnabled, reason="DD4hep is not set up")
def test_odd():
    from acts.examples.dd4hep import DD4hepGeometryService, DD4hepDetector

    dd4hepConfig = DD4hepGeometryService.Config(
        xmlFileNames=["thirdparty/OpenDataDetector/xml/OpenDataDetector.xml"]
    )
    detector = DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()
    matDeco = acts.JsonMaterialDecorator(
        rConfig=config,
        jFileName="thirdparty/OpenDataDetector/config/odd-material-mapping-config.json",
        level=acts.logging.ERROR,
    )

    geo, _ = detector.finalize(dd4hepConfig, matDeco)

    assert count_surfaces(geo) == 18824


def test_aligned_detector():
    detector, geo, deco = acts.examples.AlignedDetector.create()

    assert detector is not None
    assert geo is not None
    assert deco is not None

    assert count_surfaces(geo) == 18728
