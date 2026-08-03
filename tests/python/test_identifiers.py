import pytest

from mlsnpdr.identifiers import make_module_uid, parse_module_uid


def test_module_uid_normalizes_legacy_spelling():
    assert (
        make_module_uid("PhysicalPPIN", "louvain", "c3", "m10")
        == "physicalppin__louvain__C3__M10"
    )


def test_module_uid_round_trip():
    identity = parse_module_uid("physicalppin__louvain__C3__M10")
    assert identity.legacy_module_id == "physicalPPIN|Louvain|C3|M10"
    assert identity.uid == "physicalppin__louvain__C3__M10"


def test_invalid_network_fails():
    with pytest.raises(ValueError, match="Unknown network"):
        make_module_uid("unknown", "WF", "C1", "M1")

