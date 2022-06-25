import pytest
from ..pyfastg import validate_seq_and_compute_gc as f


def test_basic():
    assert f("ACGTA") == 2 / 5
    assert f("AAATATATAAATTTT") == 0
    assert f("CCCC") == 1
    assert f("CGCGGCGCGGGGCCC") == 1

    assert f("AUCG") == 1 / 2
    assert f("UUUUU") == 0

    # Empty sequences are allowed, per the FASTG spec:
    # "Note that empty base_strings are allowed."
    assert f("") == 0


def test_invalid_char():
    for bad in (
        "acgta",
        "ACGtA",
        "uuuu",
        "asdofijodsif",
        " ",
        "acgtttttt\t",
        "ac\ngt",
        "123",
    ):
        with pytest.raises(ValueError) as ei:
            f(bad)
        assert str(ei.value) == (
            'Sequence "{}" contains character(s) not in the alphabet '
            "{{A, C, G, T, U}}."
        ).format(bad)
