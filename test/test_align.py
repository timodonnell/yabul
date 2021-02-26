# Test sequence alignment routines
# This test is meant to run using pytest, https://docs.pytest.org/en/stable/
#
import os

from yabul import align_pair


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_align_pair_protein_global():
    query_seq = "TESTEEEE"
    reference_seq = "ETESTQE"

    result = align_pair(query_seq, reference_seq)
    assert len(result.query) == len(result.reference)
    assert result.query.find("TEST") == result.reference.find("TEST")
    print("\n")
    print(result)


def test_align_pair_protein_local():
    query_seq = "TESTEEEE"
    reference_seq = "ETESTQE"

    result = align_pair(query_seq, reference_seq, local=True)
    assert len(result.query) == len(result.reference)
    assert result.query.find("TEST") == result.reference.find("TEST")
    assert result.query.startswith("TEST")
    assert result.reference.startswith("TEST")
    print("\n")
    print(result)


def test_align_pair_dna_local():
    query_seq = "CCAGATATACGGAA"
    reference_seq = "CGTATACGCAA"

    result = align_pair(
        query_seq, reference_seq, local=True, substitution_matrix="dnafull")
    print("\n")
    print(result)
    assert len(result.query) == len(result.reference)
    assert result.query.find("TATA") == result.reference.find("TATA")
    assert result.query.startswith("TATA")
    assert result.reference.startswith("TATA")


def test_align_pair_dna_global():
    query_seq = "CCAGATATACGGAA"
    reference_seq = "CGTATACGCAA"

    result = align_pair(
        query_seq, reference_seq, substitution_matrix="dnafull")
    print("\n")
    print(result)
    assert len(result.query) == len(result.reference)
    assert result.query.find("TATA") == result.reference.find("TATA")



