# Test fasta reading and writing
# This test is meant to run using pytest, https://docs.pytest.org/en/stable/
#
import shutil
import os
import tempfile

from numpy.testing import assert_array_equal

from yabul import read_fasta, write_fasta


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_read_protein_basic():
    df = read_fasta(os.path.join(DATA_DIR, "protein.fasta"))
    assert len(df) == 4
    assert df.loc["TEST_PROTEIN", "description"] == "TEST_PROTEIN Just a test"
    assert df.loc["TEST_PROTEIN", "sequence"] == "TEST"


def test_read_and_write():
    for test_file in ["protein.fasta", "transcripts.fasta"]:
        df = read_fasta(os.path.join(DATA_DIR, test_file))

        temp_dir = tempfile.mkdtemp()
        try:
            # Uncompressed write
            out = os.path.join(temp_dir, "out.fasta")
            write_fasta(out, df.set_index("description").sequence.items())
            df2 = read_fasta(out)
            assert_array_equal(df, df2)

            # Compressed write
            out_gz = os.path.join(temp_dir, "out.fasta.gz")
            write_fasta(out_gz, df.set_index("description").sequence.items())
            # Check it really did compress:
            assert os.path.getsize(out_gz) < os.path.getsize(out)
            df2 = read_fasta(out_gz)
            assert_array_equal(df, df2)
        finally:
            shutil.rmtree(temp_dir)



