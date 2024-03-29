# Test fasta reading and writing
# This test is meant to run using pytest, https://docs.pytest.org/en/stable/
#
import shutil
import os
import tempfile
import time

from numpy.testing import assert_array_equal

from yabul import read_fasta, write_fasta


DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
print(DATA_DIR)


def test_read_protein_basic():
    df = read_fasta(os.path.join(DATA_DIR, "protein.fasta"))
    assert len(df) == 4
    assert df.loc["TEST_PROTEIN", "description"] == "TEST_PROTEIN Just a test"
    assert df.loc["TEST_PROTEIN", "sequence"] == "TEST"

def test_read_protein_uniprot():
    start = time.time()
    df = read_fasta(os.path.join(DATA_DIR, "human.uniprot.one_per_gene.fasta.gz"))
    duration = time.time() - start
    print("Read time:", duration)
    assert len(df) == 20610
    
    df2 = read_fasta(
        os.path.join(DATA_DIR, "human.uniprot.one_per_gene.fasta.gz"),
        progress=False)
    assert_array_equal(df, df2)


    df3 = read_fasta(
        os.path.join(DATA_DIR, "human.uniprot.one_per_gene.fasta.gz"),
        chunksize=500,
        progress=False)
    assert_array_equal(df, df3)

    for num in [5, 499, 500, 501, 502]:
        df4 = read_fasta(
            os.path.join(DATA_DIR, "human.uniprot.one_per_gene.fasta.gz"),
            max_sequences=num,
            chunksize=500,
            progress=False)
        assert len(df4) == num
        assert_array_equal(df.head(num), df4)


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



EXPECTED_TRANSCRIPT = "NM_001126115.2"
EXPECTED_TRANSCRIPT_DESCRIPTION = "NM_001126115.2 Homo sapiens tumor protein p53 (TP53), transcript variant 5, mRNA"
EXPECTED_TRANSCRIPT_SEQUENCE = """
TCCTACAGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTG
GGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATG
ACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGC
ATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGT
GGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAAC
AGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTA
ATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGA
AGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCC
AACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCC
GTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGG
GAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGC
CATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGACATTCTCCACTTCTTGTTCCCCACT
GACAGCCTCCCACCCCCATCTCTCCCTCCCCTGCCATTTTGGGTTTTGGGTCTTTGAACCCTTGCTTGCA
ATAGGTGTGCGTCAGAAGCACCCAGGACTTCCATTTGCTTTGTCCCGGGGCTCCACTGAACAAGTTGGCC
TGCACTGGTGTTTTGTTGTGGGGAGGAGGATGGGGAGTAGGACATACCAGCTTAGATTTTAAGGTTTTTA
CTGTGAGGGATGTTTGGGAGATGTAAGAAATGTTCTTGCAGTTAAGGGTTAGTTTACAATCAGCCACATT
CTAGGTAGGGGCCCACTTCACCGTACTAACCAGGGAAGCTGTCCCTCACTGTTGAATTTTCTCTAACTTC
AAGGCCCATATCTGTGAAATGCTGGCATTTGCACCTACCTCACAGAGTGCATTGTGAGGGTTAATGAAAT
AATGTACATCTGGCCTTGAAACCACCTTTTATTACATGGGGTCTAGAACTTGACCCCCTTGAGGGTGCTT
GTTCCCTCTCCCTGTTGGTCGGTGGGTTGGTAGTTTCTACAGTTGGGCAGCTGGTTAGGTAGAGGGAGTT
GTCAAGTCTCTGCTGGCCCAGCCAAACCCTGTCTGACAACCTCTTGGTGAACCTTAGTACCTAAAAGGAA
ATCTCACCCCATCCCACACCCTGGAGGATTTCATCTCTTGTATATGATGATCTGGATCCACCAAGACTTG
TTTTATGCTCAGGGTCAATTTCTTTTTTCTTTTTTTTTTTTTTTTTTCTTTTTCTTTGAGACTGGGTCTC
GCTTTGTTGCCCAGGCTGGAGTGGAGTGGCGTGATCTTGGCTTACTGCAGCCTTTGCCTCCCCGGCTCGA
GCAGTCCTGCCTCAGCCTCCGGAGTAGCTGGGACCACAGGTTCATGCCACCATGGCCAGCCAACTTTTGC
ATGTTTTGTAGAGATGGGGTCTCACAGTGTTGCCCAGGCTGGTCTCAAACTCCTGGGCTCAGGCGATCCA
CCTGTCTCAGCCTCCCAGAGTGCTGGGATTACAATTGTGAGCCACCACGTCCAGCTGGAAGGGTCAACAT
CTTTTACATTCTGCAAGCACATCTGCATTTTCACCCCACCCTTCCCCTCCTTCTCCCTTTTTATATCCCA
TTTTTATATCGATCTCTTATTTTACAATAAAACTTTGCTGCCA
""".replace("\n", "").strip()


def test_read_transcripts():
    df = read_fasta(os.path.join(DATA_DIR, "transcripts.fasta"))
    assert df.loc[EXPECTED_TRANSCRIPT, "description"] == EXPECTED_TRANSCRIPT_DESCRIPTION
    assert df.loc[EXPECTED_TRANSCRIPT, "sequence"] == EXPECTED_TRANSCRIPT_SEQUENCE
