from LibCharm import IO

def test_load_file():
    print(IO.load_file('tests/test_sequence.fasta', file_format="fasta"))
