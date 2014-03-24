from LibCharm.Sequence import Sequence
from LibCharm import IO

seq = IO.load_file('tests/test_sequence.fasta', file_format="fasta")


def test_sequence_1():
    assert Sequence(seq, 83333, 4227)


def test_sequence_2():
    assert Sequence(seq, 83333, 4227, strong_stop=False)


def test_sequence_3():
    assert Sequence(seq, 83333, 4227, use_frequency=True)


def test_sequence_4():
    assert Sequence(seq, 83333, 4227, use_replacement_table=False)


def test_sequence_5():
    assert Sequence(seq, 83333, 4227, lower_alternative=False)


def test_sequence_6():
    assert Sequence(seq, 83333, 4227, strong_stop=False, use_frequency=True)


def test_sequence_7():
    assert Sequence(seq, 83333, 4227, strong_stop=False, use_frequency=True, use_replacement_table=False)


def test_sequence_8():
    assert Sequence(seq, 83333, 4227, strong_stop=False, use_frequency=True, use_replacement_table=False,
                    lower_alternative=False)