from LibCharm.CodonUsageTable import CodonUsageTable


def test_codonusagetable_fraction():
    assert CodonUsageTable('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333&aa=1&style=N')


def test_codonusagetable_frequency():
    assert CodonUsageTable('http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=83333&aa=1&style=N',
                           use_frequency=True)