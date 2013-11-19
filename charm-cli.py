#!/usr/bin/env python

"""charm-cli.py: Simple command line interface for CHarm."""

import argparse

import libcharm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-o', '--output', type=str, help='path to output file')
    parser.add_argument('origin', type=int, help='species id of origin organism taken from '
                                                 '\'http://www.kazusa.or.jp/codon\' (e.g. \'83333\' for E. coli K12)')
    parser.add_argument('host', type=int, help='species id of host organism taken from '
                                               '\'http://www.kazusa.or.jp/codon\' (e.g. \'83333\' for E. coli K12)')
    parser.add_argument('input', type=str, help='input file in FASTA format')
    args = parser.parse_args()

    sequence = libcharm.Sequence(libcharm.open_input_file(args.input), args.origin, args.host)

    harmonized_codons = sequence.get_harmonized_codons()
    verify_sequence = sequence.verify_harmonized_sequence()

    print('SUMMARY:\n')
    if verify_sequence:
        print('Success! Translation of harmonized and original sequence match!')
    else:
        print('ERROR: Translations of harmonized and original sequence DO NOT match!')
    print('Harmonized codons: {}\n'.format(len(harmonized_codons)))
    print('Codon-harmonized sequence:\n\n{}'.format(sequence.harmonized_sequence))
    exit(0)


if __name__ == "__main__":
    main()


