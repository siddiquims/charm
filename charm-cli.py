#!/usr/bin/env python

"""charm-cli.py: Simple command line interface for CHarm."""

import argparse
import logging

from libcharm import LibCHarm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-o', '--output', type=str, help='path to output file')
    parser.add_argument('-f', '--frequency', action='store_true', help='use frequency/1000 instead of fraction')
    parser.add_argument('origin', type=int, help='species id of origin organism taken from '
                                                 '\'http://www.kazusa.or.jp/codon\' (e.g. \'83333\' for E. coli K12)')
    parser.add_argument('host', type=int, help='species id of host organism taken from '
                                               '\'http://www.kazusa.or.jp/codon\' (e.g. \'83333\' for E. coli K12)')
    parser.add_argument('input', type=str, help='input file in FASTA format')
    args = parser.parse_args()

    logger = logging.getLogger('charm-cli')
    logger.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    try:
        fh = logging.FileHandler('charm-cli.log', 'w')
        fh.setLevel(logging.INFO)
        logger.addHandler(fh)
    except IOError as e:
        logger.warning(
            'WARNING: Cannot create log file! Run charm-cli from a directory to which you have write access.')
        logger.warning(e)
        pass

    charm = LibCHarm()

    sequence = LibCHarm.Sequence(charm.open_input_file(args.input), args.origin, args.host, args.frequency)

    harmonized_codons = sequence.get_harmonized_codons()
    verify_sequence = sequence.verify_harmonized_sequence()

    logger.info('SUMMARY:\n')
    if verify_sequence:
        text = 'Success! Translation of harmonized and original sequence match:\n\n' \
               '{}\n'.format(sequence.harmonized_translated_sequence)
        logger.info(text)
    else:
        logger.error('ERROR: Translations of harmonized and original sequence DO NOT match!')
    logger.info('Harmonized codons: {}\n'.format(len(harmonized_codons)))

    table_header = '{:<10} {:^3} {:^4}    {:^4} {:^7} {:>6} {:<7} {:>6}'.format('position', 'aa', 'orig', 'new',
                                                                                'initial', 'final', 'origin', 'target')
    logger.info(table_header)

    for c in sequence.codons:
        if str(c['original']) != str(c['new']):
            line = '{:<10} {:^3} {:<4} -> {:<4} {:<5.2f} -> {:<3.2f}  {:<5.2f} -> {:<3.2f}'.format(c['position'],
                                                                                                   c['aa'],
                                                                                                   c['original'],
                                                                                                   c['new'],
                                                                                                   c['initial_df'],
                                                                                                   c['final_df'],
                                                                                                   c['origin_f'],
                                                                                                   c['target_f'])
        else:
            line = '{:<10} {:^3} {:<12} {:<5.2f}          {:<5.2f} -> {:<3.2f}'.format(c['position'], c['aa'],
                                                                                       c['original'], c['initial_df'],
                                                                                       c['origin_f'], c['target_f'])

        logger.info(line)

    logger.info('\nCodon-harmonized sequence:\n\n{}'.format(sequence.harmonized_sequence))
    exit(0)


if __name__ == "__main__":
    main()


