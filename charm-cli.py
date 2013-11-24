#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""charm-cli.py: Simple command line interface for CHarm."""

import argparse
import logging

try:
    import matplotlib

    matplotlib.use('Agg')
    matplotlib.rc('font', **{'sans-serif': 'DejaVu Sans',
                             'serif': 'DejaVu Serif',
                             'family': 'sans-serif'})
    import matplotlib.pyplot
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)
try:
    import numpy
except ImportError as e:
    print('ERROR: {}'.format(e.msg))
    exit(1)

from libcharm import LibCHarm


def autolabel(rects, ax, labels, vertical=True):
    if vertical:
        rotation = 'vertical'

    if len(labels) == len(rects):
        heights = []
        for rect in rects:
            height = rect.get_height()
            heights.append(height)

        max_height = max(heights)

        for rect in rects:
            i = rects.index(rect)
            label = labels[i]
            height = rect.get_height()
            if height > 0:
                y = 1.05 * height
            else:
                y = 0.02 * max_height
            ax.text(rect.get_x() + rect.get_width() / 2., y, str(label),
                    ha='center', va='bottom', rotation=rotation, size='x-small')


def plot_codon_usage(sequence, ax):#, prefix=None):

    x1 = x2 = numpy.arange(len(sequence.codons))
    bar_width = 0.5
    xlabels = []

    origin_f = []
    target_f = []

    for c in sequence.codons:
        origin_f.append(c['origin_f'])
        target_f.append(c['target_f'])
        xlabels.append(c['aa'])

    origin_f = numpy.array(origin_f)
    target_f = numpy.array(target_f)

    p1 = ax.bar(x1, origin_f, color='b', width=bar_width)
    p2 = ax.bar(x2 + (0.5 * bar_width), target_f, color='r', width=bar_width)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='both', direction='out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.set_xticks(x1 + bar_width / 2)
    ax.set_xticklabels(xlabels)
    ax.set_xlabel('amino acid')

    if sequence.use_frequency:
        ax.set_ylabel('codon usage [frequency/1000]')
    else:
        ax.set_ylabel('codon usage [fraction]')
    ax.legend((p1, p2), ('Origin organism', 'Host organism'), loc=2, bbox_to_anchor=(1, 1))
    ax.hlines(sequence.lower_threshold, 0, len(x1), colors='k', linestyles='dotted', **{'linewidth': 2})

    if not sequence.use_frequency:
        major_locator = matplotlib.ticker.MultipleLocator(0.1)
        minor_locator = matplotlib.ticker.MultipleLocator(0.01)
    else:
        major_locator = matplotlib.ticker.MultipleLocator(10)
        minor_locator = matplotlib.ticker.MultipleLocator(1)

    ax.yaxis.set_major_locator(major_locator)
    ax.yaxis.set_minor_locator(minor_locator)


def plot_codon_usage_differences(sequence, ax):#, prefix=None):
    x1 = numpy.arange(len(sequence.codons))

    bar_width = 0.8
    xlabels = []

    df = []
    bar_labels = []

    for c in sequence.codons:
        df.append(c['final_df'])
        xlabels.append(c['aa'])
        label = u'{} → {}'.format(c['original'], c['new'])
        bar_labels.append(label)

    df = numpy.array(df)

    p1 = ax.bar(x1, df, color='b', width=bar_width)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='both', which='both', direction='out')
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    ax.set_xticks(x1 + bar_width / 2)
    ax.set_xticklabels(xlabels, **{'family': 'monospace'})
    ax.set_xlabel('amino acid')

    ax.set_ylabel(r'Differential codon usage $f_{origin} - f_{host}$')

    if not sequence.use_frequency:
        major_locator = matplotlib.ticker.MultipleLocator(0.05)
        minor_locator = matplotlib.ticker.MultipleLocator(0.01)
    else:
        major_locator = matplotlib.ticker.MultipleLocator(10)
        minor_locator = matplotlib.ticker.MultipleLocator(1)

    ax.yaxis.set_major_locator(major_locator)
    ax.yaxis.set_minor_locator(minor_locator)

    autolabel(p1, ax, bar_labels, vertical=True)


def plot(sequence, prefix=None):
    if prefix:
        filename = '{}_charm_results.svg'.format(prefix)
    else:
        filename = 'charm_results.svg'

    fig, axarr = matplotlib.pyplot.subplots(2, figsize=(50, 20), dpi=300)

    plot_codon_usage(sequence, axarr[0])
    plot_codon_usage_differences(sequence, axarr[1])

    matplotlib.pyplot.savefig(filename, format='svg', orientation='landscape', papertype='a4')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true', help='increase output verbosity')
    parser.add_argument('-p', '--prefix', type=str, help='prefix for output files')
    parser.add_argument('-f', '--frequency', action='store_true', help='use frequency/1000 instead of fraction')
    parser.add_argument('-t', '--threshold', type=float,
                        help='Lower threshold of codon usage. Defaults to 0.1 and 5 for fraction and frequency respectively')
    parser.add_argument('-to', '--translation_table_origin', type=int,
                        help='id of translation table; Default is: standard genetic code = 1; '
                             'id corresponds to \'trans_table\' '
                             'on http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi')
    parser.add_argument('-th', '--translation_table_host', type=int,
                        help='id of translation table; Default is: standard genetic code = 1; '
                             'id corresponds to \'trans_table\' '
                             'on http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi')
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
        if args.prefix:
            log_filename = '{}_charm-cli.log'.format(args.prefix)
        else:
            log_filename = 'charm-cli.log'
        fh = logging.FileHandler(log_filename, 'w')
        fh.setLevel(logging.INFO)
        logger.addHandler(fh)
    except IOError as e:
        logger.warning('WARNING: Cannot create log file! Run charm-cli from a directory to '
                       'which you have write access.')
        logger.warning(e.msg)
        pass

    charm = LibCHarm()

    if args.translation_table_origin:
        translation_table_origin = args.translation_table_origin
    else:
        translation_table_origin = 1

    if args.translation_table_host:
        translation_table_host = args.translation_table_host
    else:
        translation_table_host = 1

    if args.threshold:
        lower_threshold = args.threshold
    elif args.frequency:
        lower_threshold = 5
    else:
        lower_threshold = 0.1

    sequence = charm.Sequence(charm.open_input_file(args.input), args.origin, args.host,
                              translation_table_origin=translation_table_origin,
                              translation_table_host=translation_table_host,
                              use_frequency=args.frequency,
                              lower_threshold=lower_threshold)

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

    #    logger.info('Plotting data. Resulting files will be prefixed with \'{}\''.format(args.prefix))

    plot(sequence, args.prefix)

    exit(0)


if __name__ == "__main__":
    main()


