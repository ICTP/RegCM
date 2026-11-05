#!/usr/bin/env python3
"""
Assemble the CORDEX chunks still in segments at the end of an experiment.

pycordexer writes each month as a segment and assembles a chunk (1 year
sub-daily, 5 years daily, 10 years monthly) as soon as all its records are
there. The chunks at the edges of the experiment are shorter and are never
complete: run this script once the whole experiment has been processed to
assemble them, and any other chunk left in segments. A chunk with records
missing inside it is not assembled and is reported.

Usage: cordex_finalize.py [--first-year Y] [--last-year Y] [--dry-run] dir
"""

import argparse
import logging
import sys

from utilities import chunks


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('output_dir',
                        help='The CORDEX output directory of pycordexer')
    parser.add_argument('--first-year', type=int, default=None,
                        help='First year of the experiment')
    parser.add_argument('--last-year', type=int, default=None,
                        help='Last year of the experiment')
    parser.add_argument('--dry-run', action='store_true',
                        help='List the chunks, write nothing')
    parser.add_argument('-v', '--verbosity',
                        choices=['debug', 'info', 'warning'], default='info')
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.verbosity.upper()),
                        format='%(levelname)s: %(message)s')
    written = chunks.finalize(args.output_dir, args.first_year,
                              args.last_year, args.dry_run)
    for item in written:
        print(item if not args.dry_run else
              'would assemble {1} {2}-{3} in {0}'.format(*item))
    left = chunks.remaining_segments(args.output_dir)
    if left and not args.dry_run:
        print('{} segments not assembled (see the errors above):'
              .format(len(left)), file=sys.stderr)
        for path in left:
            print('  ' + path, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
