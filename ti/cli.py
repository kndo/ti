#!/usr/bin/env python

import os

import click
import yaml

from .validate import valid_inputs
from .validate import valid_mo
from .ti import compute_ti


@click.command()
@click.argument(
    'infile',
    type=click.Path(exists=True),
)
@click.option(
    '-o',
    '--outfile',
    type=click.Path(),
    help='Output file',
)
@click.option(
    '-a',
    'mo_i_a',
    default='H,L',
    help='Molecular orbital(s) on monomer A',
)
@click.option(
    '-b',
    'mo_j_b',
    default='H,L',
    help='Molecular orbital(s) on monomer B',
)
def main(infile, outfile, mo_i_a, mo_j_b):

    with open(infile) as f:
        inputs = yaml.load(f, Loader=yaml.FullLoader)

    # TODO: Make all error message uniform
    valid_inputs(inputs)
    valid_mo(mo_i_a)
    valid_mo(mo_j_b)

    if not outfile:
        basename = os.path.basename(infile)
        cols = basename.split('.')[:-1]
        rootname = '.'.join(cols)
        outfile = rootname + '.out'

    compute_ti(mo_i_a, mo_j_b, outfile, **inputs)


if __name__ == '__main__':
    main()
