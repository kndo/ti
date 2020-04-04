#!/usr/bin/env python

import click
import yaml

from .inputs import validate
from .ti import compute_ti
from .ti import get_homo_number
from .ti import get_interaction_matrix
from .ti import get_mo_coeffs
from .ti import get_mo_number


@click.command()
@click.argument(
    'infile',
    type=click.Path(exists=True)
)
@click.option(
    '-o',
    '--outfile',
    type=click.Path(),
    help='Output file'
)
@click.option(
    '--outfile-format',
    type=click.Choice(['log', 'json', 'yaml'], case_sensitive=False),
    default='log',
    help='Output file format'
)
@click.option(
    '-a',
    'mo_i_a',
    default='H,L',
    help='Molecular orbital(s) on molecule A'
)
@click.option(
    '-b',
    'mo_j_b',
    default='H,L',
    help='Molecular orbital(s) on molecule B'
)
def main(infile, outfile, outfile_format, mo_i_a, mo_j_b):

    with open(infile) as f:
        inputs = yaml.load(f, Loader=yaml.FullLoader)
    validate(inputs)
    # validate_mo(mo_i_a)
    # validate_mo(mo_j_b)

    # ti = TransferIntegral(**inputs)
    # ti.compute(mo_i_a, mo_j_b)
    # ti.write(output)

    compute_ti(mo_i_a, mo_j_b, **inputs)


if __name__ == '__main__':
    main()
