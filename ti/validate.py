import sys

from voluptuous import Any
from voluptuous import Invalid
from voluptuous import IsFile
from voluptuous import Required
from voluptuous import REMOVE_EXTRA
from voluptuous import Schema


def valid_inputs(inputs):
    schema = Schema({
        Required('program'): Any('g09'),
        Required('fock'): IsFile(),
        Required('overlap'): IsFile(),
        Required('mo_a'): IsFile(),
        Required('mo_b'): IsFile(),
        Required('mo_type'): Any('fchk', 'punch'),
        Required('log_a'): IsFile(),
        Required('log_b'): IsFile(),
    }, extra=REMOVE_EXTRA)

    try:
        schema(inputs)
    except Invalid as err:
        err_path = [str(i) for i in err.path]
        print(f'Invalid field in input file: {err.msg} ({".".join(err_path)}).')
        sys.exit(1)


def valid_mo(mo):
    mos = mo.split(',')

    for mo in mos:
        # TODO: Use regex to check that MO is 'H-n' or 'L+m'
        pass
