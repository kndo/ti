"""Compute the transfer integral between molecular orbitals (MOs) {i} on
monomer A and MOs {j} on monomer B.

If {i} and {j} are not specified, the default values are the HOMO (H) and
LUMO (L) on A and B. Valid MO names are e.g. "H", "L", "H-1" and "L+2", etc.
Multiple values are specified using a comma with no spaces, e.g. "H-1,H,L,L+2".

The required fields in the input file are:
    - program:        software that was used for the QM calculations
    - fock:           file containing fock matrix elements
    - overlap:        file containing overlap matrix elements
    - mo_a:           file containing MO coefficients of monomer A
    - mo_b:           file containing MO coefficients of monomer B
    - mo_type:        file type of the MO coefficients files
    - log_a:          log file from monomer A calculation
    - log_b:          log file from monomer B calculation

The overlap values are determined from: Sij = <Ca_i|O|Cb_j>, where O is the
overlap matrix and Ca_i and Cb_j are the matrices containing the MO coefficients
of monomers A and B, respectively.

Likewise, the energy values are determined using the Fock matrix:
    - Ei        = <Ca_i|F|Ca_i>
    - Ej        = <Cb_j|F|Cb_j>
    - Jij (Eij) = <Ca_i|F|Cb_j>

NB: Ca is an [n, n] matrix and Cb is an [m, m] matrix, where n and m are the
    no. of MO coefficients of monomers A and B, respectively. O and F, however,
    are [n+m, n+m] matrices. Therefore, only a subset of those matrices are used
    in computing Ei, Ej, Sij, Jij, etc.

    For example:
        in Ei = <Ca_i|F|Ca_i>,
        the F is actually F[:n,:n]

        in Ej = <Cb_j|F|Cb_j>,
        the F is actually F[n:,n:]

        in Jij = <Ca_i|F|Cb_j>,
        the F is actually F[:n,n:]

Using Sij, Ei, Ej, Jij from above, the transfer integral (Jij_eff) is
calculated from:

    Jij_eff = (Jij - 0.5*(Ei+Ej)*Sij) / (1 - Sij*Sij)
"""

import itertools
import os
import sys

import numpy as np


HARTREE_TO_EV = 27.2114

def get_homo_number(log_file, program):
    n = 0

    with open(log_file, 'r') as f:

        if program == 'g09':
            for line in f:
                if 'Alpha  occ. eigenvalues' in line:
                    cols = line.split()
                    n += len(cols[4:])
                elif 'Alpha virt. eigenvalues' in line:
                    break

    return n


def get_mo_number(mo_name, homo_number):
    n = homo_number

    if mo_name.startswith('H'):
        if '-' in mo_name:
            offset = int(mo_name.split('-')[-1])
            n -= offset
    elif mo_name.startswith('L'):
        n += 1  # lumo level
        if '+' in mo_name:
            offset = int(mo_name.split('+')[-1])
            n += offset

    return n


def get_mo_coeffs(mo_file, mo_file_type, program):
    n_mo = 0
    coeffs = []

    with open(mo_file, 'r') as f:

        if program == 'g09':

            if mo_file_type == 'punch':
                n_coeffs_per_line = 5
                coeff_len = 15  # characters; use b/c no separation between coeffs

                for line in f.readlines()[1:]:  # Ignore first line
                    if 'Alpha MO' in line:
                        n_mo += 1
                    else:
                        for coeff in grouper(line, coeff_len):
                            if coeff != '\n':
                                coeff = coeff.replace('D', 'E')
                                coeff = float(coeff)
                                coeffs.append(coeff)

            elif mo_file_type == 'fchk':
                is_coeff_line = False
                for line in f:
                    if is_coeff_line:
                        cols = line.split()
                        for coeff in cols:
                            coeffs.append(coeff)
                        if len(coeffs) == n_coeffs:
                            break

                    elif 'Number of basis functions' in line:
                        cols = line.split()
                        n_mo = int(cols[-1])
                        n_coeffs = n_mo * n_mo

                    elif 'Alpha MO coefficients' in line:
                        is_coeff_line = True  # The following lines contain coeffs
                        continue

    C = np.matrix(coeffs)
    C = C.reshape(n_mo, n_mo)

    return C


def get_interaction_matrix(matrix_file, n_mo_a, n_mo_b, program):
    ndim = n_mo_a + n_mo_b
    zeros = np.zeros((ndim, ndim))
    M = np.matrix(zeros)

    with open(matrix_file, 'r') as f:

        if program == 'g09':
            numel = []
            for line in f.readlines()[3:]:  # Ignore first lines
                cols = line.split()
                for coeff in cols:
                    numel.append(float(coeff))

            k = -1
            # Form matrix by filling bottom diagonal elements first
            for i in range(ndim):
                for j in range(i+1):
                    k += 1
                    M[i,j] = numel[k]

            # Form top diagonal elements by symmetry
            for i in range(ndim):
                for j in range(i+1, ndim):
                    M[i,j] = M[j,i]

    return M


# def mo_sign(C, Cr):
#     CrT = np.transpose(Cr)
#     sign = C*CrT

#     m, n = sign.shape
#     for i in range(m):
#         for j in range(n):
#             if sign[i,j] < 0:
#                 sign[i,j] = -1.0
#             elif sign[i,j] > 0:
#                 sign[i,j] = 1.0

#     return sign


output_template = """INPUT FILES:
  Fock matrix = {fock}
  Overlap matrix = {overlap}
  Mon. A MOs = {mo_a}
  Mon. A log = {log_a}
  Mon. B MOs = {mo_b}
  Mon. B log = {log_b}

MONOMER A:
  No. of MOs = {n_mo_a}
  HOMO = {n_homo_a} ; LUMO = {n_lumo_a}

MONOMER B:
  No. of MOs = {n_mo_b}
  HOMO = {n_homo_b} ; LUMO = {n_lumo_b}

TRANSFER INTEGRAL:
  MO(A)  MO(B)
  i      j           Sij    Ei [eV]    Ej [eV]  Jij [meV]  Ei_eff [eV]  Ej_eff [eV]  Jij_eff [meV]"""


def compute_ti(mo_i_a, mo_j_b, outfile,
               program, fock, overlap, mo_a, mo_b, mo_type, log_a, log_b):

    if os.path.isfile(outfile):
        print('Outfile %s already exists. Delete [y/n]?' % outfile)
        ans = sys.stdin.readline()
        if not ans.startswith('y'):
            print('Exiting ...')
            sys.exit(1)
        print('')

    fout = open(outfile, 'w')

    def print2(line):
        print(line)
        print(line, file=fout)

    n_homo_a = get_homo_number(log_a, program)
    n_homo_b = get_homo_number(log_b, program)

    mo_i_a = mo_i_a.split(',')
    mo_j_b = mo_j_b.split(',')

    mo_name_a = {}
    mo_num_a = []
    for name in mo_i_a:
        num = get_mo_number(name, n_homo_a)
        mo_name_a[num] = name
        mo_num_a.append(num)

    mo_name_b = {}
    mo_num_b = []
    for name in mo_j_b:
        num = get_mo_number(name, n_homo_b)
        mo_name_b[num] = name
        mo_num_b.append(num)

    Ca = get_mo_coeffs(mo_a, mo_type, program)
    Cb = get_mo_coeffs(mo_b, mo_type, program)

    n = len(Ca)
    m = len(Cb)

    F = get_interaction_matrix(fock, n, m, program)
    O = get_interaction_matrix(overlap, n, m, program)

    CaT = np.transpose(Ca)
    CbT = np.transpose(Cb)

    # <Ca|O|Cb>
    Sab = Ca*O[:n,n:]*CbT

    # <Ca|F|Ca>
    Ea = Ca*F[:n,:n]*CaT
    Ea *= HARTREE_TO_EV

    # <Cb|F|Cb>
    Eb = Cb*F[n:,n:]*CbT
    Eb *= HARTREE_TO_EV

    # <Ca|F|Cb>
    Eab = Ca*F[:n,n:]*CbT
    Eab *= HARTREE_TO_EV

    output = output_template.format(
        fock=fock,
        overlap=overlap,
        mo_a=mo_a,
        mo_b=mo_b,
        log_a=log_a,
        log_b=log_b,
        n_mo_a=n,
        n_mo_b=m,
        n_homo_a=n_homo_a,
        n_lumo_a=n_homo_a+1,
        n_homo_b=n_homo_b,
        n_lumo_b=n_homo_b+1,
    )
    print2(output)

    # Use dummy vars p and q, b/c matrices start at index 0
    for p in mo_num_a:
        i = p - 1

        for q in mo_num_b:
            j = q - 1

            Sij = Sab[i,j]
            Ei = Ea[i,i]
            Ej = Eb[j,j]
            Jij = Eab[i,j]

            Ei_eff = 0.5 * ((Ei+Ej) - 2*Jij*Sij + (Ei-Ej)*np.sqrt(1-Sij*Sij)) / (1-Sij*Sij)
            Ej_eff = 0.5 * ((Ei+Ej) - 2*Jij*Sij - (Ei-Ej)*np.sqrt(1-Sij*Sij)) / (1-Sij*Sij)

            Jij_eff = (Jij - 0.5*(Ei+Ej)*Sij) / (1-Sij*Sij)

            print2('  %-6s %-6s %8.5f %10.5f %10.5f %10.3f %12.5f %12.5f %14.3f' % \
                   (mo_name_a[p], mo_name_b[q], Sij, Ei, Ej, Jij*1e3, Ei_eff, Ej_eff, Jij_eff*1e3))

    print2('')


def grouper(line, nchars):
    it = iter(line)
    while True:
       chunk = tuple(itertools.islice(it, nchars))
       if not chunk:
           return
       yield ''.join(chunk)
