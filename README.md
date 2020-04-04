# Transfer Integral

Compute the transfer integral between molecular orbitals (MOs) {i} on
monomer A and MOs {j} on monomer B.

If {i} and {j} are not specified, the default values are the HOMO (H) and
LUMO (L) on A and B. Valid MO names are e.g. "H", "L", "H-1" and "L+2", etc.
Multiple values can be specified using a comma e.g. "H-1,H,L,L+2"

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

    NOTE: Ca is an [n, n] matrix and Cb is an [m, m] matrix, where n and m
          are the no. of MO coefficients of monomers A and B, respectively.
          O and F, however, are [n+m, n+m] matrices. Therefore, only a subset
          of those matrices are used in computing Ei, Ej, Sij, Jij, etc.

          E.g., in Ei = <Ca_i|F|Ca_i>,
                the F is actually F[:n,:n]

                in Ej = <Cb_j|F|Cb_j>,
                the F is actually F[n:,n:]

                in Jij = <Ca_i|F|Cb_j>,
                the F is actually F[:n,n:]

Using Sij, Ei, Ej, Jij from above, the transfer integral (Jij_eff) is
calculated from:

    Jij_eff = (Jij - 0.5*(Ei+Ej)*Sij) / (1 - Sij*Sij)
