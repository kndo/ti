# Transfer Integral

`ti` is a program that calculates the transfer integral between the molecular
orbitals (MOs) {i} on monomer A and {j} on monomer B.

It takes in an YAML input file that specifies all the necessary log and output
files from a QM calculation, such as those from Gaussian09.

Usage:
```
ti in/A_B_0.yaml -o out/A_B_0.out
```

Output:
```
INPUT FILES:
  Fock matrix = dim/A_B_0.fm
  Overlap matrix = dim/A_B_0.om
  Mon. A MOs = mon/A.mo
  Mon. A log = mon/A.log
  Mon. B MOs = mon/B_0.mo
  Mon. B log = mon/B_0.log

MONOMER A:
  No. of MOs = 46
  HOMO = 8 ; LUMO = 9

MONOMER B:
  No. of MOs = 46
  HOMO = 8 ; LUMO = 9

TRANSFER INTEGRAL:
  MO(A)  MO(B)
  i      j           Sij    Ei [eV]    Ej [eV]  Jij [meV]  Ei_eff [eV]  Ej_eff [eV]  Jij_eff [meV]
  H      H      -0.02422   -7.43996   -7.43996    236.467     -7.43860     -7.43860         56.306
  H      L      -0.00000   -7.43996   -0.22039      0.000     -7.43996     -0.22039         -0.000
  L      H      -0.00000   -0.22039   -7.43996      0.000     -0.22039     -7.43996          0.000
  L      L      -0.10621   -0.22039   -0.22039    164.697     -0.20521     -0.20521        142.901
```

Included in this repo is an example set of calculations for the ethylene (H2C=CH2)
dimer, where the monomers are rotated from 0 to 90 degress w.r.t. each other.
Monomer and dimer calculations are in `/mon` and `/dim`, respectively.
Input and output files for `ti` are in `/in` and `/out`, respectively.
