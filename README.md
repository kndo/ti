# Transfer Integral

Ti is a program for calculating the transfer integral between the molecular
orbitals (MOs) on two neighboring molecules (commonly known as a dimer complex).

If the MOs are not specified, it defaults to the HOMO (H) and LUMO (L) for both
molecules. Multiple MO values can be specified using a comma with no spaces,
e.g. 'H-2,H-1,H'.

The program takes in an YAML input file that specifies all the necessary log
and output files from a set of QM calculations, such as those from Gaussian09.

### Usage:
```
$ ti --help
Usage: ti [OPTIONS] INFILE

Options:
  -o, --outfile PATH  Output file
  -a TEXT             Molecular orbital(s) on monomer A
  -b TEXT             Molecular orbital(s) on monomer B
  --help              Show this message and exit.
```

```
$ ti in/A_B_0.yaml -a H-1,H,L -b H-1,H,L
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
  H-1    H-1     0.00082   -9.82493   -9.82493    -10.264     -9.82492     -9.82492         -2.223
  H-1    H      -0.00000   -9.82493   -7.43996      0.000     -9.82493     -7.43996          0.000
  H-1    L       0.00000   -9.82493   -0.22039     -0.000     -9.82493     -0.22039          0.000
  H      H-1    -0.00000   -7.43996   -9.82493      0.000     -7.43996     -9.82493          0.000
  H      H      -0.02422   -7.43996   -7.43996    236.467     -7.43860     -7.43860         56.306
  H      L      -0.00000   -7.43996   -0.22039      0.000     -7.43996     -0.22039         -0.000
  L      H-1    -0.00000   -0.22039   -9.82493      0.000     -0.22039     -9.82493         -0.000
  L      H      -0.00000   -0.22039   -7.43996      0.000     -0.22039     -7.43996          0.000
  L      L      -0.10621   -0.22039   -0.22039    164.697     -0.20521     -0.20521        142.901
```

Included in this repo is an example set of calculations for the ethylene
(H2C=CH2) dimer, where the monomers are rotated from 0 to 90 degress w.r.t.
each other. Monomer and dimer calculations are in `/mon` and `/dim`,
respectively. Input and output files are in `/in` and `/out`, respectively.
