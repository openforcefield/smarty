# smirky sampling of Torsions

This is an example of how to use smirky, a command line tool for sampling chemical perception of Bonds, Angles, proper or improper Torsions, or van der Waal parameters. Default smirky behaivor only requires two inputs, this is an example of all input options into smirky. 

### Input files explained

* `atom_OR_bases.smarts` - element numbers that form the base of atoms, such as `"#6"` and their associated odds
* `atom_OR_decorators.smarts` - decorators and associated odds that are combined with element numbers such as `X4` in `"[#6X3,#7]"`
* `atom_AND_decorators.smarts` - decorators and associated odds for patterns that are "AND'd" to the end of an atom for example `r5` in `"[#6X4,#7X3;r5]"`
* `bond_OR_bases.smarts` - bond bases and their associated odds, that is '-', '=', ':', or '#' typically
* `bond_AND_decorators.smarts` - bond decorators that can be "AND'd" in a bond, such as '@' in `"[#6r6]-,:;@[#7r6]"`
* `atom_odds_forTorsions.smarts` - keywords or indices for atoms in torsions and odds of making changes to them
* `bond_odds_forTorsions.smarts` - keywords or indices for bonds in torsions and odds of making changes to them
* `initial_Torsions.smarts` - SMIRKS patterns for initial patterns
* `substitutions.smarts` - SMIRKS patterns and the short hand they can be replaced with

### Command line call

```
smirky --molecules AlkEthOH_test_filt1_ff.mol2 \
    --typetag Torsion \
    --atomORbases atom_OR_bases.smarts \
    --atomORdecors atom_OR_decorators.smarts \
    --atomANDdecors atom_AND_decorators.smarts \
    --bondORbase bond_OR_bases.smarts \
    --bondANDdecors bond_AND_decorators.smarts \
    --atomOddsFile atom_odds_forTorsions.smarts \
    --bondOddsFile bond_odds_forTorsions.smarts \
    --initialtypes initial_Torsions.smarts \
    --substitutions substitutions.smarts \
    --smirff forcefield/Frosst_AlkEthOH.ffxml \
    --iteratorsion 1000 \
    --temperature 0.001 \ 
    --verbose True \
    --output output
```

### Output files created
* output.log - detailed log of each iteration, changes made and if it was accepted or rejected
* output.csv - a "trajectory" file that describes the torsions at each iteration
* output.pdf - plot showing the overall score vs iteration
* output_results.smarts - smarts file showing the file SMIRKS and their matched results

### Detailed output explained

Here is a segment of output.log with explaination of what happens in a smirky simulation

##### Match initial input

Type initial parameters 
```
INDEX      TORSIONS  MOLECULES   TYPE NAME: SMIRKS
    1 :          0          0 | 0: [*:1]~[*:2]~[*:3]~[*:4]
    2 :       1737         42 | C-C: [*:1]~[#6:2]~[#6:3]~[*:4]
    3 :        438         42 | C-O: [*:1]~[#6:2]~[#8:3]~[*:4]
TOTAL :       2175         42
```
Remove elements that are not used in this molecule set (remember AlkEthOH only has carbon, oxygen, and hydrogen)
```
removing unused element ([#5]) from list
removing unused element ([#7]) from list
removing unused element ([#9]) from list
removing unused element ([#14]) from list
removing unused element ([#15]) from list
removing unused element ([#16]) from list
removing unused element ([#17]) from list
removing unused element ([#35]) from list
removing unused element ([#53]) from list
```
##### Comparing to SMIRNOFF99Frosst

Use the forcefield tools to type all molecules with SMIRNOFF reference.
Compare reference types to initial parameter types

```
Creating labeler from forcefield/Frosst_AlkEthOH.ffxml...
Creating graph matching current types with reference types...
Graph creation took 0.304 s
Computing maximum weight match...
Maximum weight match took 0.001 s
PROPOSED:
Torsion type matches:
0: [*:1]~[*:2]~[*:3]~[*:4]                                       no match
C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                                   matches                           t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]:      574 Torsion    types matched
C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                                   matches                         t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]:      156 Torsion    types matched
730 / 2175 total Torsions match (33.563 %)
```
Show current statistics before sampling begins
```
INDEX      TORSIONS  MOLECULES   TYPE NAME: SMIRKS                                  REF TYPE: SMIRKS                                   FRACTION OF REF TYPED MOLECULES MATCHED
    1 :          0          0 | 0: [*:1]~[*:2]~[*:3]~[*:4]
    2 :       1737         42 | C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                     t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]                 574 /     574 (100.000%)
    3 :        438         42 | C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                     t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]               156 /     156 (100.000%)
TOTAL :       2175         42 |                                                        730 /     2175 match (33.563 %)
```

##### Example move to generate a new Torsion

Create a new torsion, in this case by changing the 4th atom from generic (*) to an oxygen not bound to hydrogen (`#8H0`)

```
Iteration 1 / 1000
Attempting to create new subtype: '4778' ([*:1]~[#6:2]~[#6:3]~[#8!H0:4]) from parent type 'C-C' ([*:1]~[#6:2]~[#6:3]~[*:4])
    Probability of making this environment is 0.004 %Proposal is valid...
```
Compare proposed types to the SMIRNOFF reference types
```
Creating graph matching current types with reference types...
Graph creation took 0.176 s
Computing maximum weight match...
Maximum weight match took 0.001 s
PROPOSED:
Torsion type matches:
C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                                   matches                           t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]:      574 Torsion    types matched
C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                                   matches                         t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]:      156 Torsion    types matched
4778: [*:1]~[#6:2]~[#6:3]~[#8!H0:4]                              matches                         t0012: [#1:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]:      190 Torsion    types matched
920 / 2175 total Torsions match (42.299 %)
```
##### Using temperature and score to accept or reject move
Use change in score and temperature to calculate the probability of accepting the move.
A move with an increased score will always be accepted, the higher the temperature the
more probable a move with a decreased score will be accepted
```
Proposal score: 730 >> 920 : log_P_accept = 8.73563e+01
Accepted.
INDEX      TORSIONS  MOLECULES   TYPE NAME: SMIRKS                                  REF TYPE: SMIRKS                                   FRACTION OF REF TYPED MOLECULES MATCHED
    1 :       1436         42 | C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                     t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]                 574 /     574 (100.000%)
    2 :        438         42 | C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                     t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]               156 /     156 (100.000%)
    3 :        301         42 | 4778: [*:1]~[#6:2]~[#6:3]~[#8!H0:4]                t0012: [#1:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]               190 /     307 ( 61.889%)
TOTAL :       2175         42 |                                                        920 /     2175 match (42.299 %)

```
Hierarchy shows which parent types lead to the generation of child types
```
Torsion type hierarchy:
    C-C ([*:1]~[#6:2]~[#6:3]~[*:4])
        4778 ([*:1]~[#6:2]~[#6:3]~[#8!H0:4])
    C-O ([*:1]~[#6:2]~[#8:3]~[*:4])
```
##### Final Iteration in this example
```
Iteration 999 / 1000
Attempting to destroy type 1876 : [#1:1]~[#6:2]~[#6:3]~[#1:4]...
Proposal is valid...
Creating graph matching current types with reference types...
Graph creation took 0.249 s
Computing maximum weight match...
Maximum weight match took 0.004 s
PROPOSED:
Torsion type matches:
C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                                   matches                           t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]:      574 Torsion    types matched
C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                                   matches                         t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]:      156 Torsion    types matched
4808: [*:1]~[#6:2]~[#8:3]~[#1!X4:4]                              matches                          t0002: [a,A:1]-[#6X4:2]-[#8X2:3]-[#1:4]:      101 Torsion    types matched
8090: [#6!H1:1]~[#6:2]~[#8:3]~[#1!X4:4]                          matches                         t0006: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#1:4]:       87 Torsion    types matched
7751: [*:1]~[#6:2]~[#6:3]~[#6:4]                                 matches                         t0001: [a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]:      146 Torsion    types matched
1068: [#6!H3:1]~[#6:2]~[#6:3]~[#6:4]                             matches                       t0007: [#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]:      131 Torsion    types matched
6774: [#1H0:1]~[#6:2]~[#6:3]~[#6:4]                              matches                         t0005: [#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]:      552 Torsion    types matched
8025: [#6:1]~[#6:2]~[#8:3]~[#6!H3:4]                             matches                       t0008: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]:       66 Torsion    types matched
1813 / 2175 total Torsions match (83.356 %)
Proposal score: 2120 >> 1813 : log_P_accept = -1.41149e+02
Rejected.
INDEX      TORSIONS  MOLECULES   TYPE NAME: SMIRKS                                  REF TYPE: SMIRKS                                   FRACTION OF REF TYPED MOLECULES MATCHED
    1 :        334         42 | C-C: [*:1]~[#6:2]~[#6:3]~[*:4]                     t0012: [#1:1]-[#6X4:2]-[#6X4:3]-[#8X2:4]               307 /     307 (100.000%)
    2 :        168         30 | C-O: [*:1]~[#6:2]~[#8:3]~[*:4]                     t0003: [a,A:1]-[#6X4:2]-[#8X2:3]-[!#1:4]               156 /     156 (100.000%)
    3 :        117         42 | 4808: [*:1]~[#6:2]~[#8:3]~[#1!X4:4]                t0002: [a,A:1]-[#6X4:2]-[#8X2:3]-[#1:4]                101 /     101 (100.000%)
    4 :         87         42 | 8090: [#6!H1:1]~[#6:2]~[#8:3]~[#1!X4:4]            t0006: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#1:4]                87 /     103 ( 84.466%)
    5 :        146         40 | 7751: [*:1]~[#6:2]~[#6:3]~[#6:4]                   t0001: [a,A:1]-[#6X4:2]-[#6X4:3]-[a,A:4]               146 /     146 (100.000%)
    6 :        131         37 | 1068: [#6!H3:1]~[#6:2]~[#6:3]~[#6:4]               t0007: [#6X4:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]             131 /     131 (100.000%)
    7 :        552         40 | 6774: [#1H0:1]~[#6:2]~[#6:3]~[#6:4]                t0005: [#1:1]-[#6X4:2]-[#6X4:3]-[#6X4:4]               552 /     552 (100.000%)
    8 :         66         30 | 8025: [#6:1]~[#6:2]~[#8:3]~[#6!H3:4]               t0008: [#6X4:1]-[#6X4:2]-[#8X2:3]-[#6X4:4]              66 /      66 (100.000%)
    9 :        574         42 | 1876: [#1:1]~[#6:2]~[#6:3]~[#1:4]                  t0004: [#1:1]-[#6X4:2]-[#6X4:3]-[#1:4]                 574 /     574 (100.000%)
TOTAL :       2175         42 |                                                       2120 /     2175 match (97.471 %)

Torsion type hierarchy:
    C-C ([*:1]~[#6:2]~[#6:3]~[*:4])
        7751 ([*:1]~[#6:2]~[#6:3]~[#6:4])
            1068 ([#6!H3:1]~[#6:2]~[#6:3]~[#6:4])
            6774 ([#1H0:1]~[#6:2]~[#6:3]~[#6:4])
        1876 ([#1:1]~[#6:2]~[#6:3]~[#1:4])
    C-O ([*:1]~[#6:2]~[#8:3]~[*:4])
        4808 ([*:1]~[#6:2]~[#8:3]~[#1!X4:4])
            8090 ([#6!H1:1]~[#6:2]~[#8:3]~[#1!X4:4])
        8025 ([#6:1]~[#6:2]~[#8:3]~[#6!H3:4])

```

## More smirky tests

The results from smirky tests get a bit bulky so we are not storing them on github. 
More extensive tests of both SMARTY and SMIRKY are available online:
Zanette, Camila et al. (2018), Supporting Information: Toward learned chemical perception of force field typing rules, v3, UC Irvine Dash, Dataset, [https://doi.org/10.7280/D1CD4C](https://doi.org/10.7280/D1CD4C)
