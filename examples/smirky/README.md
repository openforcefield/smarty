# smirky sampling of Torsions

This is an example of how to use smirky, a command line tool for sampling chemical perception of Bonds, Angles, proper or improper Torsions, or van der Waal parameters. Default smirky behaivor only requires two inputs, this is an example of all input options into smirky. Files in this directory include:

* `atom_OR_bases.smarts` - element numbers that form the base of atoms, such as `"#6"` and their associated odds
* `atom_OR_decorators.smarts` - decorators and associated odds that are combined with element numbers such as `X4` in `"[#6X3,#7]"`
* `atom_AND_decorators.smarts` - decorators and associated odds for patterns that are "AND'd" to the end of an atom for example `r5` in `"[#6X4,#7X3;r5]"`
* `bond_OR_bases.smarts` - bond bases and their associated odds, that is '-', '=', ':', or '#' typically
* `bond_AND_decorators.smarts` - bond decorators that can be "AND'd" in a bond, such as '@' in `"[#6r6]-,:;@[#7r6]"`
* `atom_odds_forTorsions.smarts` - keywords or indices for atoms in torsions and odds of making changes to them
* `bond_odds_forTorsions.smarts` - keywords or indices for bonds in torsions and odds of making changes to them
* `initial_Torsions.smarts` - SMIRKS patterns for initial patterns
* `substitutions.smarts` - SMIRKS patterns and the short hand they can be replaced with

In this example, the following was called from the command line:

```
smirky --molecules AlkEtOH_test_filt1_ff.mol2 \
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
    --smirff forcefield/Frosst_AlkEtOH.ffxml \
    --iteratorsion 1000 \
    --temperature 0.001 \ 
    --verbose True \
    --output output
```

The following files were created:
* output.log - detailed log of each iteration, changes made and if it was accepted or rejected
* output.csv - a "trajectory" file that describes the torsions at each iteration
* output.pdf - plot showing the overall score vs iteration
* output_results.smarts - smarts file showing the file SMIRKS and their matched results

Final score was 97.5%

The results from smirky tests get a bit bulky so we are no storing them on github. 
We maintain a Google Drive directory, 
(smirky_testing)[https://drive.google.com/drive/folders/0BwF2-3puCvfESmFGNTQ4SGlsYnc?usp=sharing]
that is public if you want to see more examples. 
Please note the Google Drive directory is still a work in progress, more documentation
is on its way. 
