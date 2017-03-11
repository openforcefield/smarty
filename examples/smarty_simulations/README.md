# Example application of SMARTY atom type sampler to recover parm99 typing of alkanes, ethers, and alcohols

These are example outputs for a variety of smarty uses. Each example is listed below with the associated command line call.
Each example has the three output files with the title of the example as the name:
* `*.csv` - example trajectory file, a csv file that is readable with the `score\_util.py` methods
* `*.log` - stored commandline output for that simulation
* `*.pdf` - plot showing the score verses iteration for the simulation

These are only examples of how to use smarty. All input files are those included in the smarty package 
available at `smart/data/`, the utility here allows those files to be used in simulations.

## AlkEthOH

Typical smarty behavior with the AlkEthOH molecule set
with combinatorial decorators and sampling all atoms

```
smarty --basetypes atomtypes/basetypes.smarts \
    --initialtypes atomtypes/basetypes.smarts \
    --decorators atomtypes/new-decorators.smarts \
    --molecules AlkEthOH_test_filt1_tripos.mol2 \
    --reference AlkEthOH_test_filt1_ff.mol2 \
    --iterations 1000 \
    --temperature 0.01 \
    --trajectory AlkEthOH.csv \
    --plot AlkEthOH.pdf >> AlkEthOH.log
```

**Example Output** 
this output shows how smarty is used to sample atomtypes 
and compared to the parm@frosst typed reference molecules

##### Initializing smarty:
```
Loading molecules from '/Users/bannanc/anaconda/lib/python2.7/site-packages/smarty/data/molecules/AlkEthOH_test_filt1_tripos.mol2'...
42 molecules read
0.006 s elapsed
Loading molecules from '/Users/bannanc/anaconda/lib/python2.7/site-packages/smarty/data/molecules/AlkEthOH_test_filt1_ff.mol2'...
42 molecules read
0.006 s elapsed
Sampling all atomtypes
```
Store bond types that are used in these molecules
```
USED BOND TYPES:
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :        803         42 |                                                           singly                                -
    2 :          0          0 |                                                           doubly                                =
    3 :          0          0 |                                                           triply                                #
    4 :          0          0 |                                                         aromatic                                :
TOTAL :        803         42
```
Type molecules with base types and store those with matches
```
MATCHED BASETYPES:
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :        464         42 |                                                       c_hydrogen                             [#1]
    2 :        232         42 |                                                         c_carbon                             [#6]
    3 :          0          0 |                                                       c_nitrogen                             [#7]
    4 :        107         42 |                                                         c_oxygen                             [#8]
    5 :          0          0 |                                                       c_fluorine                             [#9]
    6 :          0          0 |                                                    c_phosphorous                            [#15]
    7 :          0          0 |                                                         c_sulfur                            [#16]
    8 :          0          0 |                                                       c_chlorine                            [#17]
    9 :          0          0 |                                                       c_selenium                            [#34]
   10 :          0          0 |                                                        c_bromine                            [#35]
   11 :          0          0 |                                                         c_iodine                            [#53]
TOTAL :        803         42
Removing basetype '[#7]' ('c_nitrogen'), which is unused.
Removing basetype '[#9]' ('c_fluorine'), which is unused.
Removing basetype '[#15]' ('c_phosphorous'), which is unused.
Removing basetype '[#16]' ('c_sulfur'), which is unused.
Removing basetype '[#17]' ('c_chlorine'), which is unused.
Removing basetype '[#34]' ('c_selenium'), which is unused.
Removing basetype '[#35]' ('c_bromine'), which is unused.
Removing basetype '[#53]' ('c_iodine'), which is unused.
```
Type molecules with initial types and store the ones that are used
```
MATCHED INITIAL TYPES:
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS
    1 :        464         42 |                                                       c_hydrogen                             [#1]
    2 :        232         42 |                                                         c_carbon                             [#6]
    3 :          0          0 |                                                       c_nitrogen                             [#7]
    4 :        107         42 |                                                         c_oxygen                             [#8]
    5 :          0          0 |                                                       c_fluorine                             [#9]
    6 :          0          0 |                                                    c_phosphorous                            [#15]
    7 :          0          0 |                                                         c_sulfur                            [#16]
    8 :          0          0 |                                                       c_chlorine                            [#17]
    9 :          0          0 |                                                       c_selenium                            [#34]
   10 :          0          0 |                                                        c_bromine                            [#35]
   11 :          0          0 |                                                         c_iodine                            [#53]
TOTAL :        803         42
Removing initial atom type '[#7]', as it matches no atoms
Removing initial atom type '[#9]', as it matches no atoms
Removing initial atom type '[#15]', as it matches no atoms
Removing initial atom type '[#16]', as it matches no atoms
Removing initial atom type '[#17]', as it matches no atoms
Removing initial atom type '[#34]', as it matches no atoms
Removing initial atom type '[#35]', as it matches no atoms
Removing initial atom type '[#53]', as it matches no atoms
```
Use bi-partite scoring sceme to score current atomtypes against reference
```
Creating graph matching current atom types with reference atom types...
Graph creation took 0.008 s
Computing maximum weight match...
Maximum weight match took 0.001 s
```
Initial types and which reference they are paired with and initial score (67.746 %)
```
Atom type matches:
c_hydrogen                                                       matches       HC :      244 atoms matched
c_carbon                                                         matches       CT :      232 atoms matched
c_oxygen                                                         matches       OH :       68 atoms matched
544 / 803 total atoms match (67.746 %)
```
##### Example move in chemical space
```
Iteration 16 / 1000
Attempting to create new subtype: '[#1]' (c_hydrogen) -> '[#1$(*~[#6])]' (c_hydrogen any c_carbon )
Proposal is valid...
```
Score proposed atomtypes against reference
```
Creating graph matching current atom types with reference atom types...
Graph creation took 0.007 s
Computing maximum weight match...
Maximum weight match took 0.001 s
PROPOSED:
Atom type matches:
c_hydrogen                                                       matches       HO :       68 atoms matched
c_carbon                                                         matches       CT :      232 atoms matched
c_oxygen                                                         matches       OH :       68 atoms matched
c_hydrogen any c_carbon                                          matches       HC :      244 atoms matched
612 / 803 total atoms match (76.214 %)
```
##### Accepting or Rejecting a Move
A move that leads to an increased score will always be accepted.
A move with a decrease has a probability of being accepted depending on the temperature.
A 0.0 temperature will lead lead to a complete optimizer where only moves leading to an increased score are accepted,
however these can get stuck in local optima. By using a non-zero temperature we allow more moves to be accepted
and a larger chemical space to be explored.
```
Proposal score: 544 >> 612 : log_P_accept = 8.46824e+00
Accepted.
```
Score by reference atomtype
```
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :         68         42 |                                                       c_hydrogen                             [#1]       HO               68 /               68 (100.000%)
    2 :        232         42 |                                                         c_carbon                             [#6]       CT              232 /              232 (100.000%)
    3 :        107         42 |                                                         c_oxygen                             [#8]       OH               68 /               68 (100.000%)
    4 :        396         42 |                                         c_hydrogen any c_carbon                     [#1$(*~[#6])]       HC              244 /              244 (100.000%)
TOTAL :        803         42 |                                                                                                         612 /      803 match (76.214 %)
```
Atomtype hierarchy shows which parent type a child descends from
```
Atom type hierarchy:
    [#6]
    [#8]
    [#1]
        [#1$(*~[#6])]
```
##### Final iteration of this simulation
```
Iteration 999 / 1000
Attempting to destroy atom type [#6] : c_carbon...
Destruction rejected for atom type [#6] because this is a generic type which was initially populated.
Rejected.
INDEX        ATOMS  MOLECULES                                                          TYPE NAME                           SMARTS REF TYPE        FRACTION OF REF TYPED MOLECULES MATCHED
    1 :        291         42 |                                                       c_hydrogen                             [#1]       HC              244 /              244 (100.000%)
    2 :        232         42 |                                                         c_carbon                             [#6]       CT              232 /              232 (100.000%)
    3 :         39         30 |                                                         c_oxygen                             [#8]       OS               39 /               39 (100.000%)
    4 :         68         42 |                                         c_hydrogen any c_oxygen                     [#1$(*~[#8])]       HO               68 /               68 (100.000%)
    5 :         27         21 | c_hydrogen any c_carbon any c_carbon (any c_oxygen) (singly c_oxygen) [#1$(*~[#6](-[#8])(~[#8])~[#6])]       H2               27 /               33 ( 81.818%)
    6 :         78         25 | c_hydrogen any c_carbon any c_carbon (any c_oxygen) (singly c_hydrogen) [#1$(*~[#6](-[#1])(~[#8])~[#6])]       H1               78 /              116 ( 67.241%)
    7 :         68         42 |                                         c_oxygen any c_hydrogen                     [#8$(*~[#1])]       OH               68 /               68 (100.000%)
TOTAL :        803         42 |                                                                                                         756 /      803 match (94.147 %)

Atom type hierarchy:
    [#1]
        [#1$(*~[#8])]
        [#1$(*~[#6](-[#8])(~[#8])~[#6])]
        [#1$(*~[#6](-[#1])(~[#8])~[#6])]
    [#8]
        [#8$(*~[#1])]
    [#6]
Maximum score achieved: 0.99
```

## Hydrogen

This is an example of how to implement the elemental sampler for smarty
you only need to add the `--element` option. In this case instead of considering
all atoms, we only sample atom types for hydrogen. 
This allows for more efficient testing of the smarty tool as we can 
focus on the chemical perception sampling around one element. 
In the AlkEthOH, there is only 1 carbon and 2 oxygens, so the 5 hydrogen types
are the best example of this behavior.  

```
smarty --element 1 \
    --basetypes atomtypes/basetypes.smarts \
    --initialtypes atomtypes/basetypes.smarts \
    --decorators atomtypes/new-decorators.smarts \
    --molecules AlkEthOH_test_filt1_tripos.mol2 \
    --reference AlkEthOH_test_filt1_ff.mol2 \
    --iterations 1000 \
    --temperature 0.01 \
    --trajectory Hydrogen.csv \
    --plot Hydrogen.pdf >> Hydrogen.log
```

## Simple-Decorators 

With the simple decorator option new atomtypes are generated by ANDing 
decorator SMARTS patterns to the end of a parent atomtype.
This method is not capable of even getting the complexity in the AlkEthOH
molecule set as it does not allow for beta substitution from the primary atom.
 
```
smarty --basetypes atomtypes/basetypes.smarts \
    --initialtypes atomtypes/basetypes.smarts \
    --decorators atomtypes/decorators.smarts \
    --substitutions atomtypes/replacements.smarts \
    --molecules AlkEthOH_test_filt1_tripos.mol2 \
    --reference AlkEthOH_test_filt1_ff.mol2 \
    --iterations 1000 \
    --temperature 0.01 \
    --trajectory Simple-decorators.csv \
    --plot Simple-decorators.pdf \
    --decoratorbehavior simple-decorators >> Simple-decorators.log
```

## More smarty tests
We have done more extensive testing of this tool, but the results are 
a bit bulky to keep on GitHub. We maintain a public (Google Drive Directory)[https://drive.google.com/drive/folders/0BwF2-3puCvfEeWNuNnlsTm1CTlU?usp=sharing]
with these results. Please note it is a work in progress so documentation is on going.
