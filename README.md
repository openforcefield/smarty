[![Build Status](https://travis-ci.org/openforcefield/smarty.svg?branch=master)](https://travis-ci.org/openforcefield/smarty?branch=master)
[![DOI](https://zenodo.org/badge/60921138.svg)](https://zenodo.org/badge/latestdoi/60921138)

# `smarty`: Exploring Bayesian atom type sampling

This is a simple example of how Bayesian atom type sampling using reversible-jump Markov chain Monte Carlo (RJMCMC) [1] over SMARTS types might work.

All tools for implementation of the SMIRNOFF in OpenMM have been moved to the [openforcefield repository](https://github.com/openforcefield/openforcefield)

## Manifest

* `examples/` - some toy examples - look here to get started
* `smarty/` - simple toolkit illustrating the use of RJMCMC to sample over SMARTS-specified atom types and SMIRKS-specified bonded and non-bonded parameter types.
* `devtools/` - continuous integration and packaging scripts and utilities
* `oe_license.txt.enc` - encrypted OpenEye license for continuous integration testing
* `.travis.yml` - travis-ci continuous integration file
* `utilities/` - some utility functionality relating to the project, specifically testing the speed of ChemicalEnvironments for sampling in SMIRKY.

## Prerequisites

Install [miniconda](http://conda.pydata.org/miniconda.html) first. On `osx` with `bash`, this is:
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:${PATH}""
```

You must first install the OpenEye toolkit:
```
pip install -i https://pypi.anaconda.org/OpenEye/simple OpenEye-toolkits
```

You can then use conda to install smarty:
```
conda config --add channels omnia
conda install -c omnia smarty
```

## Installation

Install `smarty` from the `smarty/` directory with:
```bash
pip install .
```
If you modify the `smarty` source code (rather than the examples), reinstall with
```bash
pip install . --upgrade
```

## Documentation


## SMARTY atom type sampler

Check out the example in `examples/smarty/`:

Atom types are specified by SMARTS matches with corresponding parameter names.

First, we start with a number of initial "base types" which are essentially indestructible (often generic) atom types, specified in `atomtypes/basetypes.smarts`:
```
% atom types
[#1]    hydrogen
[#6]    carbon
[#7]    nitrogen
[#8]    oxygen
[#9]    fluorine
[#15]   phosphorous
[#16]   sulfur
[#17]   chlorine
[#35]   bromine
[#53]   iodine
```
Note that lines beginning with `%` are comment lines.

We also specify a number of starting types, "initial types" which can be the same or different from the base types. These follow the same format, and `atomtypes/basetypes.smarts` can be reused unless alternate behavior is desired (such as starting from more sophisticated initial types).

We have two sampler options for SMARTY which differ in how focused the sampling is. The original sampler samples over all elements/patterns at once, whereas the elemental sampler focuses on sampling only one specific element. The principle of sampling is the same; the only change is in which elements we sample over. To sample only over a single element, such as oxygen, for example, we use the elemental sampler to focus on that element.


### Generating New SMARTS patterns

There are two options for how to change SMARTS patterns when creating new atom types.
One is using combinatorial decorators (default) and the other is using simple decorators (`--decoratorbehavior=simple-decorators`). However, it should be noted that we have found the simple decorators insufficient at distinguishing atomtypes even for the most simple sets of molecules.

**Combinatorial Decorators**

The first option (combinatorial-decorator) attempt to create the new atomtype adding an Alpha or Beta substituent to a basetype or an atomtype.
This decorators are different from the simple-decorator option and do not have atom types or bond information on it.
The new decorators are listed in `AlkEthOH/atomtypes/new-decorators.smarts` and `parm@frosst/atomtypes/new-decorators.smarts`:

 ```
 % total connectivity
 X1             connections-1
 X2             connections-2
 X3             connections-3
 X4             connections-4
 % total-h-count
 H0             total-h-count-0
 H1             total-h-count-1
 H2             total-h-count-2
 H3             total-h-count-3
 % formal charge
 +0             neutral
 +1             cationic+1
 -1             anionic-1
 % aromatic/aliphatic
 a              aromatic
 A              aliphatic
 ```
Each decorator has a corresponding string token (no spaces allowed!) that is used to create human-readable versions of the corresponding atom types.

For example, we may find the atom type ```[#6]&H3``` which is `carbon total-h-count-3` for a C atom bonded to three hydrogens.

**Simple Decorators**
The second option (simple-decorators) attempts to split off a new atom type from a parent atom type by combining (via an "and" operator, `&`) the parent atom type with a "decorator".
The decorators are listed in `AlkEthOH/atomtypes/decorators.smarts` or `parm@frosst/atomtypes/decorators.smarts`:
```
% bond order
$([*]=[*])     double-bonded
$([*]#[*])     triple-bonded
$([*]:[*])     aromatic-bonded
% bonded to atoms
$(*~[#1])      hydrogen-adjacent
$(*~[#6])      carbon-adjacent
$(*~[#7])      nitrogen-adjacent
$(*~[#8])      oxygen-adjacent
$(*~[#9])      fluorine-adjacent
$(*~[#15])     phosphorous-adjacent
$(*~[#16])     sulfur-adjacent
$(*~[#17])     chlorine-adjacent
$(*~[#35])     bromine-adjacent
$(*~[#53])     iodine-adjacent
% degree
D1             degree-1
D2             degree-2
D3             degree-3
D4             degree-4
D5             degree-5
D6             degree-6
% valence
v1             valence-1
v2             valence-2
v3             valence-3
v4             valence-4
v5             valence-5
v6             valence-6
% total-h-count
H1             total-h-count-1
H2             total-h-count-2
H3             total-h-count-3
% aromatic/aliphatic
a              atomatic
A              aliphatic
```
This option also has the corresponding string tokens.

Newly proposed atom types are added to the end of the list.
After a new atom type is proposed, all molecules are reparameterized using the new set of atom types.
Atom type matching proceeds by trying to see if each SMARTS match can be applied working from top to bottom of the list.
This means the atom type list is hierarchical, with more general types appearing at the top of the list and more specific subtypes appearing at the bottom.

If a proposed type matches zero atoms, the RJMCMC move is rejected.

Currently, the acceptance criteria does not include the full Metropolis-Hastings acceptance criteria that would include the reverse probability.  This needs to be added in.

### Elemental Decomposition

The input option `--element` allows a user to specify which atoms types to sample based on atomic number. The default input is 0 (corresponding to no specified atomic number) and will attempt to match all atom types. If an element number is given (i.e. `--element=1` for hydrogen) only atoms with that atomic number are considered. Specifying an element number does not affect any other smarty behavior.

Finally, here is a complete list of input options for smarty. Under `usage` all bracketed parameters are optional.
```
Usage:     Sample over atom types, optionally attempting to match atom types in a reference typed set of molecules.

    usage: smarty --basetypes smartsfile --initialtypes smartsfile
            --decorators smartsfile --molecules molfile
            [--element atomicnumber --substitutions smartsfile --reference molfile
            --decoratorbehavior combinatorial-decorators/simple-decorators
            --iterations niterations --temperature temperature --trajectory trajectorfile
            --plot plotfile]

    example:
    python smarty --basetypes=atomtypes/basetypes.smarts --initialtypes=atomtypes/initialtypes.smarts \
            --decorators=atomtypes/decorators.smarts --substitutions=atomtypes/substitutions.smarts \
            --molecules=molecules/zinc-subset-tripos.mol2.gz --reference=molecules/zinc-subset-parm@frosst.mol2.gz \
            --iterations 1000 --temperature=0.1


Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -e ELEMENT, --element=ELEMENT
                        By default the element value is 0 corresponding to
                        sampling all atomtypes. If another atomic number is
                        specified only atoms with that atomic number are
                        sampled (i.e. --element=8 will only sample atomtypes
                        for oxygen atoms).
  -b BASETYPES, --basetypes=BASETYPES
                        Filename defining base or generic atom types as SMARTS
                        atom matches; these are indestructible and normally
                        are elemental atom types.
  -f BASETYPES, --initialtypes=BASETYPES
                        Filename defining initial (first) atom types as SMARTS
                        atom matches.
  -d DECORATORS, --decorators=DECORATORS
                        Filename defining decorator atom types as SMARTS atom
                        matches.
  -s SUBSTITUTIONS, --substitutions=SUBSTITUTIONS
                        Filename defining substitution definitions for SMARTS
                        atom matches (OPTIONAL).
  -r REFMOL, --reference=REFMOL
                        Reference typed molecules for computing likelihood
                        (must match same molecule and atom ordering in
                        molecules file) (OPTIONAL).
  -m MOLECULES, --molecules=MOLECULES
                        Small molecule set (in any OpenEye compatible file
                        format) containing 'dG(exp)' fields with experimental
                        hydration free energies.
  -i ITERATIONS, --iterations=ITERATIONS
                        MCMC iterations.
  -t TEMPERATURE, --temperature=TEMPERATURE
                        Effective temperature for Monte Carlo acceptance,
                        indicating fractional tolerance of mismatched atoms
                        (default: 0.1). If 0 is specified, will behave in a
                        greedy manner.
  -l TRAJECTORY_FILE, --trajectory=TRAJECTORY_FILE
                        Name for trajectory file output, trajectory saves only
                        changes to the list of 'atomtypes' for each iteration.
                        If the file already exists, it is overwritten.
  -p PLOT_FILE, --plot=PLOT_FILE
                        Name for output file of a plot of the score versus
                        time. If not specified, none will be written. If
                        provided, needs to use a file extension suitable for
                        matplotlib/pylab. Currently requires a trajectory file
                        to be written using -l or --trajectory.
  -x DECORATOR_BEHAVIOR, --decoratorbehavior=DECORATOR_BEHAVIOR
                        Choose between simple-decorators or combinatorial-
                        decorators (default = combinatorial-decorators).
```

---

## smirky

Check out examples in `examples/smirky/`:

This tool can sample any chemical environment type relevant to SMIRNOFFs, that is atoms, bonds, angles, and proper and improper torsions, one at a time
Scoring is analous to smarty (explained above), but uses a SMIRNOFF with existing parameters as a reference insteady of atomtyped molecules.

Input for this tool can require up to four different file types
* MOLECULES - any file that are readable in openeye, mol2, sdf, oeb, etc.
* ODDSFILES - File with the form "smarts     odds" for the different decorator or bond options
* SMARTS - .smarts file type with the form "smarts/smirks      label/typename"
* REFERENCE - a SMIRNOFF file with reference atoms, bonts, angles, torsions, and impropers

```
Usage:     Sample over fragment types (atoms, bonds, angles, torsions, or impropers)
    optionally attempting to match created types to an established SMIRNOFF.
    For all files left blank, they will be taken from this module's
    data/odds_files/ subdirectory.

    usage smirky --molecules molfile --typetag fragmentType
            [--atomORbases AtomORbaseFile --atomORdecors AtomORdecorFile
            --atomANDdecors AtomANDdecorFile --bondORbase BondORbaseFile
            --bondANDdecors BondANDdecorFile --atomIndexOdds AtomIndexFile
            --bondIndexOdds BondIndexFile --replacements substitutions
            --initialFragments initialFragments --SMIRNOFF referenceSMIRNOFF
            --temperature float --verbose verbose
            --iterations iterations --output outputFile]

    example:
    smirky -molecules AlkEthOH_test_filt1_ff.mol2 --typetag Angle



Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -m MOLECULES, --molecules=MOLECULES
                        Small molecule set (in any OpenEye compatible file
                        format) containing 'dG(exp)' fields with experimental
                        hydration free energies. This filename can also be an
                        option in this module's data/molecules sub-directory
  -T TYPETAG, --typetag=TYPETAG
                        type of fragment being sampled, options are 'VdW',
                        'Bond', 'Angle', 'Torsion', 'Improper'
  -e ODDFILES, --atomORbases=ODDFILES
                        Filename defining atom OR bases and associated
                        probabilities. These are combined with atom OR
                        decorators in SMIRKS, for example in
                        '[#6X4,#7X3;R2:2]' '#6' and '#7' are atom OR bases.
                        (OPTIONAL)
  -O ODDFILES, --atomORdecors=ODDFILES
                        Filename defining atom OR decorators and associated
                        probabilities. These are combined with atom bases in
                        SMIRKS, for example in '[#6X4,#7X3;R2:2]' 'X4' and
                        'X3' are ORdecorators. (OPTIONAL)
  -A ODDFILES, --atomANDdecors=ODDFILES
                        Filename defining atom AND decorators and associated
                        probabilities. These are added to the end of an atom's
                        SMIRKS, for example in '[#6X4,#7X3;R2:2]' 'R2' is an
                        AND decorator. (OPTIONAL)
  -o ODDFILES, --bondORbase=ODDFILES
                        Filename defining bond OR bases and their associated
                        probabilities. These are OR'd together to describe a
                        bond, for example in '[#6]-,=;@[#6]' '-' and '=' are
                        OR bases. (OPTIONAL)
  -a ODDFILES, --bondANDdecors=ODDFILES
                        Filename defining bond AND decorators and their
                        associated probabilities. These are AND'd to the end
                        of a bond, for example in '[#6]-,=;@[#7]' '@' is an
                        AND decorator.(OPTIONAL)
  -D ODDSFILE, --atomOddsFile=ODDSFILE
                        Filename defining atom descriptors and probabilities
                        with making changes to that kind of atom. Options for
                        descriptors are integers corresponding to that indexed
                        atom, 'Indexed', 'Unindexed', 'Alpha', 'Beta', 'All'.
                        (OPTIONAL)
  -d ODDSFILE, --bondOddsFile=ODDSFILE
                        Filename defining bond descriptors and probabilities
                        with making changes to that kind of bond. Options for
                        descriptors are integers corresponding to that indexed
                        bond, 'Indexed', 'Unindexed', 'Alpha', 'Beta', 'All'.
                        (OPTIONAL)
  -s SMARTS, --substitutions=SMARTS
                        Filename defining substitution definitions for SMARTS
                        atom matches. (OPTIONAL).
  -f SMARTS, --initialtypes=SMARTS
                        Filename defining initial (first) fragment types as
                        'SMIRKS    typename'. If this is left blank the
                        initial type will be a generic form of the given
                        fragment, for example '[*:1]~[*:2]' for a bond
                        (OPTIONAL)
  -r REFERENCE, --smirff=REFERENCE
                        Filename defining a SMIRNOFF force fielce used to
                        determine reference fragment types in provided set of
                        molecules. It may be an absolute file path, a path
                        relative to the current working directory, or a path
                        relative to this module's data subdirectory (for built
                        in force fields). (OPTIONAL)
  -i ITERATIONS, --iterations=ITERATIONS
                        MCMC iterations.
  -t TEMPERATURE, --temperature=TEMPERATURE
                        Effective temperature for Monte Carlo acceptance,
                        indicating fractional tolerance of mismatched atoms
                        (default: 0.1). If 0 is specified, will behave in a
                        greedy manner.
  -p OUTPUT, --output=OUTPUT
                        Filename base for output information. This same base
                        will be used for all output files created. If None
                        provided then it is set to 'typetag_temperature'
                        (OPTIONAL).
  -v VERBOSE, --verbose=VERBOSE
                        If True prints minimal information to the commandline
                        during iterations. (OPTIONAL)
``

## The SMIRNOFF force field format

The SMIRNOFF force field format is documented [here](https://github.com/openforcefield/openforcefield/blob/master/The-SMIRNOFF-force-field-format.md).
It was previously avaialbe in this repository, but has been moved.
SMIRNOFF99Frosst, a version of SMIRNOFF mirroring the parameters found in the parm@Frosst force field, is now housed in its own [repository](https://github.com/openforcefield/smirnoff99Frosst).
`forcefield.py` and other modules required to implement the SMIRNOFF format for simulations in OpenMM have also been moved. These scripts and examples on how to use them can be found at [openforcefield/openforcefield](https://github.com/openforcefield/openforcefield).

## References

[1] Green PJ. Reversible jump Markov chain Monte Carlo computation and Bayesian model determination. Biometrika 82:711, 1995.
http://dx.doi.org/10.1093/biomet/82.4.711
