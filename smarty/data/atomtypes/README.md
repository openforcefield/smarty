# Atom type SMARTS components

## Formats
smarts files are used as input for the smarty sampler
there are a variety of types, detailed below. All follow
the same general format.

Comments beginning with `%` are ignored throughout the file.
Each line has the format
```
<SMARTS> <typename>
```
where `<SMARTS>` is an [OpenEye SMARTS string](https://docs.eyesopen.com/toolkits/cpp/oechemtk/SMARTS.html) and `<typename>` is a human-readable typename associated with that atom type.

Atom type definitions are hierarchical, with the last match in the file taking precedence over earlier matches.

### Initial and Base types

These are both used to initialize the smarty sampler.
`basetypes` are considered more generic. 
These are the atomtypes used to create new atomtypes.
See the file `basetypes.smarts`.

`initial` types can be more complex 
for example the files
`initialtypes.smarts` or `initiali\_AlkEthOH.smarts`

Best practices should have base and initial types that are listed from most to
least general

### Simple and Combinatorial Decorators

A `decorators` file contains a list of SMARTS

In smarty, when using simple decorators, the new atomtypes are created only
by ANDing the decorator SMARTS component to the parent atomtype (using the `&` operator).
The human-readable `<decoratorname>` is appended (with a space) to the parent name to keep a human-readable annotation of the proposed child atom type.


Example simple decorators are in *`decorators.smarts`* and are typically more complicated as they must include all 
ways of generating new atomtypes

Combinatorial decorators use a more complex set of rules to generate new SMARTS strings. 
In this case, bonded atoms are found in the basetypes, so only "non-bonding decorators" need to be 
in the decorator file. 
For exampl see *`new-decorators.smarts`* 

### Substitutions

It is often convenient to define various tokens that are substituted for more sophisticated SMARTS expressions.

For example, we could define some elemental substitutions along with some substitutions for halogens:
```
% elements
[#9]    fluorine
[#17]   chlorine
[#35]   bromine
[#53]   iodine

% halogens
[$smallhals,$largehals]     halogen
[$fluorine,$chlorine]       smallhals
[$bromine,$iodine]          largehals
```

The [`OESmartsLexReplace`](http://docs.eyesopen.com/toolkits/python/oechemtk/OEChemFunctions/OESmartsLexReplace.html) function is used to implement these replacements.

## Manifest
* `basetypes.smarts` - basetypes file with elemental atom types - this is a good choice to begin with
* `initial.smarts` - basetypes file with more sophisticated atom types
* `initial\_AlkEthOH.smarts` - the "answer" SMARTS strings for the AlkEthOH molecule set
* `decorators.smarts` - `decorators` file with a variety of decorators
* `decorators-simple.smarts` - minimal `decorators` file for testing
* `new-decorators.smarts` - decorators file without bond information (new modular framework)
* `substitutions.smarts` - minimal `substitutions` file
