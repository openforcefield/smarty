#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smirk.py

Example illustrating a scheme to create and destroy bonds, angles, torsions and impropers
automatically using chemical environments to produce SMIRKS

AUTHORS

John Chodera <john.chodera@choderalab.org>, Memorial Sloan Kettering Cancer Center.
Caitlin Bannan <bannanc@uci.edu>, UC Irvine
Additional contributions from the Mobley lab, UC Irvine, including David Mobley, and Camila Zanette
and from the Codera lab, Josh Fass

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import networkx

import time

from . import AtomTyper
from score_utils import load_trajectory
from score_utils import scores_vs_time
from environment import *
from forcefield_labeler import *
from utils import *

# ==============================================================================
# PRIVATE SUBROUTINES
# ==============================================================================

def _getNewLabel(current, lowLim=1000, highLim=10000, maxIt = 1000000):
    # TODO: write doc string
    label = random.randint(lowLim, highLim)
    it = 0
    while label in current:
        label = random.randint(lowLim, highLim)
        it += 1
        if it > maxIt:
            return None
    return label

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class TypeSampler(object):
    """
    SMIRKS sampler for atoms, bonds, angles, torsions, and impropers.
    """
    def __init__(self, molecules, typetag, elementList, ORdecorators,
            ANDdecorators, replacements = None,  initialtypes = None,
            SMIRFF = None, temperature = 0.1, verbose = False):
        """
        Initialize a fragment type sampler
        For VdW, Bond, Angle, Torsion, or Improper

        Parameters
        -----------
        molecules : list of OEmol objects, required
            List of molecules used for sampling
        typetag : string, required
            Must 'Bond', 'Angle', 'Torsion', 'Improper', or 'VdW'
            'VdW' is for single labeled atom
        elementList : list of strings, required
            List of elements or combined elements
        ORdecorators: list of strings, required
            List of decorators that can be combined directly with an atom
            for example: for [#6X4, #8X2] 'X4' and 'X2' are ORdecorators
        ANDdecorators: list of strings, required
            List of decorators that are AND'd to the end of an atom
            for example: in [#6,#7,#8;H0;+0] 'H0' and '+0' are ANDdecorators
        replacements: list of the form [short hand, smarts], optional
        initialtypes: list of chemical environments, optional
            if None, the typetag is used to make an empty environment, such as [*:1]~[*:2] for a bond
        SMIRFF: string, optional
            file with the SMIRFF you wish to compare fragment typing with
        temperature : float, optional, default=0.1
            Temperature for Monte Carlo acceptance/rejection
        verbose : bool, optional, default=False
            If True, verbose output will be printed.

        Notes
        -----

        """
        # Save properties that remain unchanged
        self.verbose = verbose
        self.ORdecorators = ORdecorators
        self.ANDdecorators = ANDdecorators
        self.temperature = temperature
        self.SMIRFF = SMIRFF
        self.replacements = replacements

        if not typetag.lower() in ['bond', 'angle', 'torsion','improper','vdw']:
            raise Exception("Error typetag %s is not recognized, please use 'Bond', 'Angle', 'Torsion', 'Improper', or 'VdW' ")
        self.typetag = typetag
        self.forcetype = self.get_force_type(self.typetag)

        # Save bond list to use throughout
        self.bondORset = ['-', '=', '#', ':', '!-', '!=', '!#', '!:']
        self.bondANDset = ['@', '!@']

        # get molecules and add explicit hydrogens
        self.molecules = copy.deepcopy(molecules)
        for mol in self.molecules:
            OEAddExplicitHydrogens(mol)

        # if no initialtypes specified make empty bond
        self.emptyEnv = self.emptyEnvironment(self.typetag)
        if initialtypes == None:
            self.envList = [copy.deepcopy(self.emptyEnv)]
        else:
            self.envList = [copy.deepcopy(initialtype) for initialtype in initialtypes]

        self.typeLabels = []
        for env in self.envList:
            if env.label == None:
                env.label = _getNewLabel(self.typeLabels)
            self.typeLabels.append(env.label)

        # TODO: decide how to handle parent dictionary with environment objects.

        # Make typelist to fit method set up
        typelist = [[env.asSMIRKS(), env.label] for env in self.envList]
        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        if self.verbose: self.show_type_statistics(typelist, typecounts, molecule_typecounts)

        # Store elements in the molecules
        self.elements = []

        for element in elementList:
            # remove [ and ]
            e = element.replace(']','')
            e = e.replace('[','')
            self.elements.append(e)


        # Compute total types being sampled
        self.total_types = 0.0
        smirks = self.emptyEnv.asSMIRKS()
        for mol in self.molecules:
            matches = self.get_SMIRKS_matches(mol, smirks)
            self.total_types += len(matches)

        # Store reference molecules
        # TODO: update how to handle reference molecules
        self.reference_types = []
        self.current_atom_matches = None
        self.reference_typed_molecules = dict()
        self.reference_typename_dict = dict()
        if self.SMIRFF is not None:
            # get labeler for specified SMIRFF
            self.labeler = ForceField_labeler(get_data_filename(self.SMIRFF))
            # if verbose = True here it prints matches for every type for  every molecule!
            labels = self.labeler.labelMolecules(self.molecules, verbose = False)

            # save the type we are considering
            self.ref_labels = [l[self.forcetype] for l in labels]

            # get smiles to key reference typed molecule dictionary
            smiles = [OEMolToSmiles(mol) for mol in molecules]
            # Extract list of reference SMIRKS types present in molecules
            for idx, label_set in enumerate(self.ref_labels):
                smile = smiles[idx]
                self.reference_typed_molecules[smile] = {}
                for (indices, pid, smirks) in label_set:
                    self.reference_typename_dict[pid] = smirks
                    self.reference_typed_molecules[smile][tuple(indices)] = pid

            self.reference_types = [[smirks, pid] for pid, smirks in self.reference_typename_dict.items()]
            # Compute current atom matches
            [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
            # Count atom types.
            self.reference_type_counts = { pid : 0 for (smirks, pid) in self.reference_types }
            for label_set in self.ref_labels:
                for (atom_indices, pid, smirks) in label_set:
                    self.reference_type_counts[pid] += 1
            if self.verbose: self.show_type_statistics(typelist, typecounts, molecule_typecounts, self.type_matches)
        return

    def get_force_type(self, typetag):
        """
        Uses typetag to get the Force type key word used to read the SMIRFF file

        Parameters
        ----------
        typetag: string, required
            'vdw', 'bond', 'angle', 'torsion', 'improper'
            indicates the type of system being sampled
        """
        if typetag.lower() == 'vdw':
            return 'NonbondedForce'
        if typetag.lower() == 'bond':
            return 'HarmonicBondForce'
        if typetag.lower() == 'angle':
            return 'HarmonicAngleForce'
        if typetag.lower() == 'torsion':
            return 'PeriodicTorsionForce'
        # TODO: what is the force word for impropers?
        if typetag.lower() == 'improper':
            return 'Improper'
        return None

    def emptyEnvironment(self, typetag):
        """
        Returns an empty atom, bond, angle, torsion or improper

        Parameters
        -----------
        typetag: string, required
            'vdw', 'bond', 'angle', 'torsion', 'improper'
            indicates the type of system being sampled
        """
        if typetag.lower() == 'vdw':
            return AtomChemicalEnvironment()
        if typetag.lower() == 'bond':
            return BondChemicalEnvironment()
        if typetag.lower() == 'angle':
            return AngleChemicalEnvironment()
        if typetag.lower() == 'torsion':
            return TorsionChemicalEnvironment()
        if typetag.lower() == 'improper':
            return ImproperChemicalEnvironment()
        return None

    def get_SMIRKS_matches(self, mol, smirks):
        """
        Gets atom indices for a smirks string in a given molecule

        Parameters
        ----------
        mol : an OpenEye OEMol object
        smirks : a string for the SMIRKS string being parsed

        Returns
        --------
        matches: list of tuples
            atom indices for labeled atom in the smirks
        """
        if self.replacements is not None:
            smirks = OELexReplace(smirks, self.replacements)

        qmol = OEQMol()
        if not OEParseSmarts(qmol, smirks):
            raise Exception("Error parsing SMIRKS %s" % smirks)

        # ValenceDict was written to handle symmetry in SMIRKS
        matches = ValenceDict()

        # then require non-unique matches
        unique = False
        ss = OESubSearch(qmol)

        for match in ss.Match(mol, unique):
            indices = dict()
            for ma in match.GetAtoms():
                patMap = ma.pattern.GetMapIdx()
                # if patMap == 0, then it's an unidexed atom
                if patMap != 0:
                    indices[patMap-1] = ma.target.GetIdx()

            indices = [indices[idx] for idx in range(len(indices))]
            matches[indices] = ''

        return matches.keys()

    def get_typed_molecules(self, typelist):
        """
        Creates a dictionary assigning a typename
        for each set of atom indices in each molecule

        Parameters
        ----------
        typelist: list of tuples in the form (smirks, typename)

        Returns
        -------
        typeDict: embedded dictionary
            keys: SMILES string for each molecule
                keys: tuple of indices assigned a parameter type
        """
        # TODO: write doc string
        typeDict = dict()
        for mol in self.molecules:
            smiles = OEMolToSmiles(mol)
            typeDict[smiles] = {}
            for [smirks, typename] in typelist:
                matches = self.get_SMIRKS_matches(mol, smirks)
                for match in matches:
                    typeDict[smiles][match] = typename

        return typeDict

    def best_match_reference_types(self, typelist):
        """
        Determine best match for each parameter with reference atom types

        Parameters
        ----------
        typelist : list of list with form [smarts, typename]

        Returns
        -------
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference atomtypes, along with number of atoms equivalently typed in reference molecule set.
        total_atom_type_matches : int
            The total number of correspondingly typed atoms in the reference molecule set.

        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.

        """
        if self.SMIRFF is None:
            if self.verbose: print('No reference SMIRFF specified, so skipping likelihood calculation.')
            return None

        # Create bipartite graph (U,V,E) matching current atom types U with reference atom types V via edges E with weights equal to number of atoms typed in common.
        if self.verbose: print('Creating graph matching current types with reference types...')
        initial_time = time.time()
        import networkx as nx
        graph = nx.Graph()

        # Get current atomtypes and reference atom types
        current_typenames = [ typename for (smirks, typename) in typelist ]
        reference_typenames = [ typename for (smirks, typename) in self.reference_types ]
        # check that current atom types are not in reference atom types
        if set(current_typenames) & set(reference_typenames):
            raise Exception("Current and reference type names must be unique")
        # Add current types
        for typename in current_typenames:
            graph.add_node(typename, bipartite=0)
        # add reference types
        for typename in reference_typenames:
            graph.add_node(typename, bipartite=1)
        # Add edges.
        types_in_common = dict()
        for current_typename in current_typenames:
            for reference_typename in reference_typenames:
                types_in_common[(current_typename,reference_typename)] = 0

        current_typed_molecules = self.get_typed_molecules(typelist)

        for smile, indexDict in current_typed_molecules.items():
            for indices, current_typename in indexDict.items():
                reference_typename = self.reference_typed_molecules[smile][indices]
                types_in_common[(current_typename, reference_typename)] += 1

        for current_typename in current_typenames:
            for reference_typename in reference_typenames:
                weight = types_in_common[(current_typename,reference_typename)]
                graph.add_edge(current_typename, reference_typename, weight=weight)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Graph creation took %.3f s' % elapsed_time)

        # Compute maximum match
        if self.verbose: print('Computing maximum weight match...')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Maximum weight match took %.3f s' % elapsed_time)

        # Compute match dictionary and total number of matches.
        type_matches = list()
        total_type_matches = 0
        for current_typename in current_typenames:
            if current_typename in mate:
                reference_typename = mate[current_typename]
                counts = graph[current_typename][reference_typename]['weight']
                total_type_matches += counts
                type_matches.append( (current_typename, reference_typename, counts) )
            else:
                type_matches.append( (current_typename, None, None) )

        # Report on matches
        if self.verbose:
            print("PROPOSED:")
            self.show_type_matches(typelist, type_matches)

        return (type_matches, total_type_matches)

    def show_type_matches(self, typelist, type_matches):
        """
        Show pairing of current to reference atom types.

        Parameters
        ----------
        typelist: list of [smirks, typenames]
        type_matches : list of (current_typename, reference_typename, counts)
            List of atom type matches.

        Returns
        --------
        fraction_matched_atoms, the fractional count of matched atoms

        """
        print('Atom type matches:')
        total_type_matches = 0
        current_dict = dict()
        for (current_typename, reference_typename, counts) in type_matches:
            current_dict[current_typename] = (reference_typename, counts)

        for [current_smirks, current_typename] in typelist:
            (reference_typename, counts) = current_dict[current_typename]
            current_combo = "%s: %s" % (current_typename, current_smirks)
            if reference_typename is not None:
                reference_smirks = self.reference_typename_dict[reference_typename]
                reference_combo = "%s: %s" % (reference_typename, reference_smirks)
                print('%-64s matches %64s: %8d %-10s types matched' % (current_combo, reference_combo, counts, self.typetag))
                total_type_matches += counts
            else:
                print('%-64s no match' % (current_combo))

        fraction_matched = float(total_type_matches) / float(self.total_types)
        print('%d / %d total %ss match (%.3f %%)' % (total_type_matches, self.total_types, self.typetag, fraction_matched * 100))

        return fraction_matched

    def PickAType(self, envList):
        """
        Takes a list of chemical environments and returns one random entry
        """
        type_index = random.randint(0, len(envList)-1)
        return envList[type_index]

    def sample_types(self):
        """
        Perform one step of sampling the current set of chemical environments

        """
        # Copy current atomtypes for proposal.
        proposed_envList = copy.deepcopy(self.envList)
        ntypes = len(proposed_envList)

        valid_proposal = True

        if random.random() < 0.5:
            # Pick an atom type to destroy.
            atomtype_index = random.randint(0, natomtypes-1)
            (atomtype, typename) = proposed_atomtypes[atomtype_index]
            if self.verbose: print("Attempting to destroy atom type %s : %s..." % (atomtype, typename))
            # Reject deletion of (populated) base types as we want to retain
            # generics even if empty
            if [atomtype, typename] in self.used_basetypes:
                if self.verbose: print("Destruction rejected for atom type %s because this is a generic type which was initially populated." % atomtype )
                return False

            # Delete the atomtype.
            proposed_atomtypes.remove([atomtype, typename])

            # update proposed parent dictionary
            for parent, children in proposed_parents.items():
                if atomtype in [at for [at, tn] in children]:
                    children += proposed_parents[atomtype]
                    children.remove([atomtype, typename])

            del proposed_parents[atomtype]

            # Try to type all molecules.
            try:
                self.type_molecules(proposed_atomtypes, proposed_molecules)
            except AtomTyper.TypingException as e:
                # Reject since typing failed.
                if self.verbose: print("Typing failed; rejecting.")
                valid_proposal = False
        else:
            if self.decorator_behavior == 'simple-decorators':
                # Pick an atomtype to subtype.
                atomtype_index = random.randint(0, natomtypes-1)
                # Pick a decorator to add.
                decorator_index = random.randint(0, ndecorators-1)
                # Create new atomtype to insert by appending decorator with 'and' operator.
                (atomtype, atomtype_typename) = self.atomtypes[atomtype_index]
                (decorator, decorator_typename) = self.decorators[decorator_index]
                result = re.match('\[(.+)\]', atomtype)
                proposed_atomtype = '[' + result.groups(1)[0] + '&' + decorator + ']'
                proposed_typename = atomtype_typename + ' ' + decorator_typename
                if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atomtype, atomtype_typename, decorator, decorator_typename, proposed_atomtype, proposed_typename))

                # Update proposed parent dictionary
                proposed_parents[atomtype].append([proposed_atomtype, proposed_typename])
                # Hack to make naming consistent with below
                atom1smarts, atom1typename = atomtype, atomtype_typename

            else:
                # combinatorial-decorators
                nbondset = len(self.bondset)
                # Pick an atomtype
                atom1type = self.PickAnAtom(self.unmatched_atomtypes)
                atom1smarts, atom1typename = atom1type
                # Check if we need to add an alfa or beta substituent
                if self.HasAlpha(atom1type):
                    # Has alpha
                    bondset_index = random.randint(0, nbondset-1)
                    atom2type = self.PickAnAtom(self.used_basetypes)
                    if random.random() < 0.5 or atom1type[0][2] == '1': # Add Beta Substituent Atom randomly or when it is Hydrogen
                        proposed_atomtype, proposed_typename = self.AddBetaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type)
                    else: # Add another Alpha Substituent if it is not a Hydrogen
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type, first_alpha = False)
                    if self.verbose: print("Attempting to create new subtype: -> '%s' (%s)" % (proposed_atomtype, proposed_typename))
                else:
                    # Has no alpha
                    if random.random() < 0.5:
                        # Add a no-bond decorator
                        decorator_index = random.randint(0, ndecorators-1)
                        decorator = self.decorators[decorator_index]
                        proposed_atomtype, proposed_typename = self.AtomDecorator(atom1type, decorator)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], decorator[0], decorator[1], proposed_atomtype, proposed_typename))
                    else:
                        bondset_index = random.randint(0, nbondset-1)
                        atom2type = self.PickAnAtom(self.used_basetypes)
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, self.bondset[bondset_index], atom2type, first_alpha = True)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], proposed_atomtype, proposed_typename))


                # Update proposed parent dictionary
                proposed_parents[atom1type[0]].append([proposed_atomtype, proposed_typename])

            proposed_parents[proposed_atomtype] = []

            # Check that we haven't already determined this atom type isn't matched in the dataset.
            if proposed_atomtype in self.atomtypes_with_no_matches:
                if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                return False

            # Check if proposed atomtype is already in set.
            existing_atomtypes = set()
            for (a, b) in self.atomtypes:
                existing_atomtypes.add(a)
            if proposed_atomtype in existing_atomtypes:
                if self.verbose: print("Atom type already exists; rejecting to avoid duplication.")
                valid_proposal = False

            # Check for valid proposal before proceeding.
            if not valid_proposal:
                return False

            # Insert atomtype immediately after.
            proposed_atomtypes.insert(natomtypes, [proposed_atomtype, proposed_typename]) # Insert in the end (hierarchy issue)
            # Try to type all molecules.
            try:
                # Type molecules.
                self.type_molecules(proposed_atomtypes, proposed_molecules)
                # Compute updated statistics.
                [proposed_atom_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_atomtypes, proposed_molecules)
                # Reject if new type is unused.
                if (proposed_atom_typecounts[proposed_typename] == 0):
                    # Reject because new type is unused in dataset.
                    if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                    valid_proposal = False
                    # Store this atomtype to speed up future rejections
                    self.atomtypes_with_no_matches.add(proposed_atomtype)
                # Reject if parent type is now unused, UNLESS it is a base type
                if (proposed_atom_typecounts[atom1typename] == 0) and (atom1smarts not in self.basetypes_smarts):
                    # Reject because new type is unused in dataset.
                    if self.verbose: print("Parent type '%s' (%s) now unused in dataset; rejecting." % (atom1smarts, atom1typename))
                    valid_proposal = False
            except AtomTyper.TypingException as e:
                print("Exception: %s" % str(e))
                # Reject since typing failed.
                if self.verbose: print("Typing failed for one or more molecules using proposed atomtypes; rejecting.")
                valid_proposal = False

        # Check for valid proposal
        if not valid_proposal:
            return False

        if self.verbose: print('Proposal is valid...')

        # Accept automatically if no reference molecules
        accept = False
        if self.reference_typed_molecules is None:
            accept = True
        else:
            # Compute effective temperature
            if self.temperature == 0.0:
                effective_temperature = 1
            else:
                effective_temperature = (self.total_atoms * self.temperature)

            # Compute likelihood for accept/reject
            (proposed_atom_type_matches, proposed_total_atom_type_matches) = self.best_match_reference_types(proposed_atomtypes, proposed_molecules)
            log_P_accept = (proposed_total_atom_type_matches - self.total_atom_type_matches) / effective_temperature
            print('Proposal score: %d >> %d : log_P_accept = %.5e' % (self.total_atom_type_matches, proposed_total_atom_type_matches, log_P_accept))
            if (log_P_accept > 0.0) or (numpy.random.uniform() < numpy.exp(log_P_accept)):
                accept = True

        # Accept or reject
        if accept:
            self.atomtypes = proposed_atomtypes
            self.molecules = proposed_molecules
            self.parents = proposed_parents
            self.atom_type_matches = proposed_atom_type_matches
            self.total_atom_type_matches = proposed_total_atom_type_matches
            return True
        else:
            return False

    def compute_type_statistics(self, typelist):
        """
        Compute statistics for numnber of molecules assigned each type.

        Parameters
        ----------
        typelist: list of lists with the form [smarts, typename]

        Returns
        -------
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounds (dict) - number of molecules that contain each fragment type

        """
        # Zero type counts by atom and molecule.
        typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        typed_molecules = self.get_typed_molecules(typelist)
        # Count number of atoms with each type.
        for molecule, indexDict in typed_molecules.items():
            typenames_in_molecule = set()
            for indices, typename in indexDict.items():
                typecounts[typename] += 1
                typenames_in_molecule.add(typename)
            for typename in typenames_in_molecule:
                molecule_typecounts[typename] += 1

        return (typecounts, molecule_typecounts)

    def show_type_statistics(self, typelist, typecounts, molecule_typecounts, type_matches=None):
        """
        Print type statistics.

        Parameters
        -----------
        typelist - list of lists with form [smirks, typename]
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounds (dict) - number of molecules that contain each fragment type
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference atomtypes, along with number of atoms equivalently typed in reference molecule set.
        """
        index = 1
        ntypes = 0

        if type_matches is not None:
            reference_type_info = dict()
            for (current_typename, reference_typename, count) in type_matches:
                reference_type_info[current_typename] = (reference_typename, count)

        # Print header
        if type_matches is not None:
            print "%5s   %10sS %10s   %-50s %-50s %30s" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS', 'REF TYPE: SMIRKS', 'FRACTION OF REF TYPED MOLECULES MATCHED')
        else:
            print "%5s   %10sS %10s   %-50s" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS')

        # Print counts
        for [smarts, typename] in typelist:
            current_combo = "%s: %s" % (typename, smarts)
            if type_matches is not None:
                (reference_typename, reference_count) = reference_type_info[typename]
                if reference_typename is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    reference_fraction = float(reference_count) / float(reference_total)
                    reference_combo = "%s: %s" % (reference_typename, self.reference_typename_dict[reference_typename])
                    print "%5d : %10d %10d | %-50s %-50s %7d / %7d (%7.3f%%)" % (index, typecounts[typename], molecule_typecounts[typename], current_combo, reference_combo, reference_count, reference_total, reference_fraction*100)
                else:
                    print "%5d : %10d %10d | %-50s" % (index, typecounts[typename], molecule_typecounts[typename], current_combo)
            else:
                print "%5d : %10d %10d | %-50s" % (index, typecounts[typename], molecule_typecounts[typename], current_combo)

            ntypes += typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)

        if type_matches is not None:
            print "%5s : %10d %10d |  %15s %32s %8d / %8d match (%.3f %%)" % ('TOTAL', ntypes, nmolecules, '', '', self.total_type_matches, self.total_types, (float(self.total_type_matches) / float(self.total_types)) * 100)
        else:
            print "%5s : %10d %10d" % ('TOTAL', ntypes, nmolecules)

        return

    def save_type_statistics(self, typelist, typecounts, molecule_typecounts, type_matches=None):
        """
        Collects atom typecount information for a csv "trajectory" output file

        Parameters
        -----------
        typelist - list of lists with form [smirks, typename]
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounds (dict) - number of molecules that contain each fragment type
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference atomtypes, along with number of atoms equivalently typed in reference molecule set.
        """
        if type_matches is not None:
            reference_type_info = dict()
            for (current_typename, reference_typename, count) in type_matches:
                reference_type_info[current_typename] = (reference_atomtype, count)

        index = 1
        output = []
        # Print counts
        # INDEX, SMARTS, PARENT INDEX, REF TYPE, MATCHES, MOLECULES, FRACTION, OUT of, PERCENTAGE
        for [smarts, typename] in typelist:
            if type_matches is not None:
                (reference_atomtype, reference_count) = reference_type_info[typename]
                if reference_atomtype is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    # Save output
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, reference_typename, typecounts[typename], molecule_typecounts[typename], reference_count, reference_total))
                else:
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))

            else:
                output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))
            index += 1
        return output

    def print_parent_tree(self, roots, start=''):
        """
        Recursively prints the parent tree.

        Parameters
        ----------
        roots = list of smarts strings to print
        """
        for r in roots:
            print("%s%s" % (start, r))
            if r in self.parents.keys():
                new_roots = [smart for [smart, name] in self.parents[r]]
                self.print_parent_tree(new_roots, start+'\t')

    def run(self, niterations, trajFile=None, plotFile=None):
        """
        Run atomtype sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations
        trajFile : str, optional, default=None
            Output trajectory filename
        plotFile : str, optional, default=None
            Filename for output of plot of score versus time

        Returns
        ----------
        fraction_matched_atoms : float
            fraction of total atoms matched successfully at end of run

        """
        self.traj = []
        for iteration in range(niterations):
            if self.verbose:
                print("Iteration %d / %d" % (iteration, niterations))

            accepted = self.sample_types()
            typelist = [[env.asSMIRKS(), env.label] for env in self.envList]
            [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
            [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)

            if trajFile is not None:
                # Get data as list of csv strings
                lines = self.save_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
                # Add lines to trajectory with iteration number:
                for l in lines:
                    self.traj.append('%i,%s \n' % (iteration, l))

            if self.verbose:
                if accepted:
                    print('Accepted.')
                else:
                    print('Rejected.')

                # Compute atomtype statistics on molecules.
                self.show_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
                print('')

                # TODO: figure out how to handle parent dictionary with chemical environments
                # Print parent tree as it is now.
                #roots = self.parents.keys()
                # Remove keys from roots if they are children
                #for parent, children in self.parents.items():
                #    child_smarts = [smarts for [smarts, name] in children]
                #    for child in child_smarts:
                #        if child in roots:
                #            roots.remove(child)

                #print("Atom type hierarchy:")
                #self.print_parent_tree(roots, '\t')

        if trajFile is not None:
            # make "trajectory" file
            if os.path.isfile(trajFile):
                print("trajectory file already exists, it was overwritten")
            f = open(trajFile, 'w')
            start = ['Iteration,Index,Smarts,ParNum,ParentParNum,RefType,Matches,Molecules,FractionMatched,Denominator\n']
            f.writelines(start + self.traj)
            f.close()

            # Get/print some stats on trajectory
            # Load timeseries
            timeseries = load_trajectory( trajFile )
            time_fractions = scores_vs_time( timeseries )
            print("Maximum score achieved: %.2f" % max(time_fractions['all']))

        #Compute final type stats
        typelist = [ [env.asSMIRKS, env.label] for env in self.envList]
        [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        fraction_matched_atoms = self.show_type_matches(typelist, self.type_matches)

        # If verbose print parent tree:
        if self.verbose:
            # TODO: update to monitor parent/child hierarchy
            print("Need to add printing for parent dictionary back in")
            #roots = self.parents.keys()
            # Remove keys from roots if they are children
            #for parent, children in self.parents.items():
            #    child_smarts = [smarts for [smarts, name] in children]
            #    for child in child_smarts:
            #        if child in roots:
            #            roots.remove(child)

            #print("Atom type hierarchy:")
            #self.print_parent_tree(roots, '\t')
        return fraction_matched_atoms
