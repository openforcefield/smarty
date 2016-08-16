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

from score_utils import load_trajectory
from score_utils import scores_vs_time
from environment import *
from forcefield import *
from utils import *

# ==============================================================================
# PRIVATE SUBROUTINES
# ==============================================================================

def _get_new_label(current, lowLim=1000, highLim=10000, maxIt = 1000000):
    """
    generates a random number that is between limits and not in current list

    Parameters
    -----------
    current: list of current labels
    lowlim: int, optional
        lower limit for random integers generated
    highlim: int, optional
        high limit for random integers generated
    maxIt: int, optional
        max iterations before None is returned for the new label

    Returns
    -------
    random integer between lowLim and highLim not in current list
    """
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
            each element can be an atomic number ('#1'),
            list of atomic numbers ('#1,#6,#7') or
            shorthand name that is in the replacements list
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
        self.types_with_no_matches = []

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
        self.emptyEnv = self.empty_environment(self.typetag)
        if initialtypes == None:
            self.envList = [copy.deepcopy(self.emptyEnv)]
        else:
            self.envList = [copy.deepcopy(initialtype) for initialtype in initialtypes]

        typeLabels = []
        for env in self.envList:
            if env.label == None:
                env.label = _get_new_label(typeLabels)
            typeLabels.append(env.label)

        self.baseTypes = copy.deepcopy(self.envList)
        # TODO: determine if base types should be a separate input from initial types

        # Compute total types being sampled
        self.total_types = 0.0
        smirks = self.emptyEnv.asSMIRKS()
        self.IndexDict = dict()
        for mol in self.molecules:
            matches  = self.get_SMIRKS_matches(mol, smirks)
            self.total_types += len(matches)
            smiles = OEMolToSmiles(mol)
            self.IndexDict[smiles] = matches

        # TODO: decide how to handle parent dictionary with environment objects.
        self.parents = dict()

        # Make typelist to fit method set up
        typelist = [[env.asSMIRKS(), env.label] for env in self.envList]

        # check that current smirks match all types in self.molecules
        if not self.check_typed_molecules(typelist):
            raise Exception("Initial types do not type all %s in the molecules" % self.typetag)

        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        if self.verbose: self.show_type_statistics(typelist, typecounts, molecule_typecounts)

        # Store elements without the [ ]
        self.elements = set()
        for element in elementList:
            e = element.replace('[','')
            e = e.replace(']','')
            for mol in self.molecules:
                matches = self.get_SMIRKS_matches(mol, '[%s:1]' % e)
                if len(matches) > 0:
                    self.elements.add(e)
        self.elements = list(self.elements)
        for e in self.elements:
            print(e)
        # Store reference molecules
        self.reference_types = []
        self.reference_typed_molecules = dict()
        self.reference_typename_dict = dict()
        if self.SMIRFF is not None:
            if self.verbose: print("Creating labeler from %s..." % self.SMIRFF)
            # get labeler for specified SMIRFF
            self.labeler = ForceField(get_data_filename(self.SMIRFF))
            # if verbose = True here it prints matches for every type for  every molecule!
            labels = self.labeler.labelMolecules(self.molecules, verbose = self.verbose)

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

            if not self.check_typed_molecules(self.reference_types):
                raise Exception("Reference types in SMIRFF (%s) do not type all %s in the molecules" % (self.SMIRFF, self.typetag))

            # Compute current type matches
            [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
            # Count  types.
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
            return 'NonbondedGenerator'
        if typetag.lower() == 'bond':
            return 'HarmonicBondGenerator'
        if typetag.lower() == 'angle':
            return 'HarmonicAngleGenerator'
        if typetag.lower() == 'torsion':
            return 'PeriodicTorsionGenerator'
        # TODO: what is the force word for impropers?
        if typetag.lower() == 'improper':
            return 'PeriodicTorsionGenerator'
        return None

    def empty_environment(self, typetag):
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
        typeDict = dict()
        for mol in self.molecules:
            smiles = OEMolToSmiles(mol)
            typeDict[smiles] = {}
            for [smirks, typename] in typelist:
                matches = self.get_SMIRKS_matches(mol, smirks)
                for match in matches:
                    typeDict[smiles][match] = typename

        return typeDict

    def check_SMIRK_matches(self, smirks):
        """
        checks that a give smirks matches at least one set of indices in the molecules
        """
        for mol in self.molecules:
            matches = self.get_SMIRKS_matches(mol, smirks)
            if len(matches) > 0:
                return True

        return False

    def check_typed_molecules(self, typelist):
        """
        given a typelist of [smirks, typename] it check that
        all types in each molecule can be typed
        """
        typed_molecules = self.get_typed_molecules(typelist)
        for smile, indicesList in self.IndexDict.items():
            for indices in indicesList:
                if indices not in typed_molecules[smile].keys():
                    return False
        return True

    def best_match_reference_types(self, typelist):
        """
        Determine best match for each parameter with reference types

        Parameters
        ----------
        typelist : list of list with form [smarts, typename]

        Returns
        -------
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of current types equivalently typed in reference molecule set.
        total_type_matches : int
            The total number of corresponding types in the reference molecule set.

        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.

        """
        if self.SMIRFF is None:
            if self.verbose: print('No reference SMIRFF specified, so skipping likelihood calculation.')
            return None, None # None for type_matches and total_type_matches

        # Create bipartite graph (U,V,E) matching current types U with reference types V via edges E with weights equal to number of types in common.
        if self.verbose: print('Creating graph matching current types with reference types...')
        initial_time = time.time()
        import networkx as nx
        graph = nx.Graph()

        # Get current types and reference types
        current_typenames = [ typename for (smirks, typename) in typelist ]
        reference_typenames = [ typename for (smirks, typename) in self.reference_types ]
        # check that current types are not in reference types
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
        Show pairing of current to reference types.

        Parameters
        ----------
        typelist: list of [smirks, typenames]
        type_matches : list of (current_typename, reference_typename, counts)
            List of type matches.

        Returns
        --------
        fraction_matched_types, the fractional count of matched types

        """
        print('%s type matches:' % self.typetag)
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

    def get_atom_OR_type(self):
        """
        returns a string that can be added to an Atom's ORtypes
        """
        # choose "element"
        element = random.choice(self.elements)
        # determine if you're going to decorate
        if random.random() < 0.5:
            return element
        else:
            decor = random.choice(self.ORdecorators)
            return element+decor

    def create_new_environment(self, env):
        """
        Given a parent environment type it creates a new child environment
        returns child environment type
        """
        new_env = copy.deepcopy(env)
        # get new label for it
        new_env.label = _get_new_label([e.label for e in self.envList])

        #Decide between changing an atom or a bond
        # TODO: determine how frequently to chose atom or bond
        if random.random() < 0.5:
            # pick a random atom
            atom = new_env.selectAtom()
            pick = random.randint(1,6)
            # TODO: determine how frequently each changetype can occur
            if pick == 1:
                # remove atom
                remove = new_env.removeAtom(atom)
                return new_env, remove

            elif pick == 2:
                # atom bound to current atom
                OR = self.get_atom_OR_type()
                new = new_env.addAtom(atom, newORtypes = [OR])
                if new == None:
                    return new_env, False
                else:
                    return new_env, True

            elif pick == 3:
                # add OR type to this atom
                OR = self.get_atom_OR_type()
                atom.addORtype(OR)
                return new_env, True

            elif pick == 4:
                # add AND type to this atom
                atom.addANDtype(random.choice(self.ANDdecorators))
                return new_env, True

            elif pick == 5:
                # Remove ORtype
                ORs = atom.getORtypes()
                if len(ORs) == 0:
                    return new_env, False
                else:
                    ORs.remove(random.choice(ORs))
                    atom.setORtypes(ORs)
                    return new_env, True
            else:
                # Remove ANDtype
                ANDs = atom.getANDtypes()
                if len(ANDs) == 0:
                    return new_env, False
                else:
                    ANDs.remove(random.choice(ANDs))
                    atom.setANDtypes(ANDs)
                    return new_env, True

        else:
            # pick a random bond
            atom1, atom2, bond = new_env.selectBond()
            pick = random.randint(1,4)
            if pick == 1:
                # Add ORtype
                bond.addORtype(random.choice(self.bondORset))
                return new_env, True

            elif pick == 2:
                # add ANDtype
                bond.addANDtype(random.choice(self.bondANDset))
                return new_env, True

            elif pick == 3:
                # Remove ORtype
                ORs = bond.getORtypes()
                if len(ORs) == 0:
                    return new_env, False
                else:
                    ORs.remove(random.choice(ORs))
                    bond.setORtypes(ORs)
                    return new_env, True
            else:
                ANDs = bond.getANDtypes()
                if len(ANDs) == 0:
                    return new_env, False
                else:
                    ANDs.remove(random.choice(ANDs))
                    bond.setANDtypes(ANDs)
                    return new_env, True

    def sample_types(self):
        """
        Perform one step of sampling the current set of chemical environments

        """
        # Copy current sets for proposal.
        proposed_envList = copy.deepcopy(self.envList)
        proposed_parents = copy.deepcopy(self.parents)
        ntypes = len(proposed_envList)

        # chose an environment from the list to focus on:
        env = random.choice(proposed_envList)

        # determine create or destroy:
        if random.random() < 0.5:
            # TODO: determine how frequently to destroy entire types
            if self.verbose: print("Attempting to destroy type %s : %s..." % (env.label, env.asSMIRKS()))

            # Reject deletion of (populated) base types as we want to retain
            # generics even if empty
            if env.label in [e.label for e in self.baseTypes]:
                if self.verbose: print("Destruction rejected for type %s because this is a generic type which was initially populated." % env.label)
                return False

            # Delete the type.
            proposed_envList.remove(env)

            # Try to type all molecules.
            typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            if not self.check_typed_molecules(typelist):
                if self.verbose: print("Typing failed; rejecting.")
                return False

            # update proposed parent dictionary
            # TODO: update parent dictionary

        else: # create new type from chosen environment
            new_env, changed = self.create_new_environment(env)

            # TODO: determine if this is allowed
            # continue getting a new environment as long as no change was made
            while not changed:
                new_env, changed = self.create_new_environment(env)

            if self.verbose: print("Attempting to create new subtype: '%s' (%s) from parent type '%s' (%s)" % (new_env.label, new_env.asSMIRKS(), env.label, env.asSMIRKS()))

            # Check the SMIRKS for new_env is valid
            qmol = OEQMol()
            smirks = new_env.asSMIRKS()
            if self.replacements is not None:
                smirks = OESmartsLexReplace(smirks, self.replacements)
            if not OEParseSmarts(qmol, smirks):
                if self.verbose: print("Type '%s' (%s) is invalid; rejecting." % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if new_env is already in types with no matches
            if new_env.asSMIRKS() in self.types_with_no_matches:
                if self.verbose: print("Type '%s' (%s) unused in dataset; rejecting." % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if proposed type is already in set.
            if new_env.asSMIRKS() in [e.asSMIRKS for e in self.envList]:
                if self.verbose: print("Type '%s' (%s) already exists; rejecting to avoid duplication." % (new_env.label, new_env.asSMIRKS()))
                return False

            # add new type to proposed list
            proposed_envList.append(new_env)
            proposed_typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            # Compute updated statistics
            [proposed_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_typelist)

            # Reject if new type matches nothing
            if proposed_typecounts[new_env.label] == 0:
                if self.verbose: print("Type '%s' (%s) unused in dataset; rejecting." % (new_env.label, new_env.asSMIRKS()))
                self.types_with_no_matches.append(new_env.asSMIRKS())
                return False

            # Reject if parent type is now unused (UNLESS parent is a base type)
            if env.label not in [e.label for e in self.baseTypes]:
                # parent not in base types
                if proposed_typecounts[env.label] == 0:
                    if self.verbose: print("Parent type '%s' (%s) now unused in dataset; rejecting." % (env.label, env.asSMIRKS()))
                    return False

            # updated proposed parent dictionary
            # TODO: update parent dictionary

        if self.verbose: print('Proposal is valid...')

        # Accept automatically if no set was provided
        if self.SMIRFF is None:
            self.envList = proposed_envList
            self.parents = proposed_parents
            return True

        # Compute effective temperature
        if self.temperature == 0.0:
            effective_temperature = 1
        else:
            effective_temperature = (self.total_types * self.temperature)

        # Compute likelihood for accept/reject
        typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]

        (proposed_type_matches, proposed_total_type_matches) = self.best_match_reference_types(typelist)
        log_P_accept = (proposed_total_type_matches - self.total_type_matches) / effective_temperature
        print('Proposal score: %d >> %d : log_P_accept = %.5e' % (self.total_type_matches, proposed_total_type_matches, log_P_accept))
        if (log_P_accept > 0.0) or (numpy.random.uniform() < numpy.exp(log_P_accept)):
            # Change accepted
            self.envList = proposed_envList
            self.parents = proposed_parents
            self.type_matches = proposed_type_matches
            self.total_type_matches = proposed_total_type_matches
            return True

        else: # Change not accpted
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
        # Zero type counts by typename and molecule.
        typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        typed_molecules = self.get_typed_molecules(typelist)
        # Count number of indice sets with each type.
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
            Best correspondence between current and reference types, along with number of fragment types equivalently typed in reference molecule set.
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
        Collects typecount information for a csv "trajectory" output file

        Parameters
        -----------
        typelist - list of lists with form [smirks, typename]
        typecounts (dict) - number of matches for each fragment type
        molecule_typecounds (dict) - number of molecules that contain each fragment type
        type_matches : list of tuples (current_typelabel, reference_typelabel, counts)
            Best correspondence between current and reference types, along with number of fragment types equivalently typed in reference molecule set.
        """
        if type_matches is not None:
            reference_type_info = dict()
            for (current_typename, reference_typename, count) in type_matches:
                reference_type_info[current_typename] = (reference_typename, count)

        index = 1
        output = []
        # Print counts
        # INDEX, SMARTS, PARENT INDEX, REF TYPE, MATCHES, MOLECULES, FRACTION, OUT of, PERCENTAGE
        for [smarts, typename] in typelist:
            if type_matches is not None:
                (reference_typename, reference_count) = reference_type_info[typename]
                if reference_typename is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    # Save output
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, reference_typename, typecounts[typename], molecule_typecounts[typename], reference_count, reference_total))
                else:
                    output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', typecounts[typename], molecule_typecounts[typename], 0, 0))

            else:
                output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, 0, 0, 'NONE', typecounts[typename], molecule_typecounts[typename], 0, 0))
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
        Run sampler for the specified number of iterations.

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
        fraction_matched : float
            fraction of total types matched successfully at end of run

        """
        self.traj = []
        for iteration in range(niterations):
            if self.verbose:
                print("Iteration %d / %d" % (iteration, niterations))

            accepted = self.sample_types()
            typelist = [[env.asSMIRKS(), env.label] for env in self.envList]
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

                # Compute type statistics on molecules.
                self.show_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
                print('')

                # TODO: figure out how to handle parent dictionary with chemical environments
                # Print parent tree as it is now.
                print("%s type hierarchy: will go HERE" % self.typetag)

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
        typelist = [ [env.asSMIRKS(), env.label] for env in self.envList]
        [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        fraction_matched = self.show_type_matches(typelist, self.type_matches)

        # If verbose print parent tree:
        if self.verbose:
            # TODO: update to monitor parent/child hierarchy
            print("%s type hierarchy: will go HERE" % self.typetag)
            #self.print_parent_tree(roots, '\t')
        return fraction_matched
