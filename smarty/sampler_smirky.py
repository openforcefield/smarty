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
import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import networkx
import time

from  smarty.environment import *
from smarty.forcefield import *
from smarty.utils import *

# Currently not used
from smarty.score_utils import load_trajectory
from smarty.score_utils import scores_vs_time

import numpy
from numpy import random
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

def _OddsToWeights(choices):
    """
    Takes a set of properties and their odds and converts to weights

    choices: tuple of two (properties, odds)
        properties: list of anything
        odds: list of numbers or None
            if None all properties are equally probable

    returns tuple of two (properties, probabilities)
    """
    if len(choices) != 2:
        raise Exception("Error: _OddsToWeights(choices) expects 2-tuple (properties, odds)")

    props = choices[0]

    if choices[1] is None:
        current_odds = numpy.ones(len(props))
    else:
        if len(choices[0]) != len(choices[1]):
            raise Exception("Error: choices with odds should consist of a 2-tuple ([properties], [odds]) where properties and odds are equal lengths")
        current_odds = numpy.array(choices[1], dtype = np.float)

    if len(props) != len(current_odds):
        raise Exception("Error: _OddsToWeights takes a 2-tuple with lists of equal length")
    weights = current_odds / np.sum(current_odds)
    return (props, weights)

def _PickFromWeightedChoices(choices):
    """
    Given a set of properties adn their relative odds

    choices: tuple of two (properties, /weights)
        properties: list of anything
        weights: list of probabilities
    """
    if len(choices) != 2:
        raise Exception("Error: PickFromWeightedChoices expects a 2-tuple, %s does not apply" % str(choices))
    if choices[1] is None or np.sum(choices[1]) != 1.0:
        choices = _OddsToWeights(choices)
    indices = range(len(choices[0]))
    pickIndex = random.choice(indices, p=choices[1])
    return choices[0][pickIndex], choices[1][pickIndex]

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class FragmentSampler(object):
    """
    SMIRKS sampler for atoms, bonds, angles, torsions, and impropers.
    """
    def __init__(self, molecules, typetag, AtomORbases, AtomORdecorators,
            AtomANDdecorators, BondORbases, BondANDdecorators,
            AtomIndexOdds = None, BondIndexOdds = None,
            replacements = None,  initialtypes = None,
            SMIRFF = None, temperature = 0.1, outputFile = None):
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

        The following parameters come in the form of tuples of two
        with ( [features], [odds or probabilities])
        if [odds or probabilities] = None all are treated equally
        odds like [1,4] will be converted to probabilities [.2, .8]
        ------------------------------------------------
        AtomORbases : list of strings and their probabilities, required
            each element can be an atomic number ('#1') or
            shorthand name that is in the replacements list
        AtomORdecorators: list of strings and their probabilities, required
            List of decorators that can be combined directly with an atom
            for example: for [#6X4, #8X2] 'X4' and 'X2' are ORdecorators
        AtomANDdecorators: list of strings and their probabilities, required
            List of decorators that are AND'd to the end of an atom
            for example: in [#6,#7,#8;H0;+0] 'H0' and '+0' are ANDdecorators
        BondORbases : list of strings and their probabilities, required
            each bond type and the probability of using it
        BondANDdecorators : of strings and their probabilitis, required
            bond decorator options such as [ ('@', 1), ('@', 1) ]
        AtomIndexOdds and BondIndexOdds have list of atom and bond indices
            options for property lists: atom indices (integer),
            or string descriptor from : indexed, unindexed, alpha, beta
            if None all atoms and bonds will be treated equally
        ------------------------------------------------
        replacements: list of the form [short hand, smarts], optional
        initialtypes: initial typelist in form [smirks, typename], optional
            if None, the typetag is used to make an empty environment, such as [*:1]~[*:2] for a bond
        SMIRFF: string, optional
            file with the SMIRFF you wish to compare fragment typing with
        temperature : float, optional, default=0.1
            Temperature for Monte Carlo acceptance/rejection
        outputFile: string, optional, default = typetag_temperature
            output name base for log files and trajectory files

        Notes
        -----

        """
        # Save properties that remain unchanged
        self.AtomORdecorators = _OddsToWeights(AtomORdecorators)
        self.AtomANDdecorators = _OddsToWeights(AtomANDdecorators)
        self.BondORbases = _OddsToWeights(BondORbases)
        self.BondORdecorators = ( [''], [1.])
        self.BondANDdecorators = _OddsToWeights(BondANDdecorators)
        self.temperature = temperature
        self.SMIRFF = SMIRFF
        self.replacements = replacements
        self.types_with_no_matches = []

        if outputFile == None:
            self.output = "%s_%.2e" % (typetag, temperature)
        else:
            self.output = outputFile

        self.log = open("%s.log" % self.output, 'w')

        if AtomIndexOdds is None:
            AtomIndexOdds = [['all'], [1]]
        if BondIndexOdds is None:
            BondIndexOdds = [ ['all'], [1]]

        self.AtomIndexOdds = _OddsToWeights(AtomIndexOdds)
        self.BondIndexOdds = _OddsToWeights(BondIndexOdds)

        if not typetag.lower() in ['bond', 'angle', 'torsion','improper','vdw']:
            raise Exception("Error typetag %s is not recognized, please use 'Bond', 'Angle', 'Torsion', 'Improper', or 'VdW' ")
        self.typetag = typetag
        self.forcetype = self.get_force_type(self.typetag)

        # get molecules and add explicit hydrogens
        self.molecules = copy.deepcopy(molecules)
        for mol in self.molecules:
            OEAddExplicitHydrogens(mol)

        # if no initialtypes specified make empty bond
        self.emptyEnv = self._makeEnvironments(self.typetag, None)[0]
        self.envList = self._makeEnvironments(self.typetag, initialtypes)
        self.baseTypes = copy.deepcopy(self.envList)

        # Compute total types being sampled
        self.total_types = 0.0
        empty_typelist = [[self.emptyEnv.asSMIRKS(), 'empty']]
        [empty_counts, empty_molecule_counts] = self.compute_type_statistics(empty_typelist)
        self.total_types = empty_counts['empty']
        self.IndexDict = self.get_typed_molecules(empty_typelist)

        # This might be better as a graph with unlabeled nodes where the nodes are environment objects?
        self.parents = dict()
        for env in self.envList:
            self.parents[env.label] = dict()
            self.parents[env.label]['children'] = list()
            self.parents[env.label]['SMIRKS'] = env.asSMIRKS()
            self.parents[env.label]['parent'] = None

        # Make typelist to fit method set up
        typelist = [[env.asSMIRKS(), env.label] for env in self.envList]

        # check that current smirks match all types in self.molecules
        if not self.check_typed_molecules(typelist):
            raise Exception("Initial types do not type all %s in the molecules" % self.typetag)

        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        self.write_type_statistics(typelist, typecounts, molecule_typecounts)

        # Store elements without the [ ]
        self.AtomORbases = ( [], [])
        elementList = _OddsToWeights(AtomORbases)
        elementList = zip(elementList[0], elementList[1])
        for (element, prob) in elementList:
            e = element.replace('[','')
            e = e.replace(']','')
            # Check if element is used in any molecule
            for mol in self.molecules:
                useE = False
                matches = self.get_SMIRKS_matches(mol, '[%s:1]' % e)
                if len(matches) > 0:
                    useE = True
                    break
            # Save the used ones only
            if useE:
                self.AtomORbases[0].append(e)
                self.AtomORbases[1].append(prob)
            else:
                self.log.write("removing unused element (%s) from list\n" % element)

        # Store reference molecules
        self.reference_types = []
        self.reference_typed_molecules = dict()
        self.type_matches = None
        self.total_type_matches = None
        self.reference_typename_dict = dict()
        if self.SMIRFF is not None:
            self.log.write("Creating labeler from %s...\n" % self.SMIRFF)
            # get labeler for specified SMIRFF
            #self.labeler = ForceField(get_data_filename(self.SMIRFF))
            self.labeler = ForceField(self.SMIRFF)
            labels = self.labeler.labelMolecules(self.molecules, verbose=False)

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
                raise Exception("Reference types in SMIRFF (%s) do not type all %ss in the molecules" % (self.SMIRFF, self.typetag))

            # Compute current type matches
            [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
            # Count  types.
            self.reference_type_counts = { pid : 0 for (smirks, pid) in self.reference_types }
            for label_set in self.ref_labels:
                for (atom_indices, pid, smirks) in label_set:
                    self.reference_type_counts[pid] += 1
            self.write_type_statistics(typelist, typecounts, molecule_typecounts, self.type_matches)
        return

    def _makeEnvironments(self, typetag, smirksList):
        """
        Given a typetag and list of SMIRKS strings
        returns a list of chemical environment objects
        of the correct type
        """
        if typetag.lower() == 'vdw':
            chemEnv = AtomChemicalEnvironment
        elif typetag.lower() == 'bond':
            chemEnv = BondChemicalEnvironment
        elif typetag.lower() == 'angle':
            chemEnv = AngleChemicalEnvironment
        elif typetag.lower() == 'torsion':
            chemEnv = TorsionChemicalEnvironment
        elif typetag.lower() == 'improper':
            chemEnv = ImproperChemicalEnvironment
        else:
            return None

        envList = list()
        if smirksList is None:
            return [chemEnv(None, 0, self.replacements)]

        for smirks, typename in smirksList:
            envList.append(chemEnv(smirks, typename, self.replacements))

        return envList

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
        if typetag.lower() == 'torsion' or typetag.lower() == 'improper':
            return 'PeriodicTorsionGenerator'
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
            smirks = OESmartsLexReplace(smirks, self.replacements)

        qmol = OEQMol()
        if not OEParseSmarts(qmol, smirks):
            raise Exception("Error parsing SMIRKS %s" % smirks)

        # Impropers have different symmetry from other fragments
        if self.typetag.lower() == 'improper':
            matches = dict()
        else:
            # ValenceDict allow for symmetric fragments
            # for example bond (1,2) is identical to (2,1)
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
            self.log.write('No reference SMIRFF specified, so skipping likelihood calculation.\n')
            return None, None # None for type_matches and total_type_matches

        # Create bipartite graph (U,V,E) matching current types U with reference types V via edges E with weights equal to number of types in common.
        self.log.write('Creating graph matching current types with reference types...\n')
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
        self.log.write('Graph creation took %.3f s\n' % elapsed_time)

        # Compute maximum match
        self.log.write('Computing maximum weight match...\n')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        self.log.write('Maximum weight match took %.3f s\n' % elapsed_time)

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

        self.log.write("PROPOSED:\n")
        self.write_type_matches(typelist, type_matches)

        return (type_matches, total_type_matches)

    def write_type_matches(self, typelist, type_matches):
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
        self.log.write('%s type matches:\n' % self.typetag)
        if type_matches is None:
            return 0
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
                self.log.write('%-64s matches %64s: %8d %-10s types matched\n' % (current_combo, reference_combo, counts, self.typetag))
                total_type_matches += counts
            else:
                self.log.write('%-64s no match\n' % (current_combo))

        fraction_matched = float(total_type_matches) / float(self.total_types)
        self.log.write('%d / %d total %ss match (%.3f %%)\n' % (total_type_matches, self.total_types, self.typetag, fraction_matched * 100))

        return fraction_matched

    def pick_an_atom(self, env):
        """
        Uses AtomIndexOdds to chose an atom
        returns an Atom object and its associated probability
        """
        choices = copy.deepcopy(self.AtomIndexOdds)
        atom = None
        while atom is None:
            descriptor,prob = _PickFromWeightedChoices(choices)
            # TODO: selectAtom has random built in, account for those odds?
            atom = env.selectAtom(descriptor)
            choices[1][choices[0].index(descriptor)] = 0
        return atom, prob

    def pick_a_bond(self, env):
        """
        Uses BondIndexOdds to chose a bond
        returns a Bond object and its associated probability
        """
        choices = copy.deepcopy(self.BondIndexOdds)
        bond = None
        while bond is None:
            descriptor, prob = _PickFromWeightedChoices(choices)
            # TODO: selectBond also uses random, not accounting for those odds
            bond = env.selectBond(descriptor)
            choices[1][choices[0].index(descriptor)] = 0
        return bond, prob

    def add_swap_delete(self, current, new, new_prob, probabilities = None):
        """
        Makes a change to a current list of properties

        Parameters
        -----------
        current: current list for that property
        new: new property for the list
        new_prob: probability for picking that new property
        probabilities: probability or Odds list for
            [add property, swap property, delete a property]

        Returns
        --------
        newList, probability
            of creating it
            returns None, 0 if no change was made
        """
        opts = [1,2,3]
        if probabilities is None:
            probabilities = [1,1,1]

        if len(current) == 0:
            # if nothing in the list then adding is the only option
            move = 1
            move_prob = 1
        elif new in current: # add and swap are not an option
            move = 3
            move_prob = 1
        else:
            move, move_prob = _PickFromWeightedChoices( (opts, probabilities))

        newList = copy.deepcopy(current)

        if move == 2 or move == 3: # swap or delete
            remove, remove_prob = _PickFromWeightedChoices((current, None))
            newList.remove(remove)
            move_prob *= remove_prob

        if move == 1 or move == 2: # add or swap property
            newList.append(new)
            move_prob *= new_prob

        return newList, move_prob

    def add_atom(self, env, atom):
        """
        Adds an atom bonded to the current atom
        New atoms have 1 ORbase and no ORdecorators or ANDdecorators

        Returns probability of creating this atom
        """
        #choose ORbase
        base, prob = _PickFromWeightedChoices( self.AtomORbases )
        new_atom = env.addAtom(atom, newORtypes = [(base, [])])
        return prob

    def isremoveable(self,env, atom):
        """
        Returns true if the atom can be removed
        """
        if env.isIndexed(atom):
            return False # indexed atom isn't removeable
        if len(atom.getANDtypes() ) > 0:
            return False # AND decorators not added to new atoms
        ORs = atom.getORtypes()
        if len(ORs) > 1: # only 1 OR type on new atoms
            return False
        if len(ORs) == 1 and len(ORs[0][1]) > 0:
            return False # only 1 OR type with no OR decorators

        bonds = env.getBonds(atom)
        if len(bonds) > 1:
            return False # atom has more than 1 neighbor
        bond = bonds[0]
        if len(bond.getANDtypes() ) > 0:
            return False
        if len(bond.getORtypes() ) > 0:
            return False

        return True

    def change_atom(self, env, atom):
        """
        Makes changes to the provided Atom object
            in the provided environment
        returns the probability of making this change
            if 0 no change was made
        """
        # assign options and probabilities:
        # start with OR base change opt = 3
        opts = [3]
        probs = [3]
        # if atom is not beta, add atom is allowed
        if not env.isBeta(atom):
            opts.append(1)
            probs.append(5)
        # check if removeable
        if self.isremoveable(env, atom):
            opts.append(0)
            probs.append(1)
        # check for optional AND decorators
        if not (len(self.AtomANDdecorators[0]) <=1 and self.AtomANDdecorators[0][0] ==''):
            opts.append(2)
            probs.append(1)
        # check for existing OR types, for OR decorator changes
        if len(atom.getORtypes()) > 0:
            if not (len(self.AtomORdecorators[0]) <= 1 and self.AtomORdecorators[0][0] == ''):
                opts.append(4)
                probs.append(5)

        # Choose a move
        move, move_prob = _PickFromWeightedChoices( (opts, probs))
        if move == 0: # remove Atom
            removed = env.removeAtom(atom, False)
            change_prob = 1.
        elif move == 1: # add Atom
            change_prob = self.add_atom(env, atom)
        elif move == 2: # Change ANDdecorators
            change_prob = self.change_ANDdecorators(atom, self.AtomANDdecorators)
        elif move == 3: # Change ORbases
            change_prob = self.change_ORbase(atom, self.AtomORbases, self.AtomORdecorators)

        else: # move == 4 Change ORdecorators
            change_prob = self.change_ORdecorator(atom, self.AtomORdecorators)

        return move_prob * change_prob

    def change_ORdecorator(self, component, decorators):
        """
        Makes changes to the decorators associated with 1 ORbase for
        a given component (just atoms in this case)
        returns the probability of making this change
        """
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        currentORs = component.getORtypes()
        if len(currentORs) > 0:
            # Pick an ORtype (base, decorators) to make changes to
            change_OR, base_prob = _PickFromWeightedChoices((currentORs, None))

            # Remove from currentORs to make changes
            currentORs.remove(change_OR)

        # get new decorators and add back into the ORtypes list
        new_decs, new_prob = self.add_swap_delete(change_OR[1], new_decor, decor_prob*base_prob, None)

        currentORs.append( (change_OR[0], new_decs))
        component.setORtypes(currentORs)
        return new_prob

    def change_ORbase(self, component, bases, decorators):
        """
        Changes a component (atom or bond)'s ORtype list from
            a given set of bases and decorators
        Returns probability of making the move
        """
        new_base, base_prob = _PickFromWeightedChoices(bases)
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        current = component.getORtypes()
        new_OR = (new_base, [new_decor])
        new_list, new_prob = self.add_swap_delete(current, new_OR, base_prob*decor_prob, None)
        component.setORtypes(new_list)
        return new_prob

    def change_ANDdecorators(self, component, decorators):
        """
        Changs a component's (atom or bond) ANDtype list from
            a given set of decorators
        Returns probability of making the change
        """
        new_decor, decor_prob = _PickFromWeightedChoices(decorators)
        current = component.getANDtypes()
        new_list, new_prob = self.add_swap_delete(current, new_decor, decor_prob, None)
        component.setANDtypes(new_list)
        return new_prob

    def change_bond(self,env,bond):
        """
        Makes changes to the Bond object
        returns probability of making change
        """
        # Can only make changes to the bond OR or AND types
        changeOR = random.choice([True, False], p = [0.7, 0.3])
        if changeOR:
            # Bonds only have OR bases (no ORdecorators)
            new_prob = self.change_ORbase(bond, self.BondORbases, self.BondORdecorators)
            return 0.7 * new_prob

        else: # change AND type
            new_prob = self.change_ANDdecorators(bond, self.BondANDdecorators)
            return new_prob * 0.3

    def create_new_environment(self, env):
        """
        Given a parent environment type it creates a new child environment
        returns child environment type and probability of creating it
        """
        new_env = copy.deepcopy(env)
        new_env.label = _get_new_label([e.label for e in self.envList])

        if len(env.getBonds()) == 0:
            choices = ([True], None)
        else:
            choices = ([True, False], [0.9, 0.1])

        # pick to make changes to an atom or bond
        change_atom, choice_prob = _PickFromWeightedChoices(choices)

        if change_atom:
            atom, atom_prob = self.pick_an_atom(new_env)
            change_prob = self.change_atom(new_env, atom)
            prob = choice_prob * atom_prob * change_prob

        else: # Change Bond
            bond, bond_prob = self.pick_a_bond(new_env)
            change_prob = self.change_bond(new_env, bond)
            prob = choice_prob * bond_prob * change_prob

        return new_env, prob

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
        if random.random() < 0.2:
            # TODO: determine how frequently to destroy entire types
            self.log.write("Attempting to destroy type %s : %s...\n" % (env.label, env.asSMIRKS()))

            # Reject deletion of (populated) base types as we want to retain
            # generics even if empty
            if env.label in [e.label for e in self.baseTypes]:
                self.log.write("Destruction rejected for type %s because this is a generic type which was initially populated.\n" % env.label)
                return False

            # Delete the type.
            proposed_envList.remove(env)

            # Try to type all molecules.
            typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            if not self.check_typed_molecules(typelist):
                self.log.write("Typing failed; rejecting.\n")
                return False

            # save parent and children for this type
            parent = proposed_parents[env.label]['parent']
            children = proposed_parents[env.label]['children']
            # update parent dicitonary
            proposed_parents[parent]['children'] += children
            proposed_parents[parent]['children'].remove(env.label)
            for c in children:
                proposed_parents[c]['parent'] = parent
            del proposed_parents[env.label]

        else: # create new type from chosen environment
            new_env, prob = self.create_new_environment(env)

            self.log.write("Attempting to create new subtype: '%s' (%s) from parent type '%s' (%s)\n" % (new_env.label, new_env.asSMIRKS(), env.label, env.asSMIRKS()))
            self.log.write("\tProbability of making this environment is %.3f %%" % prob)

            # Check the SMIRKS for new_env is valid
            qmol = OEQMol()
            smirks = new_env.asSMIRKS()
            if self.replacements is not None:
                smirks = OESmartsLexReplace(smirks, self.replacements)
            if not OEParseSmarts(qmol, smirks):
                self.log.write("Type '%s' (%s) is invalid; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if new_env is already in types with no matches
            if new_env.asSMIRKS() in self.types_with_no_matches:
                self.log.write("Type '%s' (%s) unused in dataset; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # Check if proposed type is already in set.
            if new_env.asSMIRKS() in [e.asSMIRKS() for e in self.envList]:
                self.log.write("Type '%s' (%s) already exists; rejecting to avoid duplication.\n" % (new_env.label, new_env.asSMIRKS()))
                return False

            # add new type to proposed list
            proposed_envList.append(new_env)
            proposed_typelist = [ [e.asSMIRKS(), e.label] for e in proposed_envList]
            # Compute updated statistics
            [proposed_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_typelist)

            # Reject if new type matches nothing
            if proposed_typecounts[new_env.label] == 0:
                self.log.write("Type '%s' (%s) unused in dataset; rejecting.\n" % (new_env.label, new_env.asSMIRKS()))
                self.types_with_no_matches.append(new_env.asSMIRKS())
                return False

            # Reject if any type is emptied (UNLESS it is a basetype
            base_labels = [e.label for e in self.baseTypes]
            for label, count in proposed_typecounts.items():
                if not label in base_labels:
                    if count == 0:
                        self.log.write("Fragment with typename %s now has no matches, rejecting.\n" % label)
                        return False

            # updated proposed parent dictionary
            proposed_parents[env.label]['children'].append(new_env.label)
            proposed_parents[new_env.label] = {}
            proposed_parents[new_env.label]['children'] = list()
            proposed_parents[new_env.label]['parent'] = env.label
            proposed_parents[new_env.label]['SMIRKS'] = new_env.asSMIRKS()

        self.log.write('Proposal is valid...\n')

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
        self.log.write('Proposal score: %d >> %d : log_P_accept = %.5e\n' % (self.total_type_matches, proposed_total_type_matches, log_P_accept))
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
        molecule_typecounts (dict) - number of molecules that contain each fragment type

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

    def write_type_statistics(self, typelist, typecounts, molecule_typecounts, type_matches=None):
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
            self.log.write("%5s   %10sS %10s   %-50s %-50s %30s\n" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS', 'REF TYPE: SMIRKS', 'FRACTION OF REF TYPED MOLECULES MATCHED'))
        else:
            self.log.write("%5s   %10sS %10s   %-50s\n" % ('INDEX', self.typetag.upper(), 'MOLECULES', 'TYPE NAME: SMIRKS'))

        # Print counts
        for [smarts, typename] in typelist:
            current_combo = "%s: %s" % (typename, smarts)
            if type_matches is not None:
                (reference_typename, reference_count) = reference_type_info[typename]
                if reference_typename is not None:
                    reference_total = self.reference_type_counts[reference_typename]
                    reference_fraction = float(reference_count) / float(reference_total)
                    reference_combo = "%s: %s" % (reference_typename, self.reference_typename_dict[reference_typename])
                    self.log.write("%5d : %10d %10d | %-50s %-50s %7d / %7d (%7.3f%%)\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo, reference_combo, reference_count, reference_total, reference_fraction*100))
                else:
                    self.log.write("%5d : %10d %10d | %-50s\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo))
            else:
                self.log.write("%5d : %10d %10d | %-50s\n" % (index, typecounts[typename], molecule_typecounts[typename], current_combo))

            ntypes += typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)

        if type_matches is not None:
            frac_match = (float(self.total_type_matches) / float(self.total_types))
            self.log.write("%5s : %10d %10d |  %15s %32s %8d / %8d match (%.3f %%)\n" % ('TOTAL', ntypes, nmolecules, '', '', self.total_type_matches, self.total_types, frac_match * 100))
            return frac_match
        else:
            self.log.write("%5s : %10d %10d\n" % ('TOTAL', ntypes, nmolecules))
        return 0.0

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

    def write_parent_tree(self, roots, start='', verbose=False):
        """
        Recursively writes out the parent tree.

        Parameters
        ----------
        roots = list of typenames to print with this start
        strart = string to print at beginning of current line
        verbose = boolean, print to commandline
        """
        for r in roots:
            branch_string = "%s%s (%s)" % (start,r,self.parents[r]['SMIRKS'])
            self.log.write("%s\n" % (branch_string))
            if verbose: print(branch_string)
            new_roots = self.parents[r]['children']
            self.write_parent_tree(new_roots, start+'\t', verbose)

    def run(self, niterations, verbose = False):
        """
        Run sampler for the specified number of iterations.

        Parameters
        ----------
        niterations : int
            The specified number of iterations
        verbose : boolean, optional, default = False
            if True prints information for each iteration

        Returns
        ----------
        fraction_matched : float
            fraction of total types matched successfully at end of run

        """
        self.traj = []
        for iteration in range(niterations):
            itinfo = "Iteration %d / %d" % (iteration, niterations)
            self.log.write(itinfo+'\n')
            if verbose: print(itinfo)

            accepted = self.sample_types()
            typelist = [[env.asSMIRKS(), env.label] for env in self.envList]
            [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)

            # Get data as list of csv strings
            lines = self.save_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
            # Add lines to trajectory with iteration number:
            for l in lines:
                self.traj.append('%i,%s\n' % (iteration, l))

            if accepted:
                self.log.write('Accepted.\n')
            else:
                self.log.write("Rejected.\n")

            # Compute type statistics on molecules.
            frac_match = self.write_type_statistics(typelist, typecounts, molecule_typecounts, type_matches=self.type_matches)
            if verbose: print("Current Score: %.2f %%" % (frac_match * 100.0))
            self.log.write('\n')

            # Print parent tree as it is now.
            self.log.write("%s type hierarchy: \n" % self.typetag)
            if verbose: print("Current %s type hierarchy:" % self.typetag)
            roots = [e.label for e in self.baseTypes]
            self.write_parent_tree(roots, '\t', verbose)
            self.log.write('\n\n')
            if verbose: print('')

        # Make trajectory file
        f = open("%s.csv" % self.output, 'w')
        start = ['Iteration,Index,Smarts,ParNum,ParentParNum,RefType,Matches,Molecules,FractionMatched,Denominator\n']
        f.writelines(start + self.traj)
        f.close()

        #Compute final type stats
        typelist = [ [env.asSMIRKS(), env.label] for env in self.envList]
        [self.type_matches, self.total_type_matches] = self.best_match_reference_types(typelist)
        [typecounts, molecule_typecounts] = self.compute_type_statistics(typelist)
        fraction_matched = self.write_type_matches(typelist, self.type_matches)

        self.log.write("%s type hierarchy: \n" % self.typetag)
        roots = [e.label for e in self.baseTypes]
        self.write_parent_tree(roots, '\t', verbose)
        return fraction_matched

    def write_results_smarts_file(self):
        """
        Creates a smarts file with the most recent SMIRKS patterns compared
        to references SMIRKS patterns the format follows smarts files
        lines beginning with % are comments and the final format has

        % Results for sampling (typetag) at # temperature
        SMIRKS      label
        % matched reference SMIRKS      reference label
        ...
        % Final score was # %%
        file name is *_results.smarts
            where * is the provided output base
        returns results file name
        """
        # open results file
        smarts_file = self.output+"_results.smarts"
        smarts = open(smarts_file, 'w')

        # header
        smarts.write("%% Results for sampling %ss at %.2e\n" % (self.typetag, self.temperature))
        smarts.write("%% SMIRKS patterns for final results are below\n")

        # if there is a reference SMIRFF
        if self.SMIRFF is not None:
            smarts.write("%% followed by a their matched reference SMIRKS from %s\n" % (self.SMIRFF))
            score = 100.0 * (float(self.total_type_matches) / float(self.total_types))
            smarts.write("%%Final Score was %.3f %%\n" % (score))
            smarts.write("%%\n")
            # Make dictionary to retrieve SMIRKS
            current_dict = dict()
            for env in self.envList:
                current_dict[env.label] = env.asSMIRKS()
            # loop through current type_matches
            for (current, ref, count) in self.type_matches:
                smarts.write("%-50s %-20s\n" % (current_dict[current], current))
                if ref is not None:
                    smarts.write("%% %-48s %-20s\n" % (self.reference_typename_dict[ref], ref))

        else: # no reference SMIRFF, just print current SMIRKS
            smarts.write("%% No reference SMIRFF provided\n")
            for env in self.envList:
                smarts.write("%-50s %-20s\n" % (env.asSMIRKS(), env.label))

        smarts.close()
        return smarts_file

