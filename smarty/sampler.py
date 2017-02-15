#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
smarty.py
Example illustrating a scheme to create and destroy atom types automatically using SMARTS.
AUTHORS
John Chodera <john.chodera@choderalab.org>, Memorial Sloan Kettering Cancer Center.
Additional contributions from the Mobley lab, UC Irvine, including David Mobley, Caitlin Bannan, and Camila Zanette.
The AtomTyper class is based on 'patty' by Pat Walters, Vertex Pharmaceuticals.
"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
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

import networkx as nx

import time

from smarty.atomtyper import AtomTyper
from smarty.score_utils import load_trajectory
from smarty.score_utils import scores_vs_time

#=============================================================================================
# ATOMTYPE SAMPLER
#=============================================================================================

class AtomTypeSampler(object):
    """
    Atom type sampler.
    """
    def __init__(self, molecules, basetypes_filename, initialtypes_filename, decorators_filename, replacements_filename=None, reference_typed_molecules=None, temperature=0.1, verbose=False, decorator_behavior='combinatorial-decorators', element = 0):
        """
        Initialize an atom type sampler.
        ARGUMENTS
        molecules : list of molecules for typing
            List of molecules for typing
        basetypes_filename : str
            File defining base/generic atom types (which cannot be destroyed); often these are elemental types
        initialtypes_filename :
            File defining initial atom types (which CAN be destroyed, except for those which occur in basetypes_filename
        decorators_filename : str
            File containing decorators that can be added to existing types to generate subtypes
        replacements_filename : str, optional, default=None
            If specified, SMARTS replacement definitions will be read from this file
        reference_typed_molecules : list of OEMol, optional, default=None
            List of molecules with reference types for use in Monte Carlo acceptance.
            If specified, the likelihood function will utilize the maximal number of matched atom types with these molecules.
            If not specified, no likelihood function will be employed.
        temperature : float, optional, default=0.1
            Temperature for Monte Carlo acceptance/rejection
        verbose : bool, optional, default=False
            If True, verbose output will be printed.
        decorator_behavior : string either "combinatorial-decorators" or "simple-decorators"
            simple decorators include bonded atoms as decorators
        element : iteger >= 0
            If 0 all atomtypes sampled, otherwise only atomtypes of that atomic number are sampled
        Notes
        -----
        This is just a proof of concept for chemical perception sampling.
        Scoring for purposed atomtypes is based on reference atomtypes.
        No scoring of molecular properties is performed
        """
        # store simple input information
        self.verbose = verbose
        self.decorator_behavior = decorator_behavior
        self.typetag = 'atomtype' # internal tag
        self.temperature = temperature
        self.element = element

        # Read atomtypes (initial and base) and decorators.
        self.atomtypes = AtomTyper.read_typelist(initialtypes_filename)
        self.basetypes = AtomTyper.read_typelist(basetypes_filename)
        self.decorators = AtomTyper.read_typelist(decorators_filename)
        self.replacements = AtomTyper.read_typelist(replacements_filename)

        # Store a deep copy of the molecules since they will be annotated
        self.molecules = copy.deepcopy(molecules)

        # Save bond list to use throughout
        bondset = [("-","singly"), ("=", "doubly"), ("#","triply"), (":", "aromatic")]

        # Calculate which bonds in set are used
        bond_typelist = [("[*]%s[*]" %bond, name) for (bond, name) in bondset]
        tmpmolecules = copy.deepcopy(molecules)
        self.type_molecules(bond_typelist, tmpmolecules, 0)
        [bond_typecounts, molecule_bond_typecounts] = self.compute_type_statistics( bondset, tmpmolecules, 0)
        if self.verbose:
            print("USED BOND TYPES:")
            self.show_type_statistics(bondset, bond_typecounts, molecule_bond_typecounts)

        # only same bonds that are used
        self.bondset = [ ('~', 'any') ]
        for (bond, name) in bondset:
            if bond_typecounts[name] > 0:
                self.bondset.append( (bond, name) )

        # Rename base/initial types to ensure their names are unique
        # clashes between initial and target types will cause problems
        for idx, [smarts, typename] in enumerate(self.atomtypes):
            self.atomtypes[idx] = (smarts, 'c_'+typename)
        for idx, [smarts, typename] in enumerate(self.basetypes):
            self.basetypes[idx] = (smarts, 'c_'+typename)

        # Store smarts for basetypes
        self.basetypes_smarts = [ smarts for (smarts, name) in self.basetypes ]

        # Add any base types not already there to the initial types
        initial_smarts = [ smarts for (smarts, name) in self.atomtypes ]
        for [smarts, typename] in self.basetypes:
            if smarts not in initial_smarts:
                self.atomtypes = [(smarts, typename)] + self.atomtypes
                if self.verbose: print("Added base (generic) type `%s`, name %s, to initial types." % (smarts, typename) )

        # Maintain a list of SMARTS matches without any atom type matches in the dataset
        # This is used for efficiency.
        self.atomtypes_with_no_matches = set()

        tmpmolecules = copy.deepcopy(molecules)
        self.type_molecules(self.basetypes, tmpmolecules)
        [ basetype_typecounts, molecule_basetype_typecounts] = self.compute_type_statistics( self.basetypes, tmpmolecules, self.element)
        if self.verbose:
            print("MATCHED BASETYPES:")
            self.show_type_statistics(self.basetypes, basetype_typecounts, molecule_basetype_typecounts)

        # store used basedtypes for reference later
        used_basetypes = list()

        # Track only used base types and add unused to atomtypes with no matches
        for (smarts, atom_type) in self.basetypes:
            # If this type is used, then track it
            if basetype_typecounts[atom_type] > 0:
                used_basetypes.append( ( smarts, atom_type) )
            # If unused, it will be removed and the smarts stored in the no match list
            else:
                self.atomtypes_with_no_matches.add( smarts )
                if self.verbose: print("Removing basetype '%s' ('%s'), which is unused." % (smarts, atom_type))
        # Atom basetypes to create new smart strings
        self.basetypes = copy.deepcopy(used_basetypes)

        # Type all molecules with current typelist to ensure that starting types are sufficient.
        self.type_molecules(self.atomtypes, self.molecules, self.element)
        # Compute atomtype statistics on molecules for current atomtype set
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules, self.element)
        if self.verbose:
            print("MATCHED INITIAL TYPES:")
            self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts)

        # Track only used atomtypes and add unused to atomtypes with no matches
        used_initial_atomtypes = list()
        for (smarts, atom_type) in self.atomtypes:
            if atom_typecounts[atom_type] > 0:
                used_initial_atomtypes.append( (smarts, atom_type) )
            else:
                self.atomtypes_with_no_matches.add( smarts )
                if self.verbose: print("Removing initial atom type '%s', as it matches no atoms" % smarts)
        self.atomtypes = copy.deepcopy(used_initial_atomtypes)
        self.initial_atomtypes = copy.deepcopy(used_initial_atomtypes)

        # Type molecules again with the updated atomtype list
        self.type_molecules(self.atomtypes, self.molecules, self.element)

        # These are atomtypes where not all children have been matched
        self.unmatched_atomtypes = copy.deepcopy(self.atomtypes)

        # Creat dictionary to store children of initial atom types
        self.parents = dict()
        for [smarts, typename] in self.atomtypes:
            #store empty list of chlidren for each atomtype
            self.parents[smarts] = []
        # store reverse parent dictionary with child to parent
        self.child_to_parent = self._switch_parent_dict()

        # Compute total atoms
        self.total_atoms = 0.0
        for molecule in self.molecules:
            for atom in self._GetAtoms(molecule, self.element):
                self.total_atoms += 1.0

        # Store reference molecules
        self.reference_typed_molecules = None
        self.reference_atomtypes = set()
        self.current_atom_matches = None
        if reference_typed_molecules is not None:
            self.reference_typed_molecules = copy.deepcopy(reference_typed_molecules)
            # Extract list of reference atom types
            for molecule in reference_typed_molecules:
                for atom in self._GetAtoms(molecule, self.element):
                    atomtype = atom.GetType()
                    self.reference_atomtypes.add(atomtype)
            self.reference_atomtypes = list(self.reference_atomtypes)
            # Compute current atom matches
            [self.atom_type_matches, self.total_atom_type_matches] = self.best_match_reference_types(self.atomtypes, self.molecules)
            # Count atom types.
            self.reference_atomtypes_atomcount = { atomtype : 0 for atomtype in self.reference_atomtypes }
            for molecule in reference_typed_molecules:
                for atom in self._GetAtoms(molecule, self.element):
                    atomtype = atom.GetType()
                    self.reference_atomtypes_atomcount[atomtype] += 1
        return

    def _GetAtoms(self, molecule, element = 0):
        """
        Parameters
        ----------
        molecule : OEMol
        element : integer
            if 0 looks at all atoms, otherwise only those with the given atomic number

        Returns
        -------
        iterator over the atoms based on the molecule and element number
        """
        if element > 0:
            return molecule.GetAtoms(OEHasAtomicNum(element))
        else:
            return molecule.GetAtoms()

    def best_match_reference_types(self, atomtypes, molecules):
        """
        Determine best match for each parameter with reference atom types
        Parameters
        ----------
        atomtypes :
            Current atom types
        molecules : list of OEMol
            Typed molecules, where types are stored in self.atomtypetag string data.
        Returns
        -------
        atom_type_matches : list of tuples (current_atomtype, reference_atomtype, counts)
            Best correspondence between current and reference atomtypes, along with number of atoms equivalently typed in reference molecule set.
        total_atom_type_matches : int
            The total number of correspondingly typed atoms in the reference molecule set.
        * Currently, types for reference typed molecules are accessed via atom.GetType(), while types for current typed molecules are accessed via atom.GetStringData(self.typetag).
          This should be homogenized.
        Contributor:
        * Josh Fass <josh.fass@choderalab.org> contributed this algorithm.
        """
        if self.reference_typed_molecules is None:
            if self.verbose: print('No reference molecules specified, so skipping likelihood calculation.')
            return None

        # Create bipartite graph (U,V,E) matching current atom types U with reference atom types V via edges E with weights equal to number of atoms typed in common.
        if self.verbose: print('Creating graph matching current atom types with reference atom types...')
        initial_time = time.time()
        graph = nx.Graph()

        # Get current atomtypes and reference atom types
        current_atomtypes = [ typename for (smarts, typename) in atomtypes ]
        reference_atomtypes = [ typename for typename in self.reference_atomtypes ]
        # check that current atom types are not in reference atom types
        if set(current_atomtypes) & set(reference_atomtypes):
            raise Exception("Current and reference atom types must be unique")
        # Add current atom types
        for atomtype in current_atomtypes:
            graph.add_node(atomtype, bipartite=0)
        # Add reference atom types
        for atomtype in reference_atomtypes:
            graph.add_node(atomtype, bipartite=1)
        # Add edges.
        atoms_in_common = dict()
        for current_atomtype in current_atomtypes:
            for reference_atomtype in reference_atomtypes:
                atoms_in_common[(current_atomtype,reference_atomtype)] = 0
        for (current_typed_molecule, reference_typed_molecule) in zip(molecules, self.reference_typed_molecules):
            current_atoms = self._GetAtoms(current_typed_molecule, self.element)
            reference_atoms = self._GetAtoms(reference_typed_molecule, self.element)
            for (current_typed_atom, reference_typed_atom) in zip(current_atoms, reference_atoms):
                current_atomtype = current_typed_atom.GetStringData(self.typetag)
                reference_atomtype = reference_typed_atom.GetType()
                atoms_in_common[(current_atomtype,reference_atomtype)] += 1
        for current_atomtype in current_atomtypes:
            for reference_atomtype in reference_atomtypes:
                weight = atoms_in_common[(current_atomtype,reference_atomtype)]
                graph.add_edge(current_atomtype, reference_atomtype, weight=weight)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Graph creation took %.3f s' % elapsed_time)

        # Compute maximum match
        if self.verbose: print('Computing maximum weight match...')
        initial_time = time.time()
        mate = nx.algorithms.max_weight_matching(graph, maxcardinality=False)
        elapsed_time = time.time() - initial_time
        if self.verbose: print('Maximum weight match took %.3f s' % elapsed_time)

        # Compute match dictionary and total number of matches.
        atom_type_matches = list()
        total_atom_type_matches = 0
        for current_atomtype in current_atomtypes:
            if current_atomtype in mate:
                reference_atomtype = mate[current_atomtype]
                counts = graph[current_atomtype][reference_atomtype]['weight']
                total_atom_type_matches += counts
                atom_type_matches.append( (current_atomtype, reference_atomtype, counts) )
            else:
                atom_type_matches.append( (current_atomtype, None, None) )

        # Report on matches
        if self.verbose:
            print("PROPOSED:")
            self.show_type_matches(atom_type_matches)

        return (atom_type_matches, total_atom_type_matches)

    def show_type_matches(self, atom_type_matches):
        """
        Show pairing of current to reference atom types.
        Parameters
        ----------
        atom_type_matches : list of (current_atomtype, reference_atomtype, counts)

        Returns
        -------
        fraction_matched_atoms : the fractional count of matched atoms
        """
        print('Atom type matches:')
        total_atom_type_matches = 0
        for (current_atomtype, reference_atomtype, counts) in atom_type_matches:
            if reference_atomtype is not None:
                print('%-64s matches %8s : %8d atoms matched' % (current_atomtype, reference_atomtype, counts))
                total_atom_type_matches += counts
            else:
                print('%-64s         no match' % (current_atomtype))

        fraction_matched_atoms = float(total_atom_type_matches) / float(self.total_atoms)
        print('%d / %d total atoms match (%.3f %%)' % (total_atom_type_matches, self.total_atoms, fraction_matched_atoms * 100))

        return fraction_matched_atoms


    def AtomDecorator(self, atom1type, decorator):
        """
        Given an atom and a decorator ammend the SMARTS string with that decorator

        Parameters
        -----------
        atom1type : atomtype tuple in form (smarts, typename)
        decorator : decorator being added to current atom

        Returns
        -------
        decorated atomtype as a tuple (smarts, typename)
        """
        if self.HasAlpha(atom1type):
            # decorators should go before the $ sign on the atom
            dollar = atom1type[0].find('$')
            proposed_atomtype = atom1type[0][:dollar] + decorator[0] + atom3[0][dollar:]
            proposed_typename = atom1type[1] + ' ' + decorator[1]
        else:
            # No alpha atom so the decorator goes before the ']'
            proposed_atomtype = atom1type[0][:-1] + decorator[0] + ']'
            proposed_typename = atom1type[1] + ' '  + decorator[1]

        return (proposed_atomtype, proposed_typename)

    def PickAnAtom(self, atomList):
        """
        Parameters
        ----------
        atomList : any list of tuples in the form (smarts, typename)
                   this could include decorator or bond lists

        Returns
        -------
        one random (smarts, typename) pair from given list

        This allows for continuity in the code, this method could be changed,
        and all random choices would still be made in the same way.
        It also allowed for testing which atomtypes to choose from while sampling.
        """
        return random.choice(atomList)

    def HasAlpha(self, atom1type):
        """
        Parameter
        ---------
        atom1type : an atomtype tuple (smarts, typename)

        Returns
        -------
        True if atomtype has at least 1 alpha substituent otherwise False
        """
        # TODO: check does this work if you're using replacements
        # CCB: I don't think it will work!
        if atom1type[0].find("$") != -1:
            return True
        else:
            return False

    def AddAlphaSubstituentAtom(self, atom1type, bondset, atom2type):
        """
        Adds an atom alpha to the primary atom. The new alpha substituent
        always adds to the end of the sequence of alpha atom
        so if you have '[#1$(*~[#6])]' the next alpha atom [#8] is added in
        this way '[#1$(*~[#6])$(*~[#8])]'

        Parameters
        ----------
        atom1type : current atomtype (smarts, typename)
        bondset : bondtype to connect two atoms (smarts, bondname)
        atom2type : atom to be added (smarts, typename)

        Returns
        -------
        Atomtype with new alpha substituent (smarts, typename)
        """
        proposed_atomtype = atom1type[0][:len(atom1type[0])-1] + '$(*' + bondset[0] + atom2type[0] + ')]'
        proposed_typename = atom1type[1] + ' ' + bondset[1] + ' ' + atom2type[1] + ' '
        return (proposed_atomtype, proposed_typename)

    def AddBetaSubstituentAtom(self, atom1type, bondset, atom2type):
        """
        Adds atom2type as a beta substituent bonding it to the
        first alpha atom in atom1type. If atom1type does not have
        an alpha atom this metho will call addAlphaSubstituentAtom instead.

        Parameters
        ----------
        atom1type : parent atomtype (smarts, typename)
        bondset : bond used to connect atoms (smarts, bondname)
        atom2type : atomtype being bonded in beta position (smarts, typename)

        Returns
        -------
        child atomtype as tuple (smarts, typename)

        """

        # counting '[' tells us how many atoms are in the mix
        count = atom1type[0].count('[')
        proposed_atomtype = ""
        number_brackets = 0
        # find closed alpha atom
        closeAlpha = atom1type[0].find(']')
        # This has two atoms (already has an alpha atom)
        if count == 2:
            proposed_atomtype = atom1type[0][:closeAlpha+1]
            proposed_atomtype += bondset[0] + atom2type[0] + ')]'
            proposed_typename = atom1type[1] + bondset[1] + ' ' + atom2type[1]
            if self.verbose: print("ADD FIRST BETA SUB: proposed --- %s %s" % ( str(proposed_atomtype), str(proposed_typename)))
        elif count > 2:
            # Has an alpha atom with at least 1 beta atom
            proposed_atomtype = atom1type[0][:closeAlpha+1]
            proposed_atomtype += '(' + bondset[0] + atom2type[0] + ')'
            proposed_atomtype += atom1type[0][closeAlpha+1:]
            proposed_typename = atom1type[1] + ' (' + bondset[1] + ' ' + atom2type[1] + ')'
            if self.verbose: print("ADD MORE BETA SUB: proposed --- %s %s" % ( str(proposed_atomtype), str(proposed_typename)))
        else:
            # Has only 1 atom which means there isn't an alpha atom yet, add an alpha atom instead
            proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, bondset, atom2type)
        return (proposed_atomtype, proposed_typename)


    def sample_atomtypes(self):
        """
        Perform one step of atom type sampling.
        This is done by either removing a current atomtype
        or creating a child atom type. Then the proposed
        atomtype list is scored and the move is accepted or rejected
        """
        # Copy current atomtypes for proposal.
        proposed_atomtypes = copy.deepcopy(self.atomtypes)
        proposed_molecules = copy.deepcopy(self.molecules)
        proposed_parents = copy.deepcopy(self.parents)

        if random.random() < 0.5:
            # Pick a random index and remove atomtype at that index
            remove_index = random.randit(0, len(proposed_atomtypes))
            (atomtype, typename) = proposed_atomtypes.pop(remove_index)
            if self.verbose: print("Attempting to destroy atom type %s : %s..." % (atomtype, typename))

            # Reject deletion of (populated) base types as we want to retain
            # generics even if empty
            if atomtype in self.basetypes_smarts:
                if self.verbose: print("Destruction rejected for atom type %s because this is a generic type which was initially populated." % atomtype )
                return False

            # update proposed parent dictionary
            for parent, children in proposed_parents.items():
                if atomtype in [at for (at, tn) in children]:
                    children += proposed_parents[atomtype]
                    children.remove( (atomtype, typename) )

            del proposed_parents[atomtype]

            # Try to type all molecules.
            try:
                self.type_molecules(proposed_atomtypes, proposed_molecules, self.element)
            except AtomTyper.TypingException as e:
                # Reject since typing failed.
                if self.verbose: print("Typing failed; rejecting.")
                return False
        else:
            if self.decorator_behavior == 'simple-decorators':
                # Pick an atomtype to subtype.
                atom1type = self.PickAnAtom(self.atomtypes)
                # Pick a decorator to add.
                (decorator, decorator_typename) = self.PickAnAtom(self.decorators)

                # Create new atomtype to insert by appending decorator with 'and' operator.
                result = re.match('\[(.+)\]', atom1type[0])
                proposed_atomtype = '[' + result.groups(1)[0] + '&' + decorator + ']'
                proposed_typename = atom1type[1] + ' ' + decorator_typename
                if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], decorator, decorator_typename, proposed_atomtype, proposed_typename))

            else: # combinatorial-decorators
                # Pick an atomtype
                atom1type = self.PickAnAtom(self.atomtypes)
                # Check if we need to add an alpha or beta substituent
                if self.HasAlpha(atom1type):
                    # Has alpha
                    bondtype = self.PickAnAtom(self.bondset)
                    atom2type = self.PickAnAtom(self.basetypes)
                    if random.random() < 0.5 or atom1type[0][2] == '1': # Add Beta Substituent Atom randomly or when it is Hydrogen
                        proposed_atomtype, proposed_typename = self.AddBetaSubstituentAtom(atom1type, bondtype, atom2type)
                    else: # Add another Alpha Substituent if it is not a Hydrogen
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, bondtype, atom2type)
                    if self.verbose: print("Attempting to create new subtype: -> '%s' (%s)" % (proposed_atomtype, proposed_typename))
                else:
                    # Has no alpha
                    if random.random() < 0.5: # add decorator to primary atom
                        decorator = self.PickAnAtom(self.decorators)
                        proposed_atomtype, proposed_typename = self.AtomDecorator(atom1type, decorator)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) + '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], decorator[0], decorator[1], proposed_atomtype, proposed_typename))
                    else: # add Alpha substituent
                        bondtype = self.PickAnAtom(self.bondset)
                        atom2type = self.PickAnAtom(self.basetypes)
                        proposed_atomtype, proposed_typename = self.AddAlphaSubstituentAtom(atom1type, bondtype, atom2type)
                        if self.verbose: print("Attempting to create new subtype: '%s' (%s) -> '%s' (%s)" % (atom1type[0], atom1type[1], proposed_atomtype, proposed_typename))


            # for either decorator - update proposed parent dictionary
            proposed_parents[atom1type[0]].append( (proposed_atomtype, proposed_typename) )
            proposed_parents[proposed_atomtype] = []

            # Check that we haven't already determined this atom type isn't matched in the dataset.
            if proposed_atomtype in self.atomtypes_with_no_matches:
                if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                return False

            # Check that it is a new SMARTS pattern
            if proposed_atomtype in [smarts for (smarts, typename) in self.atomtypes]:
                if self.verbose: print("Atom type '%s' (%s) is in the existing atomtype list; rejecting." % (proposed_atomtype, proposed_typename))
                return False

            # Insert atomtype immediately after.
            proposed_atomtypes.append( (proposed_atomtype, proposed_typename) )
            # Try to type all molecules.
            try:
                # Type molecules.
                self.type_molecules(proposed_atomtypes, proposed_molecules, self.element)
                # Compute updated statistics.
                [proposed_atom_typecounts, proposed_molecule_typecounts] = self.compute_type_statistics(proposed_atomtypes, proposed_molecules, self.element)
            except AtomTyper.TypingException as e:
                print("Exception: %s" % str(e))
                # Reject since typing failed.
                if self.verbose: print("Typing failed for one or more molecules using proposed atomtypes; rejecting.")
                return False

            # Reject if new type is unused.
            if (proposed_atom_typecounts[proposed_typename] == 0):
                # Reject because new type is unused in dataset.
                if self.verbose: print("Atom type '%s' (%s) unused in dataset; rejecting." % (proposed_atomtype, proposed_typename))
                # Store this atomtype to speed up future rejections
                self.atomtypes_with_no_matches.add(proposed_atomtype)
                return False

            # Reject if any type is emptied (UNLESS it is a basetype)
            for (smarts, typename) in proposed_atomtypes:
                if not smarts in self.basetype_smarts: # not a base type
                    if proposed_atom_typecounts[typename] == 0: # no matches
                        if self.verbose: print("Atomtype '%s' (%s) is now unused in dataset; rejecting." % (smarts, typename))
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
            self.child_to_parent = self._switch_parent_dict()
            self.atom_type_matches = proposed_atom_type_matches
            self.total_atom_type_matches = proposed_total_atom_type_matches
            return True
        else:
            return False

    def type_molecules(self, typelist, molecules, element = 0):
        """
        Type all molecules with the specified typelist.
        Parameters
        ----------
        typelist : list of atomtypes or tuples in the form (smarts, typename)
        molecules : list of OEMols
        element : integer 0 for all atoms or atomic number being sampled

        For every atom in each molecule the relevant typename is assigned
        so it can be accessed at atom.GetStringData(self.typetag)
        """
        # Create an atom typer.
        atomtyper = AtomTyper(typelist, self.typetag, replacements=self.replacements)

        # Type molecules.
        for molecule in molecules:
            atomtyper.assignTypes(molecule, element)

        return

    def compute_type_statistics(self, typelist, molecules, element = 0):
        """
        Compute statistics for numnber of molecules assigned each type.
        Parameters
        ----------
        typelist : atomtype list of form (smarts, typename)
        molecules : list of OEmols
        element : 0 for all atoms or atomic number being sampled
        Returns
        -------
        atom_typecounts (dict) : counts of number of atoms containing each atomtype
        molecule_typecounds (dict) : counts of number of molecules containing each atom type
        """
        # Zero type counts by atom and molecule.
        atom_typecounts = dict()
        molecule_typecounts = dict()
        for [smarts, typename] in typelist:
            atom_typecounts[typename] = 0
            molecule_typecounts[typename] = 0

        # Count number of atoms with each type.
        for molecule in molecules:
            types_in_this_molecule = set()
            for atom in self._GetAtoms(molecule, element):
                atomtype = atom.GetStringData(self.typetag)
                types_in_this_molecule.add(atomtype)
                atom_typecounts[atomtype] += 1
            for atomtype in types_in_this_molecule:
                molecule_typecounts[atomtype] += 1

        return (atom_typecounts, molecule_typecounts)

    def show_type_statistics(self, typelist, atom_typecounts, molecule_typecounts, atomtype_matches=None):
        """
        Print atom type statistics to the commandline
        Parameters
        ----------
        typelist : atomtype list of form (smarts, typename)
        atom_typecounts : dictionary result from compute_type_statistics
        molecule_typecounts : dictionary result from compute_type_statistics
        atomtype_matches : dictionary result from best_match_references_types
                           if there are reference molecules
        """
        index = 1
        natoms = 0

        if atomtype_matches is not None:
            reference_type_info = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                reference_type_info[typename] = (reference_atomtype, count)

        # Print header
        if atomtype_matches is not None:
            print("%5s   %10s %10s   %64s %32s %8s %46s" % ('INDEX', 'ATOMS', 'MOLECULES', 'TYPE NAME', 'SMARTS', 'REF TYPE', 'FRACTION OF REF TYPED MOLECULES MATCHED'))
        else:
            print("%5s   %10s %10s   %64s %32s" % ('INDEX', 'ATOMS', 'MOLECULES', 'TYPE NAME', 'SMARTS'))

        # Print counts
        for [smarts, typename] in typelist:
            if atomtype_matches is not None:
                (reference_atomtype, reference_count) = reference_type_info[typename]
                if reference_atomtype is not None:
                    reference_total = self.reference_atomtypes_atomcount[reference_atomtype]
                    reference_fraction = float(reference_count) / float(reference_total)
                    print("%5d : %10d %10d | %64s %32s %8s %16d / %16d (%7.3f%%)" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts, reference_atomtype, reference_count, reference_total, reference_fraction*100))
                else:
                    print("%5d : %10d %10d | %64s %32s" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts))
            else:
                print("%5d : %10d %10d | %64s %32s" % (index, atom_typecounts[typename], molecule_typecounts[typename], typename, smarts))

            natoms += atom_typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)

        if atomtype_matches is not None:
            print("%5s : %10d %10d |  %64s %32s %8d / %8d match (%.3f %%)" % ('TOTAL', natoms, nmolecules, '', '', self.total_atom_type_matches, self.total_atoms, (float(self.total_atom_type_matches) / float(self.total_atoms)) * 100))
        else:
            print("%5s : %10d %10d" % ('TOTAL', natoms, nmolecules))

        return

    def save_type_statistics(self, typelist, atom_typecounts, molecule_typecounts, atomtype_matches=None):
        """
        Saves the match information in format for a trajectory file
        Parameters
        ----------
        typelist : atomtype list of form (smarts, typename)
        atom_typecounts : dictionary result from compute_type_statistics
        molecule_typecounts : dictionary result from compute_type_statistics
        atomtype_matches : dictionary result from best_match_references_types
                           if there are reference molecules
        Returns
        -------
        output : string line for trajectory file
        """
        if atomtype_matches is not None:
            reference_type_info = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                reference_type_info[typename] = (reference_atomtype, count)

        index = 1
        output = []
        ntypes = 0
        # Print counts
        # INDEX, SMARTS, PARENT INDEX, REF TYPE, MATCHES, MOLECULES, FRACTION, OUT of, PERCENTAGE
        for [smarts, typename] in typelist:
            parent = str(self.child_to_parent[smarts])
            if atomtype_matches is not None:
                (reference_atomtype, reference_count) = reference_type_info[typename]
                if reference_atomtype is not None:
                    reference_total = self.reference_atomtypes_atomcount[reference_atomtype]
                    reference_fraction = float(reference_count) / float(reference_total)
                    # Save output
                    output.append("%i,'%s','%s','%s','%s',%i,%i,%i,%i" % (index, smarts, typename, parent, reference_atomtype, atom_typecounts[typename], molecule_typecounts[typename], reference_count, reference_total))
                else:
                    output.append("%i,'%s','%s','%s','%s',%i,%i,%i,%i" % (index, smarts, typename, parent, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))

            else:
                output.append("%i,'%s',%i,%i,'%s',%i,%i,%i,%i" % (index, smarts, typename, parent, 'NONE', atom_typecounts[typename], molecule_typecounts[typename], 0, 0))

            ntypes += atom_typecounts[typename]
            index += 1

        nmolecules = len(self.molecules)
        if atomtype_matches is None:
            output.append("-1,'total','all','None','all',%i,%i,0,0" % (ntypes, nmolecules))
        else:
            output.append("-1,'total','all','None','all',%i,%i,%i,%i" % (ntypes,nmolecules,self.total_atom_type_matches,self.total_atoms))
        return output

    # TODO: CCB - I think this method should be removed, I wrote it initially.
    # Our goal was to sample only types where all children hadn't been found
    # I now think this doesn't fit with the idea of sampling where we should
    # be able to find and remove good things, we shouldn't stop looking at a
    # branch of the type tree just because each child has a 100% match
    # currently the resulting list is never used for sampling
    def get_unfinishedAtomList(self, atom_typecounts, molecule_typecounts, atomtype_matches = None):
        """
        This method prunes the set of current atomtypes so that if all branches
        of a base type have been found it no longer tries extending any atom of that base type.
        """
        # Reset unmatched atom types incase something was destroyed
        self.unmatched_atomtypes = copy.deepcopy(self.atomtypes)

        # If we don't have reference matches, unmatched_atomtypes should be all current atomtypes
        if atomtype_matches is None:
            return
        else: # store counts for each atom type
            reference_counts = dict()
            for (typename, reference_atomtype, count) in atomtype_matches:
                if reference_atomtype is None:
                    reference_counts[typename] = 0
                else:
                    reference_counts[typename] = count

        # If all of a basetype and it's children match found atoms and reference remove from list
        for [base_smarts, base_typename] in self.initial_atomtypes:
            includeBase = True

            # If the number of atoms matches the references are the same for basetypes and their children
            # then we have found all reference types for that element and should stop searching that branch
            if atom_typecounts[base_typename] == reference_counts[base_typename]:
                includeBase = False
                for [child_smarts, child_name] in self.parents[base_smarts]:
                    # If any of the children's atom count and reference count don't agree then these should stay in the unmatched_atomtypes
                    if not atom_typecounts[child_name] == reference_counts[child_name]:
                        includeBase = True
                        break

            # Remove atomtypes from completed element branches
            if not includeBase:
                self.unmatched_atomtypes.remove( (base_smarts, base_typename) )
                for child in self.parents[base_smarts]:
                    self.unmatched_atomtypes.remove(child)

        return

    def _switch_parent_dict(self):
        """
        Takes the parent dictionary and returns a dictionary in the form
        {child: parent}
        """
        child_to_parent = dict()
        for smarts in self.parents.keys():
            child_to_parent[smarts] = None

        for smarts, children in self.parents.items():
            for [child_smarts, child_typename] in children:
                child_to_parent[child_smarts] = smarts

        return child_to_parent

    def print_parent_tree(self, roots, start=''):
        """
        Recursively prints the parent tree.
        Parameters
        ----------
        roots = list of smarts strings to print
        """
        for r in roots:
            print("%s%s" % (start, r))
            if r in self.parents:
                new_roots = [smart for [smart, name] in self.parents[r]]
                self.print_parent_tree(new_roots, start+'\t')


    def run(self, niterations, trajFile=None):
        """
        Run atomtype sampler for the specified number of iterations.
        Parameters
        ----------
        niterations : int
            The specified number of iterations
        trajFile : str, optional, default=None
            Output trajectory filename
        Returns
        ----------
        fraction_matched_atoms : float
            fraction of total atoms matched successfully at end of run
        """
        if trajFile is not None:
            # make "trajectory" file
            if os.path.isfile(trajFile):
                print("trajectory file already exists, it was overwritten")
            self.traj = open(trajFile, 'w')
            self.traj.write('Iteration,Index,Smarts,Typename,ParentSMARTS,RefType,Matches,Molecules,FractionMatched,Denominator\n')

        for iteration in range(niterations):
            if self.verbose:
                print("Iteration %d / %d" % (iteration, niterations))

            accepted = self.sample_atomtypes()
            [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules, self.element)
            self.get_unfinishedAtomList(atom_typecounts, molecule_typecounts, atomtype_matches = self.atom_type_matches)

            if trajFile is not None:
                # Get data as list of csv strings
                lines = self.save_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)
                # Add lines to trajectory with iteration number:
                for l in lines:
                    self.traj.write('%i,%s \n' % (iteration, l))

            if self.verbose:
                if accepted:
                    print('Accepted.')
                else:
                    print('Rejected.')

                # Compute atomtype statistics on molecules.
                self.show_type_statistics(self.atomtypes, atom_typecounts, molecule_typecounts, atomtype_matches=self.atom_type_matches)
                print('')

                # Print parent tree as it is now.
                roots = [r for r in self.child_to_parent.keys() if self.child_to_parent[r] is None]

                print("Atom type hierarchy:")
                self.print_parent_tree(roots, '\t')

        if trajFile is not None:
            self.traj.close()
            # Get/print some stats on trajectory
            # Load timeseries
            timeseries = load_trajectory( trajFile )
            time_fractions = scores_vs_time( timeseries )
            print("Maximum score achieved: %.2f" % max(time_fractions['all']))


        #Compute final type stats
        [atom_typecounts, molecule_typecounts] = self.compute_type_statistics(self.atomtypes, self.molecules, self.element)
        fraction_matched_atoms = self.show_type_matches(self.atom_type_matches)

        # If verbose print parent tree:
        if self.verbose:
            roots = self.parents.keys()
            # Remove keys from roots if they are children
            for parent, children in self.parents.items():
                child_smarts = [smarts for [smarts, name] in children]
                for child in child_smarts:
                    if child in roots:
                        roots.remove(child)

            print("Atom type hierarchy:")
            self.print_parent_tree(roots, '\t')
        return fraction_matched_atoms
