"""
Command-line driver example for SMIRKY.

"""

import sys
import string
import time

from optparse import OptionParser # For parsing of command line arguments
import smarty

import os
import math
import copy
import re
import numpy
from numpy import random

def main():
    # Create command-line argument options.
    usage_string = """\
    Sample over fragment types (atoms, bonds, angles, torsions, or impropers)
    optionally attempting to match created types to an established SMIRFF.
    For all files left blank, they will be taken from this module's
    data/odds_files/ subdirectory.

    usage %prog --molecules molfile --typetag fragmentType
            [--atomORbases AtomORbaseFile --atomORdecors AtomORdecorFile
            --atomANDdecors AtomANDdecorFile --bondORbase BondORbaseFile
            --bondANDdecors BondANDdecorFile --atomIndexOdds AtomIndexFile
            --bondIndexOdds BondIndexFile --replacements substitutions
            --initialtypes initialFragmentsFile --SMIRFF referenceSMIRFF
            --temperature float --verbose verbose
            --iterations iterations --output outputFile]

    example:
    smirky --molecules AlkEthOH_test_filt1_ff.mol2 --typetag Angle

    """
    version_string = "%prog %__version__"
    parser = OptionParser(usage=usage_string, version=version_string)

    parser.add_option("-m", "--molecules", metavar='MOLECULES',
            action="store", type="string", dest='molecules_filename', default=None,
            help="Small molecule set (in any OpenEye compatible file format) containing 'dG(exp)' fields with experimental hydration free energies. This filename can also be an option in this module's data/molecules sub-directory")
    #TODO: ask about the the dG(exp) fields?

    parser.add_option("-T", "--typetag", metavar='TYPETAG',
            action = "store", type="choice", dest='typetag',
            default=None, choices = ['VdW', 'Bond', 'Angle', 'Torsion', 'Improper'],
            help="type of fragment being sampled, options are 'VdW', 'Bond', 'Angle', 'Torsion', 'Improper'")

    parser.add_option('-e', '--atomORbases', metavar="DECORATORS",
            action='store', type='string', dest='atom_OR_bases',
            default = 'odds_files/atom_OR_bases.smarts',
            help="Filename defining atom OR bases and associated probabilities. These are combined with atom OR decorators in SMIRKS, for example in '[#6X4,#7X3;R2:2]' '#6' and '#7' are atom OR bases. (OPTIONAL)")

    parser.add_option("-O", "--atomORdecors", metavar="DECORATORS",
            action='store', type='string', dest='atom_OR_decorators',
            default = 'odds_files/atom_decorators.smarts',
            help="Filename defining atom OR decorators and associated probabilities. These are combined with atom bases in SMIRKS, for example in '[#6X4,#7X3;R2:2]' 'X4' and 'X3' are ORdecorators. (OPTIONAL)")

    parser.add_option('-A', '--atomANDdecors', metavar="DECORATORS",
            action='store', type='string', dest='atom_AND_decorators',
            default='odds_files/atom_decorators.smarts',
            help="Filename defining atom AND decorators and associated probabilities. These are added to the end of an atom's SMIRKS, for example in '[#6X4,#7X3;R2:2]' 'R2' is an AND decorator. (OPTIONAL)")

    parser.add_option('-o', '--bondORbase', metavar="DECORATORS",
            action='store', type='string', dest='bond_OR_bases',
            default='odds_files/bond_OR_bases.smarts',
            help="Filename defining bond OR bases and their associated probabilities. These are OR'd together to describe a bond, for example in '[#6]-,=;@[#6]' '-' and '=' are OR bases. (OPTIONAL)")

    parser.add_option('-a', '--bondANDdecors', metavar="DECORATORS",
            action="store", type='string', dest='bond_AND_decorators',
            default='odds_files/bond_AND_decorators.smarts',
            help="Filename defining bond AND decorators and their associated probabilities. These are AND'd to the end of a bond, for example in '[#6]-,=;@[#7]' '@' is an AND decorator.(OPTIONAL)")

    parser.add_option('-D', '--atomOddsFile', metavar="ODDSFILE",
            action="store", type="string", dest="atom_odds",
            default='odds_files/atom_index_odds.smarts',
            help="Filename defining atom descriptors and probabilities with making changes to that kind of atom. Options for descriptors are integers corresponding to that indexed atom, 'Indexed', 'Unindexed', 'Alpha', 'Beta', 'All'. (OPTIONAL)")

    parser.add_option('-d', '--bondOddsFile', metavar="ODDSFILE",
            action="store", type="string", dest="bond_odds",
            default='odds_files/bond_index_odds.smarts',
            help="Filename defining bond descriptors and probabilities with making changes to that kind of bond. Options for descriptors are integers corresponding to that indexed bond, 'Indexed', 'Unindexed', 'Alpha', 'Beta', 'All'. (OPTIONAL)")

    parser.add_option("-s", "--substitutions", metavar="SUBSTITUTIONS",
            action="store", type="string", dest='substitutions_filename',
            default=None,
            help="Filename defining substitution definitions for SMARTS atom matches. (OPTIONAL).")

    parser.add_option("-f", "--initialtypes", metavar='INITIALTYPES',
            action="store", type="string", dest='initialtypes_filename',
            default=None,
            help="Filename defining initial fragment types. The file is formatted with two columns: 'SMIRKS    typename'. For the default the initial type will be a generic form of the given fragment, for example '[*:1]~[*:2]' for a bond (OPTIONAL)")

    parser.add_option('-r', '--smirff', metavar='REFERENCE',
            action='store', type='string', dest='SMIRFF',
            default=None,
            help="Filename defining a SMIRFF force fielce used to determine reference fragment types in provided set of molecules. It may be an absolute file path, a path relative to the current working directory, or a path relative to this module's data subdirectory (for built in force fields). (OPTIONAL)")

    parser.add_option("-i", "--iterations", metavar='ITERATIONS',
            action="store", type="int", dest='iterations',
            default=150,
            help="MCMC iterations.")

    parser.add_option("-t", "--temperature", metavar='TEMPERATURE',
            action="store", type="float", dest='temperature',
            default=0.1,
            help="Effective temperature for Monte Carlo acceptance, indicating fractional tolerance of mismatched atoms (default: 0.1). If 0 is specified, will behave in a greedy manner.")

    parser.add_option("-p", "--output", metavar='OUTPUT',
            action="store", type="string", dest='outputfile',
            default=None,
            help="Filename base for output information. This same base will be used for all output files created. If None provided then it is set to 'typetag_temperature' (OPTIONAL).")

    parser.add_option('-v', '--verbose', metavar='VERBOSE',
            action='store', type='choice', dest='verbose',
            default=False, choices = ['True', 'False'],
            help="If True prints minimal information to the commandline during iterations. (OPTIONAL)")

    # Parse command-line arguments.
    (option,args) = parser.parse_args()

    # Molecules are required
    if option.molecules_filename is None:
        parser.print_help()
        parser.error("Molecules input files must be specified.")

    verbose = option.verbose == 'True'
    # Load and type all molecules in the specified dataset.
    molecules = smarty.utils.read_molecules(option.molecules_filename, verbose=verbose)

    # Parse input odds files
    atom_OR_bases = smarty.utils.parse_odds_file(option.atom_OR_bases, verbose)
    atom_OR_decorators = smarty.utils.parse_odds_file(option.atom_OR_decorators, verbose)
    atom_AND_decorators = smarty.utils.parse_odds_file(option.atom_AND_decorators, verbose)
    bond_OR_bases = smarty.utils.parse_odds_file(option.bond_OR_bases, verbose)
    bond_AND_decorators = smarty.utils.parse_odds_file(option.bond_AND_decorators, verbose)
    atom_odds = smarty.utils.parse_odds_file(option.atom_odds, verbose)
    bond_odds = smarty.utils.parse_odds_file(option.bond_odds, verbose)

    # get initial types if provided, otherwise none
    if option.initialtypes_filename is None:
        initialtypes = None
    else:
        initialtypes = smarty.AtomTyper.read_typelist(option.initialtypes_filename)

    output = option.outputfile
    if output is None:
        output = "%s_%.2e" % ( option.typetag, option.temperature)
    # get replacements
    if option.substitutions_filename is None:
        sub_file = smarty.get_data_filename('odds_files/substitutions.smarts')
    else:
        sub_file = option.substitutions_filename
    replacements = smarty.AtomTyper.read_typelist(sub_file)
    replacements = [ (short, smarts) for (smarts, short) in replacements]

    start_sampler = time.time()
    fragment_sampler = smarty.FragmentSampler(
            molecules, option.typetag, atom_OR_bases, atom_OR_decorators,
            atom_AND_decorators, bond_OR_bases, bond_AND_decorators,
            atom_odds, bond_odds, replacements, initialtypes,
            option.SMIRFF, option.temperature, output)
    # report time
    finish_sampler = time.time()
    elapsed = finish_sampler - start_sampler
    if verbose: print("Creating %s sampler took %.3f s" % (option.typetag, elapsed))

    # Make iterations
    frac_found = fragment_sampler.run(option.iterations, verbose)
    finished = time.time()
    elapsed = finished - finish_sampler
    per_it = elapsed / float(option.iterations)
    if verbose: print("%i iterations took %.3f s (%.3f s / iteration)" % (option.iterations, elapsed, per_it))
    if verbose: print("Final score was %.3f %%" % (frac_found*100.0))

    # plot results
    plot_file = "%s.pdf" % output
    traj = "%s.csv" % output
    smarty.score_utils.create_plot_file(traj, plot_file, False)
