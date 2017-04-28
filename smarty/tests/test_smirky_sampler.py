import unittest
import smarty
from openforcefield.typing.chemistry.environment import *
from openforcefield.utils.utils import read_molecules
from smarty.sampler_smirky import *
from smarty import utils
from smarty import score_utils
from operator import itemgetter, attrgetter
import openeye.oechem
from openeye.oechem import *
import copy
import sys # used to exit while testing

class TestSmirkySampler(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        """
        Initialize TestCase and then read in odds from files in smarty/data
        """
        unittest.TestCase.__init__(self,*args, **kwargs)

        self.atom_OR_bases = utils.parse_odds_file("odds_files/atom_OR_bases.smarts" , False)
        self.atom_OR_decors = utils.parse_odds_file("odds_files/atom_decorators.smarts", False)
        self.atom_AND_decors = utils.parse_odds_file("odds_files/atom_decorators.smarts", False)
        self.bond_OR_bases = utils.parse_odds_file("odds_files/bond_OR_bases.smarts", False)
        self.bond_AND_decors = utils.parse_odds_file("odds_files/bond_AND_decorators.smarts", False)
        self.atom_odds = utils.parse_odds_file("odds_files/atom_index_odds.smarts", False)
        self.bond_odds = utils.parse_odds_file("odds_files/bond_index_odds.smarts", False)
        self.molecules = read_molecules("test_filt1_tripos.mol2", False)
        self.SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
        self.outputFile = 'test_smirky'
        replacement_file = utils.get_data_filename("odds_files/substitutions.smarts")
        self.replacements = smarty.AtomTyper.read_typelist(replacement_file)
        self.replacements = [ [short, smarts] for [smarts, short] in self.replacements]

        self.correctDict = {'VdW': [ ["[#1:1]-[#6]", 'HC'], [ "[#1:1]-[#6]-[#7,#8,F,#16,Cl,Br]", 'H1'], [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H2'], [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H3'], [ "[#1:1]-[#8]", 'HO'], [ "[#6X4:1]", 'CT'], [ "[#8X2:1]", 'OS'], [ "[#8X2+0:1]-[#1]", 'OH'] ],
                'Bond': [ ["[#6X4:1]-[#6X4:2]", 'CT-CT'], [ "[#6X4:1]-[#1:2]", 'CT-H'], [ "[#8:1]~[#1:2]", 'O~H'], [ "[#6X4:1]-[#8;X2;H1:2]", "CT-OH"], [ "[#6X4:1]-[#8;X2;H0:2]", "CT-OS"] ],
                'Angle': [ [ "[a,A:1]-[#6&X4:2]-[a,A:3]", 'any-CT-any'], [ "[#1:1]-[#6&X4:2]-[#1:3]", "H-CT-H"], [ "[#6&X4:1]-[#6&X4:2]-[#6&X4:3]", 'CT-CT-CT'], [ "[#8&X2:1]-[#6&X4:2]-[#8&X2:3]", 'O-CT-O'], [ "[#6&X4:1]-[#8&X2:2]-[#1:3]", 'CT-OH-HO'], [ "[#6X4:1]-[#8X2:2]-[#6X4:3]", 'CT-OS-CT'] ],
                'Torsion': [["[a,A:1]-[#6&X4:2]-[#6&X4:3]-[a,A:4]", 'any-CT-CT-any'], [ "[a,A:1]-[#6&X4:2]-[#8&X2:3]-[#1:4]", 'any-CT-OH-HO'], [ "[a,A:1]-[#6&X4:2]-[#8&X2:3]-[!#1:4]", 'any-CT-OS-!H'], [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[#1:4]", 'H-CT-CT-H'], [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[#6&X4:4]", 'H-CT-CT-CT'], [ "[#6&X4:1]-[#6&X4:2]-[#8&X2:3]-[#1:4]", 'CT-CT-OH-HO'], [ "[#6&X4:1]-[#6&X4:2]-[#6&X4:3]-[#6&X4:4]", 'CT-CT-CT-CT'], [ "[#6&X4:1]-[#6&X4:2]-[#8&X2:3]-[#6&X4:4]", 'CT-CT-OS-CT'], [ "[#6&X4:1]-[#8&X2:2]-[#6&X4:3]-[O&X2&H0:4]", 'CT-OS-CT-OS'], [ "[#8&X2:1]-[#6&X4:2]-[#6&X4:3]-[#8&X2:4]", 'O-CT-CT-O'], [ "[#8&X2:1]-[#6&X4:2]-[#6&X4:3]-[#1:4]", 'O-CT-CT-H'], [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[O&X2:4]", 'H-CT-CT-O'] ]}

    def test_correct_fragments(self):
        """
        Test score is 100% if correct VdW, Bond, Angles, or Torsions
        from AlkEtOH are used as input to the FragmentSampler
        """

        for typetag, initialtypes in self.correctDict.items():
            sampler = FragmentSampler(self.molecules, typetag,
                    self.atom_OR_bases, self.atom_OR_decors, self.atom_AND_decors,
                    self.bond_OR_bases, self.bond_AND_decors, self.atom_odds,
                    self.bond_odds, self.replacements, initialtypes,
                    self.SMIRFF, 0.0, self.outputFile)

            fracfound = sampler.run(2)
            self.assertAlmostEqual(fracfound, 1.0, msg = "Not finding 100%% of AlkEthOH when starting from correct %s SMIRKS." % typetag)

    def test_random_sampler(self):
        """
        Test FragmentSampler runs for 10 iterations with no failures
        Test score_utils functions with the outputFile
        """
        typetag = 'Torsion'
        sampler = FragmentSampler(self.molecules, typetag, self.atom_OR_bases,
                self.atom_OR_decors, self.atom_AND_decors, self.bond_OR_bases,
                self.bond_AND_decors, self.atom_odds, self.bond_odds,
                self.replacements, None, self.SMIRFF, 0.0, self.outputFile)
        fracfound = sampler.run(10)
        # load_trajectory converts csv file to dictionary
        timeseries = score_utils.load_trajectory('%s.csv' % self.outputFile)
        # scores_vs_time converts num/den entries to fractional scores
        scores_vs_time = score_utils.scores_vs_time(timeseries)
        # test plotting function
        score_utils.create_plot_file('%s.csv' % self.outputFile, '%s.pdf' % self.outputFile)


    def test_sampler_functions(self):
        """
        Test fragment sampler functions are working
        """
        typetag = 'Angle'
        sampler = FragmentSampler(self.molecules, typetag, self.atom_OR_bases,
                self.atom_OR_decors, self.atom_AND_decors, self.bond_OR_bases,
                self.bond_AND_decors, self.atom_odds, self.bond_odds,
                self.replacements, None, self.SMIRFF, 0.0, self.outputFile)

        typetags = [ ('VdW', 'NonbondedGenerator'),
                ('Bond', 'HarmonicBondGenerator'),
                ('Angle', 'HarmonicAngleGenerator'),
                ('Torsion', 'PeriodicTorsionGenerator'),
                ('Improper','PeriodicTorsionGenerator'),
                ('None', None)]

        for (tag, expected) in typetags:
            sample_tag, edges, sym_odds = sampler.get_type_info(tag)
            self.assertEqual(sample_tag, expected, msg = "get_force_type(%s) should return %s, but %s was returned instead" % (tag, expected, sample_tag))

        # Running each method just to make sure they work
        # get environment
        env = sampler.envList[0]
        new_env, prob = sampler.create_new_environment(env)
        # check atom methods
        atom,prob = sampler.pick_an_atom(new_env)
        removeable = sampler.isremoveable(new_env,atom)
        prob = sampler.add_atom(new_env,atom)
        prob = sampler.change_atom(new_env, atom)
        atom.addORtype('#6', ['X4'])
        prob = sampler.change_ORdecorator(atom, self.atom_OR_decors)
        prob = sampler.change_ORbase(atom, self.atom_OR_bases, self.atom_OR_decors)
        prob = sampler.change_ANDdecorators(atom, self.atom_AND_decors)

        # check bond methods
        bond,prob = sampler.pick_a_bond(new_env)
        prob = sampler.change_bond(new_env, bond)
        prob = sampler.change_ORbase(bond, self.bond_OR_bases, sampler.BondORdecorators)
        prob = sampler.change_ANDdecorators(bond, self.bond_AND_decors)

    def test_no_reference_smirff(self):
        """
        Test that sampling still works with no reference SMIRFF provided
        """
        typetag = 'Bond'
        sampler = FragmentSampler(self.molecules, typetag, self.atom_OR_bases,
                self.atom_OR_decors, self.atom_AND_decors, self.bond_OR_bases,
                self.bond_AND_decors, self.atom_odds, self.bond_odds,
                self.replacements, None, None, 0.0, self.outputFile)
        fracfound = sampler.run(10)

