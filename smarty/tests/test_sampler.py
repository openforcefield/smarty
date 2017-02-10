from functools import partial
import smarty
from smarty import AtomTyper, AtomTypeSampler, score_utils
from smarty.sampler_elemental import *
from smarty.utils import get_data_filename
import unittest
from unittest import TestCase

class TestAtomTypeSampler(TestCase):
    def __init__(self, *args, **kwargs):
        """
        Initialize TestCase including files used in all smarty tests
        """
        unittest.TestCase.__init__(self, *args, **kwargs)

        self.basetypes = get_data_filename('atomtypes/basetypes.smarts')
        self.alkethoh_answers = get_data_filename('atomtypes/initial_AlkEthOH.smarts')
        self.simple_decs = get_data_filename('atomtypes/decorators.smarts')
        self.combine_decs = get_data_filename('atomtypes/new-decorators.smarts')
        self.replacements = get_data_filename('atomtypes/replacements.smarts')

        # import molecules
        self.mols_zinc = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-tripos.mol2.gz'), verbose=False)
        self.mols_zinc_ref = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-parm@frosst.mol2.gz'), verbose=False)

        self.mols_alkethoh = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'), verbose=False)
        self.mols_alkethoh_ref = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_ff.mol2'), verbose=False)


    def test_atomtyper(self):
        """
        Test atomtype sampler with simple-decorators
        """
        atomtype_sampler = smarty.AtomTypeSampler(self.mols_zinc,
                self.basetypes, self.basetypes, self.simple_decs,
                self.replacements, self.mols_zinc_ref, 0.1, False,
                'simple-decorators')
        atomtype_sampler.run(2)

    def test_atomtyper_combinatorial(self):
        """
        Test atomtype sampler with combinatorial-decorators and optional output files
        """
        atomtype_sampler = smarty.AtomTypeSampler(self.mols_zinc,
                self.basetypes, self.basetypes, self.combine_decs,
                self.replacements, self.mols_zinc_ref, 0.1, False)
        # run sampler with optional outputs
        traj = 'test_smarty.csv'
        plot = 'test_smarty.pdf'
        atomtype_sampler.run(5, traj, plot)
        # test trajectory analysis functions on smarty output
        timeseries = score_utils.load_trajectory(traj)
        scores_vs_time = score_utils.scores_vs_time(timeseries)
        score_utils.create_plot_file(traj, plot, True, False)

    def test_atomtyper_elemental(self):
        """
        Test elemental atomtype sampler for hydrogen
        """
        atomtype_sampler = smarty.AtomTypeSamplerElemental(self.mols_alkethoh,
                self.basetypes, self.basetypes, self.combine_decs,
                self.replacements, self.mols_alkethoh_ref, 0.1, False,
                'combinatorial-decorators', '1')
        # run sampler with optional outputs
        traj = 'test_smarty.csv'
        plot = 'test_smarty.pdf'
        atomtype_sampler.run(5, traj, plot)
        # test trajectory analysis functions on smarty output
        timeseries = score_utils.load_trajectory(traj)
        scores_vs_time = score_utils.scores_vs_time(timeseries)
        score_utils.create_plot_file(traj, plot, True, False)

    def test_atomtyper_AlkEthOH(self):
        """
        Test atomtype sampler with correct "answers"
        """
        atomtype_sampler = smarty.AtomTypeSampler(self.mols_alkethoh,
                self.basetypes, self.alkethoh_answers, self.combine_decs,
                self.replacements, self.mols_alkethoh_ref, 0.1, False)
        # Start sampling atom types.
        fracfound = atomtype_sampler.run(2)
        # Ensure fraction found is 1.0
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct SMARTS.")

    def test_atomtyper_elemental_AlkEthOH(self):
        """
        Test elemental sampler with correct "answers"
        """
        atomtype_sampler = smarty.AtomTypeSamplerElemental(self.mols_alkethoh,
                self.basetypes, self.alkethoh_answers, self.combine_decs,
                self.replacements, self.mols_alkethoh_ref, 0.1, False,
                'combinatorial-decorators','1')
        # Start sampling atom types.
        fracfound = atomtype_sampler.run(2)

        # Ensure fraction found is 1.0
        if fracfound < 1.0:
            raise Exception("Not finding 100% of Hydrogens of AlkEthOH when starting from"
                            " correct SMARTS.")

