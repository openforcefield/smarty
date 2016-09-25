from functools import partial
import smarty
from smarty import AtomTyper, AtomTypeSampler
from smarty.sampler_elemental import *
from smarty.utils import get_data_filename
from unittest import TestCase

class TestAtomTypeSampler(TestCase):
    def test_atomtyper(self):
        basetypes_filename = get_data_filename('atomtypes/basetypes.smarts')
        initialtypes_filename = get_data_filename('atomtypes/basetypes.smarts')
        decorators_filename = get_data_filename('atomtypes/decorators.smarts')
        replacements_filename = get_data_filename('atomtypes/replacements.smarts')
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-tripos.mol2.gz'), verbose=False)
        reference_typed_molecules = smarty.utils.read_molecules(get_data_filename('molecules/zinc-subset-parm@frosst.mol2.gz'), verbose=False)

        # Construct atom type sampler.
        atomtype_sampler = smarty.AtomTypeSampler(molecules, basetypes_filename, initialtypes_filename, decorators_filename, replacements_filename=replacements_filename, reference_typed_molecules=reference_typed_molecules, verbose=False, decorator_behavior='simple-decorators')

        # Start sampling atom types.
        atomtype_sampler.run(2)

    def test_atomtyper_AlkEthOH(self):
        # Test to make sure that AlkEthOH gets 100% at T=0 starting from correct
        # types
        basetypes_filename = get_data_filename('atomtypes/basetypes.smarts')
        initialtypes_filename = get_data_filename('atomtypes/initial_AlkEthOH.smarts')
        decorators_filename = get_data_filename('atomtypes/new-decorators.smarts')
        replacements_filename = get_data_filename('atomtypes/replacements.smarts')
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'), verbose=False)
        reference_typed_molecules = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_ff.mol2'), verbose=False)

        # Construct atom type sampler.
        atomtype_sampler = smarty.AtomTypeSampler(molecules, basetypes_filename, initialtypes_filename, decorators_filename, replacements_filename=replacements_filename, reference_typed_molecules=reference_typed_molecules, verbose=False, temperature = 0)

        # Start sampling atom types.
        fracfound = atomtype_sampler.run(2)

        # Ensure fraction found is 1.0
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct SMARTS.")

    def test_atomtyper_elemental(self):
        # Test to make sure that sampler Elemental gets 100% of the Carbons on AlkEthOH at T=0 starting from correct
        # types
        basetypes_filename = get_data_filename('atomtypes/basetypes.smarts')
        initialtypes_filename = get_data_filename('atomtypes/initial_AlkEthOH.smarts')
        decorators_filename = get_data_filename('atomtypes/new-decorators.smarts')
        replacements_filename = get_data_filename('atomtypes/replacements.smarts')
        molecules = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_tripos.mol2'), verbose=False)
        reference_typed_molecules = smarty.utils.read_molecules(get_data_filename('molecules/AlkEthOH_test_filt1_ff.mol2'), verbose=False)
        
        # Construct atom type sampler.
        atomtype_sampler = smarty.AtomTypeSamplerElemental(molecules, basetypes_filename, initialtypes_filename, decorators_filename, replacements_filename=replacements_filename, reference_typed_molecules=reference_typed_molecules, verbose=True, temperature = 0, element = "6")
        
        # Start sampling atom types.
        fracfound = atomtype_sampler.run(2)
        
        # Ensure fraction found is 1.0
        if fracfound < 1.0:
            raise Exception("Not finding 100% of Carbons of AlkEthOH when starting from"
                            " correct SMARTS.")

