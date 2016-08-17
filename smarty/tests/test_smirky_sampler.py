from unittest import TestCase
from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *
from operator import itemgetter, attrgetter
from openeye.oechem import *
import copy
import sys # used to exit while testing
# General things
ORs = ['X1', 'X2', 'X3', 'X4']
ANDs = ['+0']
mol2file = 'molecules/test_filt1_tripos.mol2'
SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
elements = ["[#%i]" %i for i in range(1,119)]
replacements = None
molecules = read_molecules(get_data_filename(mol2file), verbose = False)

class TestSmirkySampler(TestCase):
    def test_correctVdW(self):
        """
        Test 100% score with correct VdW types in AlkEtOH
        """
        typetag = 'VdW'
        
        initialList = [ ["[#1:1]-[#6]", 'HC'], 
                [ "[#1:1]-[#6]-[#7,#8,F,#16,Cl,Br]", 'H1'],
                [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H2'],
                [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H3'],
                [ "[#1:1]-[#8]", 'HO'],
                [ "[#6X4:1]", 'CT'],
                [ "[#8X2:1]", 'OS'],
                [ "[#8X2+0:1]-[#1]", 'OH'] ]
        sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, False)

        fracfound = sampler.run(2)
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct VdW SMIRKS.")

    def test_correctBonds(self):
        """
        Test 100% score with correct Bond types in AlkEtOH
        """
        initialList = [ ["[#6X4:1]-[#6X4:2]", 'CT-CT'],
                [ "[#6X4:1]-[#1:2]", 'CT-H'],
                [ "[#8:1]~[#1:2]", 'O~H'],
                [ "[#6X4:1]-[#8;X2;H1:2]", "CT-OH"],
                [ "[#6X4:1]-[#8;X2;H0:2]", "CT-OS"] ] 

        typetag = 'Bond'
        sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, False)

        fracfound = sampler.run(2)
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct Bond SMIRKS.")

    def test_correctAngles(self):
        """
        Test 100% score with correct Angles types in AlkEtOH
        """
        initialList = [ [ "[a,A:1]-[#6&X4:2]-[a,A:3]", 'any-CT-any'], 
                [ "[#1:1]-[#6&X4:2]-[#1:3]", "H-CT-H"],
                [ "[#6&X4:1]-[#6&X4:2]-[#6&X4:3]", 'CT-CT-CT'],
                [ "[#8&X2:1]-[#6&X4:2]-[#8&X2:3]", 'O-CT-O'],
                [ "[#6&X4:1]-[#8&X2:2]-[#1:3]", 'CT-OH-HO'],
                [ "[#6X4:1]-[#8X2:2]-[#6X4:3]", 'CT-OS-CT'] ]

        typetag = 'Angle'
        sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, False)

        fracfound = sampler.run(2)
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct Angle SMIRKS.")

    def test_correctTorsions(self):
        """
        Test 100% score with correct Torsions types in AlkEtOH
        """
        initialList = [["[a,A:1]-[#6&X4:2]-[#6&X4:3]-[a,A:4]", 'any-CT-CT-any'], 
                [ "[a,A:1]-[#6&X4:2]-[#8&X2:3]-[#1:4]", 'any-CT-OH-HO'],
                [ "[a,A:1]-[#6&X4:2]-[#8&X2:3]-[!#1:4]", 'any-CT-OS-!H'],
                [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[#1:4]", 'H-CT-CT-H'],
                [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[#6&X4:4]", 'H-CT-CT-CT'],
                [ "[#6&X4:1]-[#6&X4:2]-[#8&X2:3]-[#1:4]", 'CT-CT-OH-HO'],
                [ "[#6&X4:1]-[#6&X4:2]-[#6&X4:3]-[#6&X4:4]", 'CT-CT-CT-CT'],
                [ "[#6&X4:1]-[#6&X4:2]-[#8&X2:3]-[#6&X4:4]", 'CT-CT-OS-CT'],
                [ "[#6&X4:1]-[#8&X2:2]-[#6&X4:3]-[O&X2&H0:4]", 'CT-OS-CT-OS'],
                [ "[#8&X2:1]-[#6&X4:2]-[#6&X4:3]-[#8&X2:4]", 'O-CT-CT-O'],
                [ "[#8&X2:1]-[#6&X4:2]-[#6&X4:3]-[#1:4]", 'O-CT-CT-H'],
                [ "[#1:1]-[#6&X4:2]-[#6&X4:3]-[O&X2:4]", 'H-CT-CT-O'] ]
        typetag = 'Torsion'
        sampler = TypeSampler(molecules, typetag, elements, ORs, ANDs, replacements, initialList, SMIRFF, 0.0, False)

        fracfound = sampler.run(2)
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct Torsion SMIRKS.")


