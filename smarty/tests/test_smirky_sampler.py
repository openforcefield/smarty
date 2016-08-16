from smarty.environment import *
from smarty.sampler_smirky import *
from smarty.utils import *
from operator import itemgetter, attrgetter
from openeye.oechem import *
import copy
import sys # used to exit while testing

class TestSmirkySampler(TestCase):
    
    def __init__(self):
        # General things
        self.ORs = ['X1', 'X2', 'X3', 'X4']
        self.ANDs = ['+0']
        self.mol2file = 'test_filt1_tripos.mol2'
        self.SMIRFF = "forcefield/Frosst_AlkEtOH.ffxml"
        self.verbose = True
        self.elements = ["[#%i]" %i for i in range(1,119)]
        self.replacements = None
        self.molecules = read_molecules(mol2file, verbose = verbose)

    def test_correctVdW(self):
        """
        Test 100% score with correct VdW types in AlkEtOH
        """
        typetag = 'VdW'
        
        initialList = [ ["[#1:1]-[#6]]", 'HC'], 
                [ "[#1:1]-[#6]-[#7,#8,F,#16,Cl,Br]", 'H1'],
                [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H2'],
                [ "[#1:1]-[#6](-[#7,#8,F,#16,Cl,Br])(-[#7,#8,F,#16,Cl,Br])-[#7,#8,F,#16,Cl,Br]", 'H3'],
                [ "[#1:1]-[#8]", 'HO'],
                [ "[#6X4:1]", 'CT'],
                [ "[#8X2:1]", 'OS'],
                [ "[#8X2+0:1]-[#1]", 'OH'] ]

        sampler = TypeSampler(self.molecules, typetag, self.elements, self.ORs, self.ANDs, self.replacements, initialList, self.SMIRFF, 0.0, False)

        fracfound = sampler.run(2)
        if fracfound < 1.0:
            raise Exception("Not finding 100% of AlkEthOH when starting from"
                            " correct SMARTS.")
