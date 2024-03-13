#!/usr/bin/env python

#Jared Adolf-Bryfogle (jadolfbr@gmail.com)

import gzip
import os
import re
import logging
from pathlib import Path
from typing import Union, DefaultDict, Dict, Tuple, List, Any, ClassVar

from collections import defaultdict

from Bio.PDB import MMCIFParser
from Bio.PDB import PDBParser
from Bio.PDB import calc_dihedral
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.Polypeptide import PPBuilder

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from jade2.basic.restype_definitions import RestypeDefinitions
from jade2.basic.structure.PDBResInfo import PDBResInfo as PDBInfo
from jade2.basic.structure.PDBResidue import PDBResidue as ResidueRecord
from jade2.basic.numeric import *
from jade2.basic import vector1
from jade2.basic.path import *
from jade2.basic.structure.biopython_util import is_connected_to_prev, is_connected_to_next

class PoseResidue(Residue):
    def __init__(self, res: Residue, chain_id: str):
        self.id: List[Any] = []
        Residue.__init__(self, res.id, res.resname, res.segid)
        self.resnum = res.id[1]
        self.chain_id = chain_id

class BioPose(object):
    """
    This is my biopython meta class.  Because biopython's interface kinda sucks.
    This is a little cleaner.

    The other way is to sublclass each Biopython class structure, which I'm not ready to do.

    Right now, you need a path as I don't know how we would use this from sequence, etc as you do in Rosetta.
    :path: Is a path to an RCSB file.  PDB (.pdb), mmCIF(.cif), and gzipped (.gz) versions.
    """
    def __init__(self, path: Union[Path, str], model_num: int = 0):

        self.res_definitions: ClassVar[RestypeDefinitions] = RestypeDefinitions()
        self.name = os.path.basename(str(path))
        self.in_path = str(path)

        ####Setup Class###:
        struct_header:      Tuple[Structure, Dict] = self.load_from_file(path) #Bio struct, Header dictionary

        self.struct:        Structure =              struct_header[0]
        self.header:        Dict  =                  struct_header[1]
        self.all_residues:  vector1 =                self._setup_all_residues(model_num)
        self.pdb_info:      PDBInfo =                self._setup_pdb_info()


    ############ IO ###################
    def load_from_file(self, path: Union[Path, str]) -> Tuple[Structure, Dict]:
        """
        Load a file from PDB or mmCIF.  .gz is supported.

        :param path: Path to PDB or mmCIF file
        :rtype: tuple(bio.PDB.Structure.Structure, dict)
        """
        structure = None

        header = ""
        if re.search(".pdb", str(path)):
            parser = PDBParser()
            header = parser.get_header()
        else:
            #No header function for MMCIF parser
            parser = MMCIFParser()

        path = str(path).strip()
        model_id = self.name
        if os.path.basename(str(path)).split('.')[-1] == 'gz':
            logging.info("Opening Gzipped file")
            GZ = gzip.open(path, 'rt')
            structure = parser.get_structure(model_id, GZ)
            GZ.close()
        else :
            structure = parser.get_structure(model_id, path)



        return (structure, header)

    def reload_from_file(self, path: Union[Path, str], model_num: int =0):
        """
        Reload a BioPose from a file path.
        :param path: str
        :param model_num: int
        :return:
        """

        ####Setup Class###:
        struct_header: Tuple[Structure, Dict] = self.load_from_file(path)  # Bio struct, Header dictionary

        self.struct: Structure = struct_header[0]
        self.header: Dict = struct_header[1]
        self.all_residues: vector1 = self._setup_all_residues(model_num)
        self.pdb_info: PDBInfo = self._setup_pdb_info()



    ############ Getting Components #############
    def structure(self) -> Structure:
        """
        Get the Bio Structure stored in this class.
        :rtype: bio.PDB.Structure.Structure
        """
        return self.struct

    def model(self, model_num: int = 0) -> Model:
        """
        Get a Bio Model of the stored structure
        :rtype: bio.PDB.Model.Model
        """
        return self.struct[model_num]

    def chain(self, chain_id: str, model_num: int = 0) -> Chain:
        """
        Get a Bio Chain of the stored structure
        :rtype: bio.PDB.Chain.Chain
        """
        return self.struct[model_num][chain_id]

    def residue(self, resnum: int, chain_id: str, icode: str =' ', alt:str =' ', model_num: int = 0) -> PoseResidue:
        """
        Get a Bio Residue of the stored structure.
        Adds a chain_id attribute.

        :rtype: bio.PDB.Residue.Residue
        """
        res = self.struct[model_num][chain_id][(alt, resnum, icode)]

        #This is a bit hacky for the type checker, but I don't see a way to construct
        # a biopython Residue from a Residue and have all the internal Entity data intact.
        res.chain_id = chain_id
        return res


    def atom(self, atom_name: str, resnum: int, chain_id: str, icode:str =' ', alt:str =' ', model_num: int =0) -> Atom:
        """
        Get a Bio Atom of the stored structure

        :rtype: bio.PDB.Atom.Atom
        """
        return self.struct[model_num][chain_id][(alt, resnum, icode)][atom_name]


    ### Lists of Structures ###

    def models(self):
        """
        Get a list of Bio models
        :rtype: list[bio.PDB.Chain.Chain
        """
        return [c for c in self.struct]

    def chains(self, model_num=0) -> List[Chain]:
        return [c for c in self.struct[model_num].get_chains()]

    def residues(self, chain_id: str, model_num: int = 0, include_alt: bool = False) -> List[PoseResidue]:
        """
        Get residues, including or not including residues with alternate location codes - which can be a PITA
        Adds chain_id attribute to residue.

        :rtype: list[bio.PDB.Residue.Residue]
        """
        resi = []
        for res in self.chain(chain_id, model_num):
            res.chain_id = chain_id
            if res.id[0] ==' ':
                resi.append(res)
            elif include_alt:
                resi.append(res)
            else:
                continue
        return resi

    def atoms(self, resnum: int, chain_id: str, icode:str =' ', alt:str =' ', model_num: int = 0) -> List[Atom]:
        """
        Get a list of Bio Atoms
        :rtype: list[bio.PDB.Atom.Atom]
        """
        return [atm for atm in self.residue(resnum, chain_id, icode, alt, model_num)]

    def name(self) -> str:
        """
        Return the basename used to construct this BioPose.
        Same as model_id for bio pdb object.
        :return:
        """
        return self.name

    def path(self) -> str:
        """
        Return the input path used to construct the BioPose
        :return:
        """
    ############ Helper Funtions ################
    def total_residue(self) -> int:
        return len(self.all_residues)

    def phi(self, res: PoseResidue, resprev: PoseResidue) -> float:
        """
        Get the Phi Angle of i in radians
        """
        resatomlist = [i.name for i in
                       self.atoms(resnum=res.id[1], chain_id=res.chain_id, icode=res.id[2], alt=res.id[0])]
        resprevatomlist = [i.name for i in
                           self.atoms(resnum=resprev.id[1], chain_id=resprev.chain_id, icode=resprev.id[2],
                                      alt=resprev.id[0])]

        if not is_connected_to_prev(res, resprev):
            return 0.0
        if not all(x in resatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0
        if not all(x in resprevatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
            cp = resprev['C'].get_vector()
            phi = calc_dihedral(cp, n, ca, c)
            return phi
        except Exception:
            print("Could not get phi for " + repr(res))
            raise LookupError



    def psi(self, res: PoseResidue, resnext: PoseResidue) -> float:
        """
        Get the Psi Angle of i in radians
        """
        resatomlist = [i.name for i in self.atoms(resnum=res.id[1], chain_id=res.chain_id, icode=res.id[2], alt=res.id[0])]
        resnextatomlist = [i.name for i in self.atoms(resnum=resnext.id[1], chain_id=resnext.chain_id, icode=resnext.id[2], alt=resnext.id[0])]

        if not is_connected_to_next(res, resnext):
            return 0.0
        if not all(x in resatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0
        if not all(x in resnextatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
            nn = resnext['N'].get_vector()
            psi = calc_dihedral(n, ca, c, nn)
            return psi
        except Exception:
            print("Could not get psi for "+repr(res))
            raise LookupError


    def omega(self, res: PoseResidue, resnext: PoseResidue, resprev: PoseResidue, rosetta_definitions: bool = True) -> float:
        """
        Get the Omega Angle of i in radians
        Omega is defined as the dihedral angle between the peptide bond of i and i + 1, as in Rosetta.
        If rosetta_definitions are False, omega is then treated as being between i and i -1
        :rtype: float
        """
        resatomlist = [i.name for i in self.atoms(resnum=res.id[1], chain_id=res.chain_id, icode=res.id[2], alt=res.id[0])]
        resprevatomlist = [i.name for i in self.atoms(resnum=resprev.id[1], chain_id=resprev.chain_id, icode=resprev.id[2], alt=resprev.id[0])]
        resnextatomlist = [i.name for i in self.atoms(resnum=resnext.id[1], chain_id=resnext.chain_id, icode=resnext.id[2], alt=resnext.id[0])]
        if not is_connected_to_prev(res, resprev):
            return 0.0
        if not all(x in resatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0
        if not all(x in resprevatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0
        if not all(x in resnextatomlist for x in ['N', 'C', 'CA', 'O']):
            return 0.0

        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()

            if rosetta_definitions and is_connected_to_next(res, resnext):
                next_n = resnext['N'].get_vector()
                next_ca = resnext['CA'].get_vector()
                omega = calc_dihedral(ca, c, next_n, next_ca)
                return omega

            elif not rosetta_definitions and is_connected_to_prev(res, resprev):
                pre_c = resprev['C'].get_vector()
                pre_ca = resprev['CA'].get_vector()
                omega = calc_dihedral(pre_ca, pre_c, n, ca)
                return omega
            else:
                return 0.0

        except BaseException:
            print("Could not get omega for "+repr(res))
            raise LookupError

    def get_sequence(self, chain_id: str, model_num: int = 0) -> str:
        """
        Get a sequence of a chain - Not including alternate res locations

        :param chain_id: str
        :param model_num: int
        :rtype: str
        """
        if self.get_chain_length(chain_id, model_num) == 0:
            return ""

        ppb = PPBuilder(radius=float(1000))
        seq = ""
        chainmodel = ppb.build_peptides(self.struct, aa_only=False)
        for aa in chainmodel:
            seq = seq + aa.get_sequence()

        return seq

    def get_pose_sequence(self, model_num=0) -> str:
        pdb_seq = "".join(
            [three_to_one(res.get_resname()) for res in struct.get_residues()]
        )
        return pdb_seq

    def find_chain_breaks(self) -> Dict[Any, List[int]]:
        """
        Author: Simon Kelow

        Returns a list of missing residues for each chain id.
        """

        chainbreak_dict = dict()

        for chain in self.model():
            missinglist = []
            reslist = [i for i in chain.get_residues() if is_aa(i, standard=False)]
            for res in reslist[10:-10]:
                resatomlist = [i.name for i in self.atoms(resnum=res.id[1], chain_id=res.chain_id, icode=res.id[2], alt=res.id[0])]
                if is_aa(res, standard=False) and all(x in resatomlist for x in ['N', 'C', 'CA', 'O']):
                    resindex = reslist.index(res)
                    resnextindex = resindex + 1
                    resnext = reslist[resnextindex]
                    resnextatomlist = [i.name for i in self.atoms(resnum=resnext.id[1], chain_id=resnext.chain_id, icode=resnext.id[2], alt=resnext.id[0])]
                    resprevindex = resindex - 1
                    resprev = reslist[resprevindex]
                    resprevatomlist = [i.name for i in self.atoms(resnum=resprev.id[1], chain_id=resprev.chain_id, icode=resprev.id[2], alt=resprev.id[0])]
                    if is_aa(resnext, standard=False) and is_aa(resprev, standard=False) and all(x in resnextatomlist for x in ['N', 'C', 'CA', 'O']) and all(x in resprevatomlist for x in ['N', 'C', 'CA', 'O']):
                        connectedtonext = is_connected_to_next(res, resnext)
                        connectedtoprev = is_connected_to_prev(res, resprev)
                        if not connectedtonext:
                            missinglist.append(res.id[1])
                        if not connectedtoprev:
                            missinglist.append(res.id[1])
                    else:
                        missinglist.append(res.id[1])
            chainbreak_dict[chain.id] = missinglist
        return chainbreak_dict

    def get_regional_sequence(self, start: int, end: int) -> str:
        seq = ""
        for i in range(start, end+1):
            aa = self.res_definitions.get_one_letter_from_three(self.all_residues[i].get_resname())
            if aa:
                seq = seq+aa
            else:
                seq = seq+'X'
        return seq

    def get_chain_length(self, chain_id: str, model_num: int = 0) -> int:
        """
        Get the number of AA in a chain - Not including alternate res locations
        """
        return len(self.residues(chain_id, model_num))

    def get_chain_ids(self, model_num: int) -> List[Any]:
        """
        Get all chain IDS for a model.
        :rtype: list[str]
        """
        ids = []
        for chain in self.model(model_num):
            ids.append(chain.id)
        return ids

    ###### Private Members #######

    def _setup_all_residues(self, model_num: int =0) -> vector1:
        """
        Setup all residues in the model, for ease of use like a pose.
        :param model_num: int
        :rtype: list[bio.PDB.Residue.Residue]
        """
        all_residues = vector1()

        for chain_id in self.get_chain_ids(model_num):
            logging.info("ChainID: "+chain_id)
            residues = self.residues(chain_id, model_num)
            all_residues.extend(residues)

        return all_residues

    def _setup_pdb_info(self) -> PDBInfo:
        """
        Setup the PDBInfo mapping from residues to pdb iformation.
        :param model_num: int
        :rtype: PDBInfo
        """
        pdb_info = PDBInfo()
        i = 1
        for res in self.all_residues:
            s = ResidueRecord(self.res_definitions.get_one_letter_from_three(res.resname), res.id[1], res.chain_id, res.id[2])
            pdb_info.add_residue(s)

        return pdb_info

########## Testing Functions #############
def test_dihedrals(pose):
    """
    Simple Test for Dihedral output
    THIS NO LONGER WORKS (JAB)
    :param pose: BioPose
    :rtype: bool
    """
    for i in range(1, pose.total_residue()+1):

        print("\n"+str(pose.pdb_info.pose2pdb(i)))
        try:
            print("Phi: "+repr(math.degrees(pose.phi(i))))
            print("Psi: "+repr(math.degrees(pose.psi(i))))
            print("Omega:"+repr(math.degrees(pose.omega(i))))
        except Exception:
            print("Could not get dihedral for resnum "+repr(i))

    return True

def test_pdbinfo(pose):
    """
    Simple Test for pdbinfo output.
    :param pose: BioPose
    :rtype: bool
    """
    for i in range(1, pose.total_residue() +1):
        print(repr(i))
        print(pose.all_residues[i].id)
        print(pose.pdb_info.pose2pdb(i))


######## Testing Main ############
if __name__ == "__main__":
    v = vector1([1, 2, 3])
    for i in v:
        print(repr(i))

    print(len(v))

    print(repr(v[3]))

    test_pdb = os.path.join(get_testing_inputs_path(),"2j88.pdb")
    pose = BioPose(test_pdb)
    test_dihedrals(pose)






