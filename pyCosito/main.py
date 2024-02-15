from fireworks.core.rocket_launcher import rapidfire
from fireworks.features.multi_launcher import launch_multiprocess
from fireworks import Firework, LaunchPad, ScriptTask, FWorker, Workflow, FileTransferTask
from fireworks.core.firework import FiretaskBase, FWAction
from fireworks.utilities.fw_utilities import explicit_serialize


from xtb.ase.calculator import XTB
from ase import Atoms
from ase.collections import g2
import numpy as np
import pymongo
import os

import logging

# Configurar el nivel de registro de la raíz del logger
logging.getLogger().setLevel(logging.ERROR)

# Configurar el nivel de registro de todos los loggers hijos de fireworks
logging.getLogger('fireworks').setLevel(logging.ERROR)

# Configurar el formato de los mensajes de registro
logging.basicConfig(format='%(levelname)s')

# Tu código de FireWorks aquí


# Utils functions

# insertar_atom
def insertar_atom(colection, label):
    xyz = Atoms(label, positions=[(0., 0., 0.)])
    colection.insert_one({"label":label, "symb":label,"xyz": xyz.get_positions().tolist()})

# insertar_molecule
def insertar_molecule(colection, label):
    xyz = g2[label]
    colection.insert_one({"label":label, "symb":xyz.get_chemical_symbols(), "xyz": xyz.get_positions().tolist()})

#actualizamos
def actualize_molecule(colection, label, energy):
    colection.update_one({"label":label},{"$set":{"energy":energy}})

# buscamos molecula o atomo
def localize_xyz(colection, label):
    molecule = colection.find_one({"label":label})
    if molecule:
        return molecule.get("symb"), molecule.get("xyz")
    else:
        return None
# Obtenemos la energia de la base de datos
def getEnergyfromDB(colection, label):
    system = colection.find_one({"label":label})
    if system:
        return system.get("energy")
    else:
        return None
# llenamos base de deatos
def llenamos_mole_db(colection, listName):
    for item in listName:
        insertar_molecule(colection, item)

# llenamos base de datos atomos
def llenamos_atms_db(colection, listName):
    for item in listName:
        insertar_atom(colection, item)



class EnergyTask(FiretaskBase):
    def run_task(self, fw_spec):
        raise NotImplementedError

def getEnergybyFW(mongoConf, label, stype, parents=[]):
    d = dict(mongoConf)
    d["label"] = label
    d["stype"] = stype

    t1 = EnergyTaskFW(d)
    return Firework(t1, name=label+"_energy", spec={"_priority":1})

@explicit_serialize
class EnergyTaskFW(EnergyTask):
    required_params = ['label']

    def run_task(self, fw_spec):
        host = self["host"]
        port = self["port"]
        nameDB = self["nameDB"]
        dbAtoms = self["dbAtoms"]
        dbMolec = self["dbMolec"]
        label = self["label"]
        stype = self["stype"]


        client = pymongo.MongoClient(f"mongodb://{host}:{port}/")
        db = client[nameDB]
        if stype == "ATOM":
            colection = db[dbAtoms]
        else:
            colection = db[dbMolec]

        symb, xyz = localize_xyz(colection,label)
        molecule = Atoms(''.join(symb), positions=xyz)
        molecule.calc = XTB(method='GFN1-xTB')
        molecule.set_cell([13.0, 14.0, 15.0])
        molecule.set_pbc([1,1,1])
        molecule.center()
        try:
            en = molecule.get_potential_energy()
        except:
            en = None
        print("Label : {0} energy {1}".format(label,en))
        actualize_molecule(colection, label, en)

        client.close()

        return FWAction()

class AtomizationWorkFlow:
    def __init__(self, path, software='XTB', nworkers=1,
                 software_kwargs=None, moleculeList=None, name=None, host='localhost', port=27017, nameDB='myDataBase'):
        self.path = path
        self.software = software
        self.nworkers = nworkers
        self.sofware_kwargs = software_kwargs
        self.moleculeList = moleculeList
        self.nmolecule = len(self.moleculeList)
        self.atomList = self.initAtomList()
        self.natom = len(self.atomList)
        self.launchpad = LaunchPad()
        self.launchpad.reset('', require_password=False)
        if name == None:
            self.name = "AtomizationWorkFlow"
        else:
            self.name = name

        self.host = host
        self.port = port
        self.nameDB = nameDB
        self.dbAtoms='fw_atoms'
        self.dbMolec='fw_molec'

        self.mongoConf = {"host" : self.host,
                          "port" : self.port,
                          "nameDB" : self.nameDB,
                          "dbAtoms" : self.dbAtoms,
                          "dbMolec" : self.dbMolec}

        client = pymongo.MongoClient(f"mongodb://{self.host}:{self.port}/")
        db = client[self.nameDB]
        colection_atoms = db[self.dbAtoms]
        colection_molec = db[self.dbMolec]

        clean = True
        if clean == True:
            colection_atoms.drop()
            colection_molec.drop()
            colection_atoms = db[self.dbAtoms]
            colection_molec = db[self.dbMolec]

        llenamos_mole_db(colection_atoms, self.atomList)
        llenamos_mole_db(colection_molec, self.moleculeList)

        client.close()


    def initAtomList(self):
        all_symbols = [g2[self.moleculeList[mol]].get_chemical_symbols() for mol in range(self.nmolecule)]
        flat_symbols = [symbol for symbols in all_symbols for symbol in symbols]
        unique_symbols = set(flat_symbols)

        return unique_symbols

    def GetAtomizationEnergy(self):
        client = pymongo.MongoClient(f"mongodb://{self.host}:{self.port}/")
        db = client[self.nameDB]
        colection0 = db[self.dbAtoms]
        colection1 = db[self.dbMolec]

        for doc in colection1.find():
            label = doc['label']
            energyMolecule = getEnergyfromDB(colection1, label)
            symb, xyz = localize_xyz(colection1, label)

            sumEnergyAtoms = 0.
            for s in symb:
                energyAtom = getEnergyfromDB(colection0, s)
                sumEnergyAtoms += energyAtom

            try:
                atomization = energyMolecule - sumEnergyAtoms
            except:
                atomization = None
            colection1.update_one({"label":label},{"$set":{"atomizationE":atomization}})

            client.close()

    def runWF(self):
        fws = []
        fws1 = []
        fws2 = []

        print(" ATOMS")
        for atom in self.atomList:
            task = getEnergybyFW(self.mongoConf, atom, stype="ATOM")
            fws1.append(task)

        fws.extend(fws1)

        print(" Molecules")
        for mol in self.moleculeList:
            task = getEnergybyFW(self.mongoConf, mol, stype="MOLECULE", parents=fws1)
            fws2.append(task)

        fws.extend(fws2)

        workflow = Workflow(fws, name=self.name)
        self.launchpad.add_wf(workflow)

        self.launch1()


        self.GetAtomizationEnergy()

    def launch1(self):
        if self.nworkers == 1:
            rapidfire(self.launchpad, FWorker())
        else:
            launch_multiprocess(self.launchpad, FWorker(), "Error", 0, self.nworkers, 5, local_redirect=False)

