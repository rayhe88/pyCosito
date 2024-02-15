from pyCosito.main import *
import os
import shutil
import logging

# Configurar el nivel de registro de la raíz del logger
logging.getLogger().setLevel(logging.ERROR)

# Configurar el nivel de registro de todos los loggers hijos de fireworks
logging.getLogger('fireworks').setLevel(logging.ERROR)

# Configurar el formato de los mensajes de registro
logging.basicConfig(format='%(levelname)s')



list_moles = ['LiH','BeH','CH','CH2_s1A1d','CH2_s3B1d','CH3','CH4','NH','NH2',
            'NH3','OH','H2O','HF','SiH2_s1A1d','SiH2_s3B1d','SiH3','SiH4','PH2',
            'PH3','SH2','HCl','Li2','LiF','C2H2','C2H4','C2H6','CN','HCN','CO',
            'HCO','H2CO','CH3OH','N2','N2H4','NO','O2','H2O2','F2','CO2','Na2',
            'Si2','P2','S2','Cl2','NaCl','SiO','CS','SO','ClO','ClF','Si2H6',
            'CH3Cl','CH3SH','HOCl','SO2']

conf={}

myWF = AtomizationWorkFlow(os.getcwd(), software='PWDFT', nworkers=1,
                 software_kwargs=conf, moleculeList=list_moles,
                 name="MyWorkFlow")

myWF.runWF()


files = os.listdir()

for file in files:
    if file.startswith("launcher"):
        shutil.rmtree(file)

