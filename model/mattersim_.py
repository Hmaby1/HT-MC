from typing import Union

import torch
import sys
import os
import numpy as np
from io import StringIO

from ase.build import bulk
from ase.units import GPa
from ase import Atom, io

from mattersim.forcefield import MatterSimCalculator
from mattersim.applications.relax import Relaxer
from mattersim.forcefield.potential import Potential
from mattersim.datasets.utils.build import build_dataloader

from chgnet.model import StructOptimizer
from pymatgen.core import Structure

class mattersim_predict():
    def __init__(self,
                 poscar_path,
                 device:str="cuda" if torch.cuda.is_available() else "cpu",
                 load_path:str=None,
                 ):
        self.poscar_path = poscar_path
        self.folder_path = os.path.split(poscar_path)[0]
        self.contcar_path = contcar_path=os.path.join(self.folder_path,'CONTCAR')
        self.output_path = os.path.join(self.folder_path,'relaxation_output.txt')
        self.device = device
        self.load_path = load_path
        self.atom = io.read(poscar_path)

    @classmethod
    def load(self,poscar_path,
             load_path:str=None,
             device:str="cuda" if torch.cuda.is_available() else "cpu"):
        mattersim = mattersim_predict(poscar_path,device,load_path)
        if os.path.exists(mattersim.contcar_path):
            atom = io.read(mattersim.contcar_path)
        else:
            atom = mattersim.relax(mattersim.load_path)
        energy = mattersim.predict(atom,mattersim.load_path)
        mattersim.contcar_atom = atom
        mattersim.energy = energy
        return mattersim
    
    def relax(self,
              load_path=None,
              perturb:float=0.01,
              relax_step:int=500,
              fmax:float=0.01,
              *,
              optimizer:str="FIRE",
              filter:str="FrechetCellFilter",
              constrain_symmetry:bool=True,
              relax_model:str='mattersim',
              max_retries:int=5) -> Atom :
        if relax_model=='mattersim':
            attempt=0
            while True:
                try:
                    # mattersim_relax
                    contcar_path = self.contcar_path
                    contcar = self.atom
                    contcar.positions += perturb * 10 * np.random.randn(len(contcar), 3)#
                    potential=Potential.from_checkpoint(load_path=load_path,device=self.device)
                    contcar.calc = MatterSimCalculator(potential=potential,device=self.device)

                    relaxer = Relaxer(
                       optimizer=optimizer,
                       filter=filter,
                       constrain_symmetry=constrain_symmetry,
                    )
                    
                    original_stdout = sys.stdout
                    output_buffer = StringIO()
                    sys.stdout = output_buffer

                    relax_result = relaxer.relax(contcar, steps=relax_step, fmax=fmax) 
                    
                    sys.stdout = original_stdout  
                    relaxed_atoms = relax_result[1]
                    break  
                except Exception as e:
                    sys.stdout = original_stdout  
                    print(f"弛豫过程失败，错误信息：{e}")
                    attempt += 1
                    
        # except:
        elif relax_model=='chgnet':
            #chgnet_relax
            structure = Structure.from_file(self.poscar_path)
            relaxer = StructOptimizer()
            structure.perturb(perturb)
            original_stdout = sys.stdout
            output_buffer = StringIO()
            sys.stdout = output_buffer
            result = relaxer.relax(structure,steps=relax_step, fmax=fmax, verbose=True)
            sys.stdout = original_stdout
            relaxed_structure=result["final_structure"]
            relaxed_atoms=Structure.to_ase_atoms(relaxed_structure)

        with open(self.output_path, "a+") as f:
            f.write(output_buffer.getvalue())  
            self.atom_save(relaxed_atoms,self.contcar_path)
            energy=self.predict(relaxed_atoms,load_path)
            f.close()

        return relaxed_atoms
    
    def predict(self,
                Atom:Atom,
                load_path:str=None,
                property:list[str]=['energy','energy_per_atom']
                ) -> float:
        '''property=['energy','energy_per_atom','']'''
        if load_path:
            potential=Potential.from_checkpoint(load_path=load_path,device=self.device)
        else:
            potential=None
        
        atom = Atom
        atom.calc = MatterSimCalculator(potential=potential,device=self.device)
        energy = atom.get_potential_energy()
        energy_per_atom = energy/len(atom)
        with open(self.output_path, "a+") as f:
            f.write(f'\nthe final structure energy:{energy}')
            f.close()


        return energy
    
    def atom_save(self,atom:Atom,save_path:str):
        io.write(filename=save_path,images=atom)
