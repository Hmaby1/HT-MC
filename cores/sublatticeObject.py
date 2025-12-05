from __future__ import annotations
import numpy as np
import random
from prettytable.prettytable import PrettyTable
from pymatgen.core import Structure
from .blankObject import BlankObject
from utilitys.formatUtilitys import Functions
from ase.neighborlist import NeighborList

class SpecieIndexesObject(object):
    
    def __init__(self, specie: str, indexes_lst: list):
       
        self.specie = specie
        self.indexes_lst = indexes_lst
        self.max_index = max(self.indexes_lst)
        self.min_index = min(self.indexes_lst)

    def __repr__(self):
        
        table = PrettyTable(["Specie", "Min_Index", "Max_Index"])
        table.add_row([self.specie, self.min_index, self.max_index])
        print(table)
        print("Indexes:")
        Functions.format_print(self.indexes_lst)
        return ""

    def __str__(self):
        return self.__repr__()

    def whether_in(self, index: int) -> bool:
        if index in self.indexes_lst:
            return True
        else:
            return False

    def dict_index2specie(self):
        index2specie = {}
        for index in self.indexes_lst:
            index2specie[index] = self.specie
        return index2specie


class SublatticeObject(object):

    def __init__(self, specie_indexes_objects_lst: list):
        self.specie_indexes_objects_lst = specie_indexes_objects_lst
        self.species_inside_sublattice_lst = [entry.specie for entry in self.specie_indexes_objects_lst]
        self.atoms_num_inside_sublattice = sum([len(specie_indexes_object.indexes_lst) \
                                    for specie_indexes_object in self.specie_indexes_objects_lst])
        self.min_index = min([specie_indexes_object.min_index \
                        for specie_indexes_object in self.specie_indexes_objects_lst])
        self.max_index = max([specie_indexes_object.max_index \
                        for specie_indexes_object in self.specie_indexes_objects_lst])

    def __repr__(self):
        
        table = PrettyTable(["Symbols in Sublattice", "The number of total atoms inside this sublattice", "Min_Index", "Max_Index"])
        table.add_row([self.species_inside_sublattice_lst, self.atoms_num_inside_sublattice, self.min_index, self.max_index])
        print(table)
        for specie_indexes_object in self.specie_indexes_objects_lst:
            print(specie_indexes_object)
        return ""

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_file(cls, poscar_path: str, species_inside_sublattice_lst: list) -> SublatticeObject:
        
        structure = Structure.from_file(poscar_path)
        all_species_lst = [specie.symbol for specie in structure.species]
        all_indexes_lst = [i for i in range(len(all_species_lst))]

        indexes_lsts_lst = []
        for i in range(len(species_inside_sublattice_lst)):
            indexes_lst = []
            for index in range(len(all_indexes_lst)):
                if (all_species_lst[index] == species_inside_sublattice_lst[i]):
                    indexes_lst.append(index)
            indexes_lsts_lst.append(indexes_lst)
        
        specie_indexes_objects_lst = []
        
        for entry in zip(species_inside_sublattice_lst, indexes_lsts_lst):
            specie_indexes_object = SpecieIndexesObject(entry[0], entry[1])  #因读取文件不正确，元素与输入不匹配容易错误
            specie_indexes_objects_lst.append(specie_indexes_object)
        
        return_object = BlankObject()
        return_object.specie_indexes_objects_lst = specie_indexes_objects_lst
        return_object.species_inside_sublattice_lst = [specie_indexes_object.specie \
                                        for specie_indexes_object in specie_indexes_objects_lst]
        return_object.atoms_num_inside_sublattice = sum([len(specie_indexes_object.indexes_lst) \
                                    for specie_indexes_object in return_object.specie_indexes_objects_lst])
        return_object.min_index = min([specie_indexes_object.min_index \
                        for specie_indexes_object in return_object.specie_indexes_objects_lst])
        return_object.max_index = max([specie_indexes_object.max_index \
                        for specie_indexes_object in return_object.specie_indexes_objects_lst])

        return_object.__class__ = cls
        
        return return_object
        
    def choose_anthor_atom(self, first_atom_index,structure,cutoff=float(5)):

        all_atom_indexes_lst = list(np.arange(self.min_index, self.max_index + 1))

        try:
            assert ( first_atom_index in all_atom_indexes_lst )
        except Exception as e:
            pass   
            raise (e)
        
        neighbors=[]
        out=[]
        while len(neighbors) == 0:
            neighbors = self.find_neighbors(first_atom_index,structure,cutoff)

            for id in neighbors:
                i=0
                for specie_indexes_object in self.specie_indexes_objects_lst:
                    if (not specie_indexes_object.whether_in(id)):
                        i+=1
                    if i == int(len(self.species_inside_sublattice_lst)):
                        out.append(id)

        for i in out:
            if i in neighbors:
                neighbors.remove(i)

        second_atom_index = random.choice(neighbors)
        
        return second_atom_index,neighbors

    def find_neighbors(self,first_atom_index,structure,cutoff=float(5)):
        ase_stru = Structure.to_ase_atoms(structure)

        cutoffs = [float(cutoff)/2] * len(ase_stru)
        nl = NeighborList(cutoffs)
        nl.update(ase_stru)
        indices, offsets = nl.get_neighbors(int(first_atom_index))

        neighbors = indices.tolist()
        
        return neighbors
    
    def whether_in(self, index: int)  -> bool:
        mark = False
        for specie_indexes_object in self.specie_indexes_objects_lst:
            if ( specie_indexes_object.whether_in(index) ):
                mark = True
        return mark
    
    def dict_index2specie(self):
        dict_index2specie = {}
        for specie_index_object in self.specie_indexes_objects_lst:
            tmp_dict_index2specie = specie_index_object.dict_index2specie()
            for key, value in tmp_dict_index2specie.items():
                dict_index2specie[key] = value
        
        return dict_index2specie
    

class StructureSublatticeObject(object):
   
    def __init__(self, sublattice_objects_lst: list):

        self.sublattice_objects_lst = sublattice_objects_lst
        self.species_inside_system_lst = [specie for sublattice_object in sublattice_objects_lst \
                                for specie in sublattice_object.species_lst]
        self.min_index = min([sublattice_object.min_index \
                        for sublattice_object in self.sublattice_objects_lst])
        self.max_index = max([sublattice_object.max_index \
                        for sublattice_object in self.sublattice_objects_lst])
    
    def __repr__(self):
        
        for sublattice_object in self.sublattice_objects_lst:
            print(sublattice_object)
        
        return ""
    
    def __str__(self):
        return self.__repr__()
    
    @classmethod
    def from_file(cls, poscar_path: str, sublattices_symbols_lst: list) -> StructureSublatticeObject:
        sublattice_objects_lst = []
        for sublattice_symbols_lst in sublattices_symbols_lst:
            sublattice_object = SublatticeObject.from_file(poscar_path = poscar_path,
                                                        species_inside_sublattice_lst = sublattice_symbols_lst)
            sublattice_objects_lst.append(sublattice_object)

        return_object = BlankObject()
        return_object.sublattice_objects_lst = sublattice_objects_lst
        return_object.species_inside_system_lst = [specie for sublattice_object in sublattice_objects_lst \
                                            for specie in sublattice_object.species_inside_sublattice_lst]
        return_object.min_index = min([sublattice_object.min_index \
                        for sublattice_object in return_object.sublattice_objects_lst])
        return_object.max_index = max([sublattice_object.max_index \
                        for sublattice_object in return_object.sublattice_objects_lst])
        return_object.__class__ = cls

        return return_object

    def dict_index2specie(self):
        dict_index2specie = {}
        for sublattice_object in self.sublattice_objects_lst:
            tmp_dict_index2specie = sublattice_object.dict_index2specie()
            for key, value in tmp_dict_index2specie.items():
                dict_index2specie[key] = value
        
        return dict_index2specie