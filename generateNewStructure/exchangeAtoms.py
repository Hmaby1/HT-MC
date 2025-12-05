from turtle import distance
import numpy as np
import random
import os
import logging

import shutil

from pymatgen.core import Structure
from copy import deepcopy
from pymatgen.io.vasp.inputs import Poscar


from cores.sublatticeObject import StructureSublatticeObject
from logger.loggerForGenerator import LoggerForExchangeAtoms
from utilitys.mpirunContext import PwdContext

from generateNewStructure.pos_convert import poscar_convert


class ExchangeAtoms(object):
    '''
    Description
    -----------
        1. 在现有 POSCAR、CONTCAR 的基础上，交换原子，产生新的 structrue 并将其存入文件夹中
    
    Attributes
    ----------
        1. self.vasp_folder_path: str
            The folder which includes INCAR, POSCAR, POTCAR ...
        2. self.vasp_folders_path: str
            The folder which includes many vasp folders
        3. self.structure: pymatgen.core.Structure
            The current structure (before exchanging two atoms)
        4. self.structure_index: str
            The index of current structure (before exchanging two atoms)
        5. self.structure_sublattcie_object: pyMC.cores.StructureSublatticeObject
            
        6. self.log_file_path: str
            The path of the log file
    '''
    def __init__(self, poscar_path: str, sublattices_symbols_lst: list, vac_dope=False,load_CHGnet=False,vac_as="V"):
        '''
        Parameters
        ----------
            1. poscar_path: str
                可以是 POSCAR 或者是 CONTCAR，取决于你想从上一个结构的 POSCAR/CONTCAR 开始
            2. sublattice_symbols_lst:
                可以交换的 sublattice 所包含的元素
                [ ["Re", "Nb"], 
                  ["S", "Se"] ]
        '''
        self.vasp_folder_path = os.path.dirname(poscar_path)
        self.vasp_folders_path = os.path.dirname(self.vasp_folder_path)
        self.structure_index = int(os.path.split(self.vasp_folder_path)[-1])

        self.log_file_path = os.path.join(self.vasp_folder_path, "process.log")
        self.vac_dope = vac_dope
        self.load_CHGnet = load_CHGnet
        
        self.vac_as = vac_as
        if self.vac_dope:
            poscar_path = self.vac_path_deal()
        

        self.structure = Structure.from_file(poscar_path)  #from_CONTCAR
        print('5,得到路径','\n',)
        self.structure_sublattice_object = \
            StructureSublatticeObject.from_file(poscar_path=poscar_path,
                                            sublattices_symbols_lst=sublattices_symbols_lst)     #得到最大最小的原子序列，用于随机交换
        self.sublattices = sublattices_symbols_lst

    def _choose_two_atoms(self,diffusion_specie:str=None,*,with_cutoff:bool=False,cutoff:float=5.0):
        '''
        Note
        ----
            1. 选择的两个原子应该属于一个 sublattice, 并且不是同一种元素
            2. 指定一个diffusion_specie必定参与交换
        
        Return
        ------
            1. 返回两个需要交换的原子, 格式为
                （第一原子种类, 第一原子index, 第二原子种类, 第二原子index )
        '''
        #print('进行choose_two_atoms-----2-StepObject-ExchangeAtoms后'+'\n')

        all_atom_indexes_lst = list( np.arange(self.structure_sublattice_object.min_index,
                                            self.structure_sublattice_object.max_index + 1) )
        dict_index2specie = self.structure_sublattice_object.dict_index2specie()  #得到待交换原子索引字典{0:'A',1:'A',...,6:'B',...}
        all_choose = [item for sublist in self.sublattices for item in sublist]
        k = 0
        while True:
            k+=1
            print('查找次数：',k)
            try:
                first_atom_index_list,first_atom_index,first_atom_specie,dict_del_first = self.get_random_index(dict_index2specie,all_choose,diffusion_specie)
                
                if diffusion_specie is not None:
                    diffusion_specie=diffusion_specie
                elif diffusion_specie is None:
                    diffusion_specie = first_atom_specie

                if with_cutoff:           #制作cutoff新dict
                    distance_dic = self.get_all_distance(first_atom_index,dict_del_first)
                    dict2={}
                    for index,dis in distance_dic.items():
                        if dis < cutoff:      #cutoff
                            dict2[index]=dis
                else:
                    dict2=dict_del_first


                if first_atom_specie == diffusion_specie:  #检索第一第二是否有diffusion_specie
                    new_dic2=dict2
                else:
                    new_dic2={}                            
                    for index,distance in dict2.items():
                        if self.structure[index].species_string == diffusion_specie:
                            new_dic2[index]=distance
                        else:
                            pass

                second_atom_index = random.choice(list(new_dic2.keys()))
                second_atom_specie = dict_del_first[second_atom_index]
                break
            except Exception as e:
                print(f"选择失败，重新尝试。Error: {e}")
                continue  # 重新选择
        
        return first_atom_specie, first_atom_index, second_atom_specie, second_atom_index


    def get_random_index(self,dict,all_choose,diffusion_specie:str=None):
        '''
        添加代码即可指定diffusion_specie为第一交换元素
        '''
        if diffusion_specie is not None:
            diffusion_specie=diffusion_specie
        elif diffusion_specie is None:
            diffusion_specie = random.choice(all_choose)
        print(diffusion_specie)
        first_atom_index_list = [k for k, v in dict.items() if v == diffusion_specie] #指定第一选择元素即可不停交换某项原子 
        first_atom_index = random.choice(first_atom_index_list)
        dict_del_first = {k:v for k,v in dict.items() if v!=diffusion_specie} #删除1号元素序列
        return first_atom_index_list,first_atom_index,diffusion_specie,dict_del_first

    def get_all_distance(self,center_index,index_dict:dict):
        structure = deepcopy(self.structure)
        distance_dic = {}
        for index,species in index_dict.items():
            distance = structure.get_distance(center_index, index)
            distance_dic[index] = distance
        return distance_dic

    def _exchange(self,diffusion_specie:str=None,*,exchange_times:int=1,with_cutoff:bool=False):
        '''
        Return
        ------
            1. new_structure: pymatgen.core.Structure
                The new structure (after exchanging)
        '''
        new_structure = deepcopy(self.structure)
        used_indices = set()  #用于记录已经使用过的索引

        max_retries = 10
        for _ in range(exchange_times):
            retries = 0
            while retries < max_retries:
                # 选择两个原子
                first_atom_specie, first_atom_index, second_atom_specie, second_atom_index = \
                    self._choose_two_atoms(diffusion_specie,with_cutoff=with_cutoff)
                print(first_atom_index,second_atom_index,used_indices)
                # 检查索引是否重复
                if first_atom_index in used_indices or second_atom_index in used_indices:
                    #print(f"重复的索引: {first_atom_index} 或 {second_atom_index}。重新采样...")
                    retries += 1
                    continue  #重新采样
                
                # 记录使用的索引
                used_indices.add(first_atom_index)
                used_indices.add(second_atom_index)
                
                # 执行原子交换
                new_structure.replace(idx=first_atom_index, species=second_atom_specie)
                new_structure.replace(idx=second_atom_index, species=first_atom_specie) #老版本为i,新版本为idx
                print(f'交换:\n{first_atom_specie}:{first_atom_index};\n{second_atom_specie}:{second_atom_index}\n')
                
                LoggerForExchangeAtoms.log_output(level=logging.INFO,
                                                log_file_path=self.log_file_path,
                                                msg="structure{4} will exchange {0}({1}) and {2}({3}) ".format(
                                                    first_atom_specie, first_atom_index,
                                                    second_atom_specie, second_atom_index,
                                                    self.structure_index
                                                ))
                break

        return new_structure
    
    def generate_new_structure(self, new_structure_index: int, elements_str_for_vaspkit:str,pick_first_specie:str=None,with_cutoff:bool=False,exchange_times:int=1):
        '''
        Description
        -----------
            1.为新的 Structure 建立文件夹 (vasp_folder)
            2.新 vasp_folder 的index = self.index + 1
        
        Return
        ------
            1.new_poscar_path: str
                The path of new poscar_path (after exchanging)
        '''

        new_structure_index = new_structure_index           #totle steps +1 
        new_vasp_folder_path = os.path.join(self.vasp_folders_path, str(new_structure_index))
        #print('7,创建新文件夹','\n')
        try:
            os.mkdir(new_vasp_folder_path)
        except FileExistsError as e:
            LoggerForExchangeAtoms.log_output(msg=str(e),
                                            log_file_path=self.log_file_path,
                                            level=logging.ERROR)
        #print('8,写入文件 -> 6','\n')
        if self.vac_dope: #空位文件写入
            new_structure = self._exchange(pick_first_specie,with_cutoff=with_cutoff,exchange_times=exchange_times)
            new_structure,V_stru = self.Vac_del(new_structure,new_structure_index,elements_str_for_vaspkit)  #对更换位置后的poscar进行位置删除
            n_poscar = Poscar(V_stru) #有空位坐标
            n_poscar_path = os.path.join(new_vasp_folder_path, str(new_structure_index)+".POSCAR")
            n_vasp_path = os.path.join(new_vasp_folder_path, str(new_structure_index)+".vasp")
            n_poscar.write_file(n_poscar_path)
            n_poscar.write_file(n_vasp_path)
        elif not self.vac_dope:
            new_structure = self._exchange(pick_first_specie,with_cutoff=with_cutoff,exchange_times=exchange_times)
            new_structure = self.pos_sort(new_structure,elements_str_for_vaspkit)  
        #无空位文件写入
        new_poscar = Poscar(new_structure)
        new_poscar_path = os.path.join(new_vasp_folder_path, "POSCAR")
        new_poscar.write_file(new_poscar_path)

        #整理 
        if not self.vac_dope:
            pwd_path = os.getcwd()
            with PwdContext(pwd_path=pwd_path, vasp_folder_path=new_vasp_folder_path) as _:
                os.system("echo \"107\n{0}\n\" | vaspkit".format(elements_str_for_vaspkit))
                os.system("mv POSCAR POSCAR_PRI")
                os.system("mv POSCAR_REV POSCAR") 

        print('new_poscar_path:'+new_poscar_path +'\n')
        return new_poscar_path
        
    
    def Vac_del(self,new_structure,new_structure_index,elements_str_for_vaspkit): 
        '''
        根据序列号写入process文件夹中一个n.vasp用于生成POSCAR过渡文件
        '''
        process = os.path.join(self.vasp_folders_path,'process')
        if not os.path.exists(process):
            os.mkdir(process)
        nvasp = os.path.join(process,str(int(new_structure_index))+".vasp")
        vesta = os.path.join(process,str(int(new_structure_index))+"-with_vac.vasp")
        new_structure = self.pos_sort(structure=new_structure,elements_str_for_vaspkit=elements_str_for_vaspkit)
        V_stru = new_structure       #用于存于下一步n.poscar
        p = Poscar(new_structure)
        p2 = Poscar(V_stru)
        p2.write_file(vesta)    #存于process里一个交换后的n#.vasp文件 有空位
        p.write_file(nvasp)   #存于process里一个交换后的n.vasp文件  无空位
        Vac_del = poscar_convert(nvasp,'poscar')
        Vac_del.poscar_trans()       
        Vac_del_p = shutil.copy(nvasp, os.path.join(process,'POSCAR')) 
        new_structure = Structure.from_file(Vac_del_p)
        return new_structure,V_stru    
    
    def pos_sort(self,structure,elements_str_for_vaspkit):
        '''
        将元素和坐标按elements_str_for_vaspkit排序

        将空位坐标添加至末尾
        '''
        new_structure = structure
        species = new_structure.species
        coords = new_structure.frac_coords
        ele_list = elements_str_for_vaspkit.split()
        sorted_indices=[]

        for ii in range(len(ele_list)):
            sorted_indices.extend([i for i, specie in enumerate(species) 
                                if specie.symbol == ele_list[ii]])
        sorted_indices.extend([i for i, specie in enumerate(species) if self.vac_dope and specie.symbol == self.vac_as])

        sorted_species = [species[i] for i in sorted_indices]
        sorted_coords = [coords[i] for i in sorted_indices]
        
        new_structure = Structure(new_structure.lattice, sorted_species, sorted_coords)
        return new_structure
    
    def vac_path_deal(self):
        '''
        Description
        -----------
        需要准备CONTCAR\n
        利用CONTCAR -> 空位n.CONTCAR\n
        用于下一步交换\n
        return ncon_path
        '''
        con_path = os.path.join(self.vasp_folder_path,"CONTCAR")
        ncon_path = os.path.join(self.vasp_folder_path,str(self.structure_index)+".CONTCAR") 
        npos_name = str(self.structure_index)+'.POSCAR'
        npos_path = os.path.join(self.vasp_folder_path,npos_name)
        if not os.path.exists(ncon_path):
            Vac_del = poscar_convert(npos_path,'poscar')
            Vac_del.poscar_verse_trans(npos_path,con_path,ncon_path)
        return ncon_path
    
    @staticmethod #useless
    def new_contcar(npos_path,ncon_path,con_path):
        Vac_del = poscar_convert(npos_path,'poscar')
        Vac_del.poscar_verse_trans(npos_path,con_path,ncon_path)
        return con_path
    
