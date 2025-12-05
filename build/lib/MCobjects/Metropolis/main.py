'''
Descripttion: 
version: 
Author: sch
Date: 2022-03-29 19:11:07
LastEditors: sch
LastEditTime: 2022-03-31 08:54:44
'''
from doctest import FAIL_FAST
from tracemalloc import start
from MCobjects.Metropolis.strategy import Exchange
from cores.stepObject import StructureState
from cores.stepObject import StepObject
from calculators.vaspCalculators import VaspTask
from model.mattersim_ import mattersim_predict

import time
import os
import sys
project_root = os.path.dirname(os.path.abspath(__file__))
sys.path.append('d:\桌面\MC')

class Metropolis(object):
    '''
    PBS VASP as calculator 
    '''
    def __init__(self, pbs_nodefile:str, np:int,
                dxec:str):
        self.pbs_nodefile = pbs_nodefile
        self.np = np
        self.dxec = dxec


    def run(self, poscar_path:str, num_loops:int, T:float,
            sublattice_symbols_lst:list,load_CHGnet=False,load_path=None,load_model:str=None,
            from_contcar=True,
            elements_str_for_vaspkit:str=None,
            load=False,
            *,
            vac_dope = False,
            vac_as = 'V',
            open_diffusion = False,
            diffusion_specie:str=None,
            time_save:bool=True):
        '''
        Parameters
        ----------
            1. poscar_path: str
                起始 POSCAR 的路径
            2. sublattice_symbols_lst: list
                [ ["Re", "Nb"],
                  ["S", "Se"] ]
        
        Note
        ----
            1. elements_str_for_vaspkit: str
                同一个亚晶格的元素的原子必须紧挨着，e.g. 
                    - "Re Nb S Se"  ( correct )
                    - "S Se Nb Re"  ( corrent )
                    - "Re S Se Ta"  ( wrong )
        '''
        assert (elements_str_for_vaspkit is not None)
        #补丁,对初始文件进行计算
        if load_model=='chgnet':
            load_CHGnet=True
        print('一','\n')
        structure_state = StructureState(poscar_path=poscar_path,
                                         vac_dope=vac_dope,
                                         load_CHGnet=load_CHGnet,
                                         load_path=load_path)
        
        if load_model=='chgnet':
            structure_state.load_model(load_CHGnet = load_CHGnet,
                                        load_path=load_path)
        elif load_model=='mattersim':
            mattersim_predict.load(structure_state.poscar_path,
                                   load_path=load_path)

        print('二','\n')
        step_object = StepObject(current_structure_state=structure_state,
                                sublattice_symbols_lst=sublattice_symbols_lst,
                                elements_str_for_vaspkit=elements_str_for_vaspkit,
                                from_contcar=from_contcar,load=load,
                                load_CHGnet = load_CHGnet,load_path=load_path,load_model=load_model,
                                vac_dope=vac_dope,vac_as = vac_as,
                                open_diffusion=open_diffusion,
                                diffusion_specie=diffusion_specie)

        print('三','\n')
        print(step_object.current_structure_state.vasp_folder_path)
        for _ in range(num_loops):
            start_time=time.time() 
            if not load_model:
                step_object.current_structure_state.get_energy()
                E_1 = step_object.current_structure_state.energy
                print(f'--------------{E_1}--------------')
                # 计算 E_2: 利用 集群(vasp) 或 机器学习模型(ML model)
                next_structure_vasp_folder = step_object.next_structure_state.vasp_folder_path
                vasp_task = VaspTask(next_structure_vasp_folder)
                vasp_task.generate_input_files(gen_poscar=False,vac_dope=vac_dope)
                vasp_task.mpirun(pbs_nodefile=self.pbs_nodefile,
                                np=self.np,
                                dxec=self.dxec)
                #计算 E_2
                step_object.next_structure_state.get_energy()
                E_2 = step_object.next_structure_state.energy


            if load_model =='chgnet':
                if not os.path.exists(os.path.join(step_object.current_structure_state.vasp_folder_path,'relaxation_output.txt')):
                    E_1 = step_object.current_structure_state.get_CHG_energy()
                E_1 = step_object.current_structure_state.get_already_predict_energy()
                if not os.path.exists(os.path.join(step_object.next_structure_state.vasp_folder_path,'relaxation_output.txt')):
                    E_2 = step_object.next_structure_state.get_CHG_energy()
                E_2 = step_object.next_structure_state.get_already_predict_energy()

            if load_model =='mattersim':
                if not os.path.exists(os.path.join(step_object.current_structure_state.vasp_folder_path,'relaxation_output.txt')):
                    get_E1 = mattersim_predict.load(step_object.current_structure_state.poscar_path,
                                                    load_path=load_path)
                    E_1 = get_E1.energy
                E_1 = step_object.current_structure_state.get_already_predict_energy()
                if not os.path.exists(os.path.join(step_object.next_structure_state.vasp_folder_path,'relaxation_output.txt')):
                    get_E2 = mattersim_predict.load(step_object.next_structure_state.poscar_path,
                                                    load_path=load_path)
                    E_2 = get_E2.energy
                E_2 = step_object.next_structure_state.get_already_predict_energy()
                

            end_time = time.time()
            execution_time=end_time - start_time

            #书写交换概率
            exchange_mark,possibility = Exchange.mark(E_1=float(E_1), E_2=float(E_2), T=T)
            accept_path = os.path.join(step_object.next_structure_state.vasp_folder_path,'Accept.txt')
            with open(accept_path, 'a') as f:
                f.write(f'本次搜索继承概率为：{possibility:.6f}\n')
            print(f'进入循环，第{_+1}次')

            if time_save:
                time_path = os.path.join(step_object.next_structure_state.vasp_folder_path,'time_record.txt')       #储存处理
                with open(time_path, 'a') as f:
                    f.write(f'{execution_time:.6f}\n')
            
            if exchange_mark:
                step_object.walk()
            else:
                step_object.walk_anew()
