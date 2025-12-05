import os
import json
import sys
from pymatgen.core import Structure
from prettytable import PrettyTable
from pymatgen.io.vasp import Oszicar

from generateNewStructure.exchangeAtoms import ExchangeAtoms
from model.mattersim_ import mattersim_predict

from chgnet.model import CHGNet
from chgnet.model import StructOptimizer
from io import StringIO

class StructureState(object):
    """用于读取结构能量信息以及文件路径信息等

    Attributes
    ----------
    poscar_path : str
        class输入的POSCAR路径
    vasp_folder_path : str
        根据POSCAR读取其所在的VASP路径(POSCAR上级目录)
    vasp_folders_path : str
        根据POSCAR读取其所在的搜索路径(POSCAR上上级目录)
    structure_index : int
        根据VASP文件路径读取结构编号
    contcar_path : str
        根据POSCAR读取同目录下的CONTCAR路径
    vac_dope : bool
        是否属于空位掺杂结构
    load_CHGnet : bool or str
        是否加载模型。为str时将加载str路径模型
    CHG_out_path : str
        设置模型输出文件路径
    """
    def __init__(self,poscar_path: str,vac_dope = False,load_CHGnet = False,load_path=None):

        self.poscar_path = poscar_path
        self.vasp_folder_path = os.path.dirname(self.poscar_path)
        self.vasp_folders_path = os.path.dirname(self.vasp_folder_path)

        self.contcar_path = os.path.join(self.vasp_folder_path, "CONTCAR")
        
        self.vac_dope = vac_dope

        self.structure_index = int(os.path.split(self.vasp_folder_path)[-1])
        self.load_CHGnet=load_CHGnet
        self.CHG_out_path = os.path.join(self.vasp_folder_path,'relaxation_output.txt')

    def load_model(self,load_CHGnet,load_path):
        #加载初始化模型
        if load_path is not None:
            self.model = CHGNet.from_file(load_path)
            print('load trained model')
        else:
            self.model = CHGNet.load()
            print('load initial model')
        #初始化文件计算
        if load_CHGnet and (not os.path.exists(self.CHG_out_path)):
            self.get_CHG_energy()
        return


    def __repr__(self):
        table = PrettyTable(["Structure_Index", "POSCAR path"])
        table.add_row([self.structure_index, self.poscar_path])
        print(table)
        return ''
    
    def __str__(self):
        return self.__repr__()
    
    def get_energy(self):
        oszicar_path = os.path.join(self.vasp_folder_path, "OSZICAR")
        energy = float(Oszicar(oszicar_path).final_energy)
        
        self.energy = energy

        return energy
    
    def get_already_predict_energy(self):
        with open(self.CHG_out_path, "r") as f:
            content = f.readlines()
            line=content[-1].split(":")
        energy_CHG = line[-1]
        return energy_CHG

    def get_CHG_energy(self):

        structure = Structure.from_file(self.poscar_path)
        relaxer = StructOptimizer()
        structure.perturb(0.1)

        original_stdout = sys.stdout
        output_buffer = StringIO()
        sys.stdout = output_buffer
        result = relaxer.relax(structure, verbose=True) #分子弛豫优化，默认step = 500
        sys.stdout = original_stdout

        with open(self.CHG_out_path, "w") as f: #保存优化过程
            f.write(output_buffer.getvalue())  
            f.close()

        relaxed_structure=result["final_structure"]
        relaxed_structure.to(os.path.join(self.vasp_folder_path,"CONTCAR"),'poscar')
        #保存优化结构
        predic = self.model.predict_structure(relaxed_structure)
        energy_CHG = predic['e']

        with open(self.CHG_out_path, "a") as f:
            f.write(f'\nthe final structure energy:{energy_CHG*relaxed_structure.num_sites}')
       
        energy_CHG = self.get_already_predict_energy()
        return energy_CHG

    def get_mattersim_energy(self):

        energy_sim = self.get_already_predict_energy()
        return energy_sim

class StepObject(object):
    """用于执行结构交换的生成过程，并记录交换后的路径等信息。

    Attributes
    ----------
    current_structure_state : str
        class输入的POSCAR路径
    sublattice_symbols_lst : list
        初始设置交换晶格位点的元素列表
    elements_str_for_vaspkit : str
        vasp计算的元素与顺序
    from_contcar : bool
        新结构是否来自于上一步CONTCAR
    load : bool
        是否是读取模式
    load_CHGnet : bool
        是否加载CHGNet模型
    load_path : str
        训练模型路径(若有)
    vac_dope : bool
        是否属于空位掺杂结构
    vac_as : str
        作为空位的元素符号
    open_diffusion : bool
        是否以5埃为半径模拟短程扩散
    diffusion_specie : str
        必须参与交换的元素
    exchange_times : int
        每次结构的交换原子对数

    Note
    ----
    默认为非读取模式:load=False,搜索从文件夹0开始;读取模式下可自定义初始文件夹及序号

    进行一次完整的交换原子流程,生成下一步文件夹以及POSCAR/.poscar
    """
    def __init__(self, current_structure_state:StructureState, sublattice_symbols_lst: list,
                elements_str_for_vaspkit:str,
                from_contcar: bool,load: bool = True,
                load_CHGnet = False,load_path=None,load_model:str=None,
                vac_dope = False,vac_as = 'V',
                open_diffusion:bool=False,
                diffusion_specie:str=None,
                exchange_times:int=1
                ):
        self.sublattice_symbols_lst = sublattice_symbols_lst
        self.from_contcar = from_contcar
        self.current_structure_state = current_structure_state
        self.elements_str_for_vaspkit = elements_str_for_vaspkit
        
        self.load_CHGnet = load_CHGnet
        self.load_path = load_path

        self.vac_dope = vac_dope
        self.vac_as = vac_as

        self.open_diffusion = open_diffusion
        self.exchange_times = exchange_times
        # 日志文件的路径
        vasp_folders_path = self.current_structure_state.vasp_folders_path
        self.step_log_path = os.path.join(vasp_folders_path, "steps.log")  #总目录下生成步数文件
        if not os.path.exists(self.step_log_path):
            os.system("touch {0}".format(self.step_log_path))

        if not load:
            self.exchange_steps = 0
            self.total_steps = 0
        else:
            self.total_steps, self.exchange_steps = self.load_info()
        self.diffusion_specie = diffusion_specie
        self.next_structure_state = self._get_next_state()
        self.save_info()
    
    def _get_next_state(self):
        if self.from_contcar:
            exchanger = ExchangeAtoms(poscar_path=self.current_structure_state.contcar_path,            #current_structure_state发生变化
                                sublattices_symbols_lst=self.sublattice_symbols_lst,
                                vac_dope=self.vac_dope,
                                vac_as=self.vac_as)
        else:
            exchanger = ExchangeAtoms(poscar_path=self.current_structure_state.poscar_path,
                                sublattices_symbols_lst=self.sublattice_symbols_lst,
                                vac_dope=self.vac_dope,
                                vac_as=self.vac_as)

        next_poscar_path = exchanger.generate_new_structure(new_structure_index=self.total_steps+1,
                                                        elements_str_for_vaspkit=self.elements_str_for_vaspkit,
                                                        pick_first_specie=self.diffusion_specie,
                                                        with_cutoff=self.open_diffusion,
                                                        exchange_times=self.exchange_times
                                                        )



        next_structure_state = StructureState(poscar_path=next_poscar_path,
                                         vac_dope=self.vac_dope,
                                         load_CHGnet=self.load_CHGnet,
                                         load_path=self.load_path)
        
        if self.load_CHGnet:
            next_structure_state.load_model(load_CHGnet = self.load_CHGnet,
                                            load_path=self.load_path)

        return next_structure_state
    
    def __repr__(self):
        print("Current StructureState:")
        print(self.current_structure_state)
        print("Next StructureState:")
        print(self.next_structure_state)
        return ''

    def __str__(self):
        return self.__repr__()

    def save_info(self):
        '''
        Description
        -----------
            1. 每一步结束，无论是否交换原子，都需要执行这个函数 -- 保存 step 信息
        '''
        dict_steps = {"Total_steps": self.total_steps, "Exchanged_steps": self.exchange_steps}
        with open(self.step_log_path, "w") as f:
            json.dump(dict_steps, f)

    def load_info(self):
        '''
        Description
        -----------
            1. 每一步开始，都需要 load 信息
            2. `self.__init__()`, `self.walk()`, `self.walk_anew()` 自动调用 `self.load_info()`
        
        Return
        ------

        '''
        try:
            with open(self.step_log_path, "r") as f:
                dict_steps = json.load(f)
            return dict_steps["Total_steps"], dict_steps["Exchanged_steps"]

        except:
            return None

    def walk(self):
        '''
        Description
        -----------
            1. 满足条件， 
                - 在 current_structure 文件夹下，建立一个 `exchanged.txt` 文件
                - self.current_structure_state = self.next_structure_state
                - self.next_structure_satte = 新的 structure_state
        '''
        # 交换的步数加 1
        self.total_steps += 1
        self.exchange_steps += 1
        
        exchanged_txt_path = os.path.join(self.current_structure_state.vasp_folder_path, "exchanged.txt")
        os.system("touch {0}".format(exchanged_txt_path))

        self.current_structure_state = self.next_structure_state    #数据进行迭代
        self.next_structure_state = self._get_next_state()
        
        self.save_info()
    

    def walk_anew(self):
        '''
        Descroption
        -----------
            1. 不满足条件 (原子不发生交换)
                - self
        '''
        self.total_steps += 1
        self.next_structure_state = self._get_next_state()
        
        self.save_info()