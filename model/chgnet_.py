from pymatgen.core import Structure
from pymatgen.io.vasp import Vasprun

from chgnet.model.model import CHGNet
import os

from ase import io
#from ..批量整理 import find_file

def get_vasp_folder(vasp_folder_path:str,vasp_folder_lis:list = []):
    '''
    从输入的路径中检索vasprun.xml输出其所在路径列表vasp_folder_lis
    '''
    for root, dirs, files in os.walk(vasp_folder_path):
        if 'vasprun.xml' in files:
            vasp_folder_lis.append(root)
    return vasp_folder_lis

def get_energy_all(vasp_path:str,energy_dict = {}):
    '''
    从vasprun.xml中提取所有步骤以及能量
    '''
    vasprun_path = os.path.join(vasp_path,'vasprun.xml')
    vasprun = Vasprun(vasp_path)
    vasprun.ionic_steps
    return energy_dict

def get_chg_energy(structure:Structure,model:CHGNet):
    prediction = model.predict_structure(structure)
    total_energy = prediction["e"]*structure.num_sites
    #print(f'------------------ structure E: {prediction["e"]*structure.num_sites} ------------------')
    return total_energy

def get_prediction_from_vasp(vasp_folder_path: str,from_contcar:bool=True,*,load:str=None,from_:bool=False,from_filename:str=None):
    '''
    从vasp_folder中读取结构使用模型进行预测能量
    Return:
    dir_id , energy_list
    '''
    dir_id=[]
    energy_list=[]

    if load:
        model = CHGNet.from_file(load)
    else:
        model = CHGNet.load()


    for root, dirs, files in os.walk(vasp_folder_path):
    #跳出下级目录循环
        if os.path.exists(os.path.join(root,'CONTCAR')):
            break
    #得到步数列表
        for dir_name in dirs:
            try:
                dir_int = int(dir_name)
                dir_id.append(dir_int)
            except ValueError:
                print(f"目录 {dir_name} 不是有效的整数,跳过")
        dir_id.sort()
    #预测
        for id in dir_id:
            if from_contcar:
                contcar_path = os.path.join(root,str(id),'CONTCAR')
            elif from_:
                contcar_path = os.path.join(root,str(id),from_filename)
            else:
                contcar_path = os.path.join(root,str(id),'POSCAR')

            print(f"----{id}----")
            if os.path.exists(contcar_path):
                contcar = Structure.from_file(contcar_path)
                prediction = model.predict_structure(contcar)
                total_energy = prediction["e"]*contcar.num_sites
                energy_list.append(total_energy)
            else:
                print(f'no {contcar_path}')

            #print(energy_list)
    return dir_id,energy_list

def find_file(directory:str,mark:str='exchanged.txt'):
    '''
    提取directory中有标志文件的文件夹
    ------
    return:\n
        markfolders:list[folder_name]
    
    for folder in directory:
        if mark in folder:
            markfolders.append(folder_name)
    '''
    markfolders=[]
    dirs = [i for i in os.listdir(directory) if os.path.isdir(os.path.join(directory,str(i)))]
    for dr in dirs:
        path=os.path.join(directory,dr)
        content = os.listdir(path)
        if mark in content:
            markfolders.append(dr)
    return markfolders
