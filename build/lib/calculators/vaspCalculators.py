'''
Descripttion: 
version: 
Author: sch
Date: 2022-01-24 16:29:21
LastEditors: sch
LastEditTime: 2022-03-30 09:26:55
'''
import os
import time
import logging
from pymatgen.io.vasp.outputs import Oszicar

from utilitys.mpirunContext import PwdContext
from logger.loggerForVaspTask import LoggerForVaspTask


class VaspTask:
    '''
    Attributes
    ----------
        1. self.vasp_folder_path: str,
            The path of the vasp folder which includes INCAR, POSCAR... 
        2. self.vasp_folders_path: str,
            The path of folder including 0, 1, 2, ...
        3. self.pwd_path: str
            当前工作路径
    
    Log files
    ---------
        1. `process.log`:
            - 开始计算
            - 计算过程
            - 得到能量
        2. `error.log`:
            出错时建立
    '''
    TIME_INTERVAL = 30
    
    def __init__(self, vasp_folder_path: str):
        '''
        Parameters
        ----------
            1. vasp_folder_path: str,
                The path of the vasp folder which includes INCAR, POSCAR... 
        '''

        self.vasp_folder_path = vasp_folder_path
        self.vasp_folders_path = os.path.dirname(self.vasp_folder_path)
        self.pwd_path = os.getcwd()
        self.log_file_path = os.path.join(self.vasp_folder_path, "process.log")

        if not os.path.exists(self.vasp_folder_path):
            LoggerForVaspTask.log_output(log_file_path=self.log_file_path,
                                    msg="Vasp folder of structure {0} doesn't exist".format(
                                        os.path.split(self.vasp_folder_path)[-1]
                                    ),
                                    level=logging.ERROR)
            raise Exception("This structuew has been calculated!")


    def mpirun(self, pbs_nodefile: str, np: int, dxec: str):
        '''
        Description
        -----------
            Work as `mpirun -machinefile $PBS_NODEFILE -np $NP $DXEC > output` in pbs shell

        Parameters
        ----------
            1. pbs_nodefile: str,
                PBS_NODEFILE
            2. np: int,
                NP
            3. dxec: str,
                DXEC

        Note
        ----
            1. You must enter into this vasp folder, then you can `mpirun` a task
        '''
        # Enter into vasp folder and qsub, and the back to `pwd_path`
        LoggerForVaspTask.log_output(log_file_path=self.log_file_path,
                            msg="Calculating the energy of structure {0} ...".format( 
                            os.path.split(self.vasp_folder_path)[-1] ),
                            level=logging.INFO)
                
        with PwdContext(pwd_path=self.pwd_path, vasp_folder_path=self.vasp_folder_path) as _:
            os.system("mpirun -machinefile {0} -np {1} {2} > output".format(pbs_nodefile, np, dxec))
        

    def generate_input_files(self, gen_poscar:bool=False,vac_dope=False):
        '''
        Note
        ----
            1. You must enter into this vasp folder, then you can cp input files
        '''
        with PwdContext(pwd_path=self.pwd_path, vasp_folder_path=self.vasp_folder_path) as manager:

            # POSCAR
            if gen_poscar:
                self._generate_poscar()

            # POTCAR
            self._generate_potcar()
            # KPOINTS
            # INCAR
            # OPTCELL
            # runvasp.pbs
            for filename in ["KPOINTS", "INCAR", "OPTCELL"]:
                self._copy_input_file(filename=filename)
    
    def _generate_poscar(self):
        '''
        Monte-Carlo 交换得到新的 pymatgen.Structure  -->  POSCAR
        '''
        vasp_folder_0_path = os.path.join(self.vasp_folders_path, str(0))
        filename_0_path = os.path.join(vasp_folder_0_path, "CONTCAR")
        
        filename_path = os.path.join(self.vasp_folder_path, "POSCAR")
        os.system("cp -r {0} {1}".format(filename_0_path, filename_path))
        
        
    def _generate_potcar(self):
        os.system("vaspkit -task 103")

    def _copy_input_file(self, filename):
        vasp_folder_0_path = os.path.join(self.vasp_folders_path, str(0))
        filename_0_path = os.path.join(vasp_folder_0_path, filename)
        
        filename_path = os.path.join(self.vasp_folder_path, filename)
        os.system("cp -r {0} {1}".format(filename_0_path, filename_path))
    
    def wait_until_task_ends(self):
        time.sleep(self.TIME_INTERVAL)

        final_energy = self._get_energy()

        LoggerForVaspTask.log_output(log_file_path=self.log_file_path,
            msg="Energy of structure {0} = {1}, finding the next structure...".format( 
                                                    os.path.split(self.vasp_folder_path)[-1], final_energy),
            level=logging.INFO)

        return final_energy
    
    def _get_energy(self):
        oszicar_path = os.path.join(self.vasp_folder_path, "OSZICAR")
        final_energy = float(Oszicar(oszicar_path).final_energy)
        return final_energy

    '''
    def wait_until_task_ends(self):
        
        outcar_path = os.path.join(self.vasp_folder_path, "OUTCAR")
        completed_mark = False

        # 间隔 log_time_multi 后。记录一次
        with PwdContext(pwd_path=self.pwd_path, vasp_folder_path=self.vasp_folder_path):
            # 判断是否已经生成 OUTCAR
            while (True):
                if (os.path.exists(outcar_path)):
                    break
                else:
                    time.sleep(self.TIME_INTERVAL)

            # 判断 VASP 计算是否完成 (OUTCAR中是否包含 `Total CPU time used (sec):` 字段)
            while (True):
                with open(outcar_path, "r") as f:
                    if "Total CPU time used (sec):" in f.read():
                        completed_mark = True
                        break
                    else:
                        time.sleep(self.TIME_INTERVAL)

                LoggerForVaspTask.log_output(log_file_path=self.log_file_path,
                                        #msg="Can't find `Total CPU time used (sec):` in OUTCAR",
                                        msg="calculating...",
                                        level=logging.INFO)
            
            time.sleep(self.TIME_INTERVAL)

            if completed_mark:
                """
                You can custom this code at here to get anything you want.
                """
                final_energy = self._get_energy()

                LoggerForVaspTask.log_output(log_file_path=self.log_file_path,
                    msg="Energy of structure {0} = {1}, finding the next structure...".format( 
                                                            os.path.split(self.vasp_folder_path)[-1], final_energy),
                    level=logging.INFO)

        return final_energy
    '''