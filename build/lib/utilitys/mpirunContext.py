'''
Descripttion: 
version: 
Author: sch
Date: 2022-03-04 15:24:51
LastEditors: sch
LastEditTime: 2022-03-04 16:22:10
'''
import os


class PwdContext():
    '''
    Attributes
    ----------
        1. self.pwd_path: 当前工作路径
        2. self.vasp_folder_path: vasp folder路径（包含 INCAR, KPOINTS 等信息的路径）
    '''
    def __init__(self, pwd_path, vasp_folder_path):
        self.pwd_path = pwd_path
        self.vasp_folder_path = vasp_folder_path

    def __enter__(self):
        os.chdir(self.vasp_folder_path)
        return self
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        os.chdir(self.pwd_path)