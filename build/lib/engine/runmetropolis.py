'''
Descripttion: 
version: 
Author: sch
Date: 2022-03-29 20:30:48
LastEditors: sch
LastEditTime: 2022-03-29 21:53:24
'''
import sys

from pyMC.MCobjects.Metropolis.main import Metropolis


PBS_NODEFILE = sys.argv[1]
NP = sys.argv[2]
DXEC = sys.argv[3]

poscar_path = "/Users/mac/我的文件/Mycode/new/new/pymc/data/mc_1/0/POSCAR"
num_loops = 1000
T = 923

sublattice_symbols_lst = [ ["Re", "Nb"],
                            ["S", "Se"] ]

from_contcar = True
elements_str_for_vaspkit = "Re Ta S Se"
load = False


def run():
    metropolis_mc = Metropolis(pbs_nodefile=PBS_NODEFILE,
                                np=NP,
                                dxec=DXEC)
    
    metropolis_mc.run(poscar_path=poscar_path,
                    num_loops=num_loops,
                    T=T,
                    sublattice_symbols_lst=sublattice_symbols_lst,
                    from_contcar=from_contcar,
                    elements_str_for_vaspkit=elements_str_for_vaspkit,
                    load=load) 


if __name__ == "__main__":
    run()