import sys
from MCobjects.Metropolis.main import Metropolis
 
PBS_NODEFILE = '1'#sys.argv[1]
NP = '1'#sys.argv[2]
DXEC = '1'#sys.argv[3]


#1.文件路径设置
poscar_path = r"D:\Desk\新建文件夹 (3)\0\POSCAR"

#2.筛选温度设置
num_loops = 4000
T = 444

#3.MC-搜索位点
sublattice_symbols_lst = [ 
                        ["Sc","Sb"],
                             ]
elements_str_for_vaspkit = "Sc Sb Te"

#4.模型加载
load_model = 'mattersim' #/'chgnet'
load_path='MatterSim-v1.0.0-5M.pth'

#5.空位缺陷开关和其余设置
vac_dope =False
vac_as = "V"
from_contcar = True
load = False

#6.程序加载与启动
def run():
    metropolis_mc = Metropolis(pbs_nodefile=PBS_NODEFILE,
                                np=NP,
                                dxec=DXEC)
    
    metropolis_mc.run(
                    poscar_path=poscar_path,
                    num_loops=num_loops,
                    T=T,
                    sublattice_symbols_lst=sublattice_symbols_lst,
                    from_contcar=from_contcar,
                    elements_str_for_vaspkit=elements_str_for_vaspkit,
                    load=load,
                    load_model=load_model,
                    load_path = load_path,
                    vac_dope = vac_dope,
                    vac_as=vac_as,
                    time_save=True,
                    open_diffusion=False,
                    diffusion_specie=None,
                    exchange_times=1
                    ) 

if __name__ == "__main__":
    run()
