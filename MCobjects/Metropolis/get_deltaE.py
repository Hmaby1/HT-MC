from chgnet.model import CHGNet
from pymatgen.io.vasp import Vasprun
from pymatgen.core import Structure
import os,json
import matplotlib.pyplot as plt

#得到chgnet预测与vasp计算结果差值

def both(vasp_folder_path:str,model_path:str):
    chg=CHGNet.load()
    model=CHGNet.from_file(model_path)
    vasp_folder_lis = []
    for root, dirs, files in os.walk(vasp_folder_path):
        if 'vasprun.xml' in files:
            vasp_folder_lis.append(root)

    energy_vasp=[]
    energy_chg=[]
    energy_model=[]
    for folder in vasp_folder_lis:
        vasprun = Vasprun(os.path.join(folder,'vasprun.xml'))
        energy_vasp.append(float(vasprun.final_energy))
        structure = vasprun.final_structure
        stru_num = len(structure.atomic_numbers)
        energy_model.append((float(model.predict_structure(structure)['e']))*int(stru_num))
        energy_chg.append((float(chg.predict_structure(structure)['e']))*int(stru_num))
    return energy_vasp,energy_model,energy_chg

def plot_differ(energy_vasp,energy_model,energy_chg):
    y1,y2,y3 = energy_vasp,energy_model,energy_chg
    delta1=[a-b for a,b in zip(y2,y1)]
    delta2=[a-b for a,b in zip(y3,y1)]
    x = range(len(y1))

    ### 0. 全局设置
    # 0.1. 字体设置
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Times New Roman'
    plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
    plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'
    # 0.2. 刻度线朝内
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.scatter(x, delta1,      #散点图
            s=7,
            color="#FD6D5A",
            edgecolors="black",
            label='model')
    plt.scatter(x, delta2,      #散点图
            s=7,
            color="#FD6D5A",
            edgecolors="black",
            label='chg')
    plt.xlabel("Steps", fontsize=28, fontweight="bold")
    plt.ylabel("Energy (eV)", fontsize=28, fontweight="bold")

    # 3. xaxis, yaxis
    plt.xticks(fontsize=20, fontweight="bold")
    plt.yticks(fontsize=20, fontweight="bold")

    # 4. y range()
    #plt.ylim()

    # 5. 刻度线粗细
    plt.tick_params(
        width=2,        # 刻度线的粗细
        length=5,       # 刻度线的长短
        #labelsize=28   # 刻度线的字体大小
    )

    # 5. 设置坐标轴的粗细
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(1.5);###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(1.5);####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(1.5);###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(1.5);###设置右边坐标轴的粗细

    # 6. label字体
    legend_font = {"size" : 23, 
            "weight": "bold"
            }
    plt.legend(loc=(0.6, 0.75),
            prop=legend_font,
            frameon=False)

    plt.show()
    return
