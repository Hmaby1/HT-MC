# pyMC 使用说明文档

**Monte Carlo 方法在多组分材料体系中的应用（Python 包）**

---

## 1. 安装方法（Installation）

### 依赖环境（Prerequisites）

- Python >= 3.10  
- gcc (GCC) = 10.2.0  

>  **提示**：Python 和 GCC 的版本要求是为了防止在低版本下安装 MLPs 软件包出错。

### 安装步骤

```bash
# 进入 setup.py 所在目录
$ cd pyMC

# 使用可编辑模式安装
$ pip install -e .
```
## 2. 使用方法（Usage）
支持的蒙特卡洛算法：


1.DFT-MC：基于 VASP 的第一性原理计算

2.CHGNet-MC：基于 CHGNet 机器学习势能

3.MatterSim-MC：基于 MatterSim 机器学习势能

 适用于结构弛豫过程，并支持包含空位（Vacancy）的材料结构模拟。

### 2.1初始 POSCAR 文件准备
#### 步骤 1：创建目录
```bash
$ mkdir -p MC_file/0
$ cd MC_file
```
#### 步骤 2：准备 POSCAR 文件

##### 【DFT-MC】
需要准备一个已经完成结构优化的目录，至少包含以下文件：


INCAR, POSCAR, POTCAR, KPOINTS, CONTCAR, vasprun.xml
```bash
$ cp -r <your relaxed file path> 0/
```
##### 【CHGNet-MC / MatterSim-MC】
 仅需准备初始 POSCAR 文件：
 ```bash
$ cp <your POSCAR path> 0/POSCAR
```
##### 【如果结构包含空位】
##### 需在 `0` 文件夹中额外准备一个含有空位的 0.POSCAR 文件
```bash
$ cp <your vac_POSCAR path> 0/0.POSCAR
```

### 2.2. 编写 runMetropolis.py 主函数
#### 推荐的主函数格式如下：

```python
def run():
    metropolis_mc = Metropolis(
        pbs_nodefile=PBS_NODEFILE,
        np=NP,
        dxec=DXEC
    )

    metropolis_mc.run(
        poscar_path=poscar_path,
        num_loops=num_loops,
        T=T,
        sublattice_symbols_lst=sublattice_symbols_lst,
        from_contcar=from_contcar,
        elements_str_for_vaspkit=elements_str_for_vaspkit,
        load=load,
        load_model=load_model,
        load_path=load_path,
        vac_dope=vac_dope,
        vac_as=vac_as,
        time_save=True,
        open_diffusion=False,
        diffusion_specie=None,
        exchange_times=1
    )
```

#### 参数说明（Parameter Description）

Metropolis(...)：若为slurm作业系统则可设置为任意字符串
| 参数名            | 说明                     |
| -------------- | ---------------------- |
| `pbs_nodefile` | PBS 作业系统节点文件路径 |
| `np`           | 并行核数        |
| `dxec`         | 并行参数        |

run(...)：
| 参数名                        | 默认值                   | 说明                                                                                                           |
| -------------------------- | --------------------- | ------------------------------------------------------------------------------------------------------------ |
| `poscar_path`              | eg:"./MC\_file/0/POSCAR" | 初始 POSCAR 文件路径或自定义路径                                                                                         |
| `num_loops`                | eg: 1000                  | 总共执行的蒙特卡洛循环次数                                                                                                |
| `T`                        | eg:273.75                | 蒙卡模拟温度（单位：K）                                                                                                 |
| `from_contcar`             | True                  | 是否从上一步的 CONTCAR 读取结构                                                                                         |
| `sublattice_symbols_lst`   | eg:\[\["Sc,Sb"]]         | 参与交换的元素列表                                                                                                    |
| `elements_str_for_vaspkit` | eg:"Sc Sb Te"            | VASP 计算中 POSCAR 结构所用到的元素顺序字符串                                                                                |
| `load`                     | False                 | 是否启用读取模式（True 表示从某一结构继续搜索）                                                                                   |
| `load_model`               | None                  | 是否使用机器学习势能（支持 "chgnet" 或 "mattersim"）                                                                        |
| `load_path`                | None                  | MLP 模型的训练权重路径：<br>• chgnet: `None`<br>• mattersim: `'MatterSim-v1.0.0-1M.pth'` 或 `'MatterSim-v1.0.0-5M.pth'` |
| `vac_dope`                 | False                 | 若结构含有空位，设置为 True                                                                                             |
| `vac_as`                   | "V"                   | 在 0.POSCAR 中代表空位的原子符号                                                                                        |
| `time_save`                | True                  | 是否输出每一步的耗时                                                                                                   |
| `open_diffusion`           | False                 | 是否开启 5 Å 范围内的扩散模拟                                                                                            |
| `diffusion_specie`         | None                  | 若指定，表示该原子每次都参与交换过程                                                                                           |
| `exchange_times`           | 1                     | 每次交换的原子数对                                                                                                    |


> 备注：

> PBS 作业系统参数在 SLURM 系统下可以设置为任意字符串；

> 若模型中包含 vacancy，务必准备好 0/0.POSCAR 并开启 vac_dope=True；

> 当前支持的机器学习模型包括：CHGNet 和 MatterSim。


# 3.输出文件（output）
```bash
$ls
0 1 2 3 4 5 6 7 8 ... steps.log [procss] runMetropolis.py runMetropolis.sh
#每个文件夹中均含一个time_record.txt记录计算时间(s)，Accept.txt记录接收概率
#每个交换接收步数中会额外包含exchanged.txt标志文件
#若打开Vac_dope，每个文件中会包含n.vasp空位文件，交换步数会额外包含n.CONTCAR空位文件
#，目录中会额外出现process文件夹包含所有空位文件
```

#### DFT-MC
标准的vasp relaxation计算文件,
#### CHGNet-MC / MatterSim-MC
仅存在POSCAR CONTCAR relaxation_output.txt计算输出文件


# 4.运行脚本参考（example）
## runMetropolis.py
```python

import sys
from MCobjects.Metropolis.main import Metropolis
 
PBS_NODEFILE = '1'#sys.argv[1]
NP = '1'#sys.argv[2]
DXEC = '1'#sys.argv[3]


#1.文件路径设置
poscar_path = r"D:\Desk\mattersim5M\SS_mattersim\0\POSCAR"

#2.筛选温度设置
num_loops = 4000
T = 444

#3.MC-搜索位点
sublattice_symbols_lst = [ 
                        ["Sc","Sb", "V"],
                             ]
elements_str_for_vaspkit = "Sc Sb Te"

#4.模型加载
load_model = 'mattersim' #/'chgnet'
load_path='MatterSim-v1.0.0-5M.pth'

#5.空位缺陷开关和其余设置
vac_dope =True
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
```
