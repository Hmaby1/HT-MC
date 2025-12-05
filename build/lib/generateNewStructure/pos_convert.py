from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.util.typing import SpeciesLike
from typing import Sequence
import os

new_poscar_path='D:\\桌面\\pyMC\\POSCAR_write\\write\\'
new_poscar_name = 'POSCAR'
class poscar_convert:
    def __init__(self,doc_path:str,doc_type:str = 'poscar',vac_as: Sequence[SpeciesLike] = ('V') ):
        self.path = doc_path
        self.type = doc_type
        self.new_pos = os.path.dirname(self.path)

        self.vac_as = vac_as
    def __repr__(self):
        print('')
        return self
    


    def cif_trans(self):
        structure = Structure.from_file(self.path)
        poscar = Poscar(structure)
        pos_path = os.path.dirname(self.path)
        poscar.write_file(pos_path+'\\1')
        self.poscar_trans(pos_path+'\\1')
        #os.remove(pos_path+'1')
        return pos_path

    def poscar_trans(self):
        if self.type.lower() == 'cif':
            path = os.path.dirname(self.path) + '\\1'
        elif self.type.lower() == 'poscar':
            path = self.path
        elif self.type.lower() == 'pos':
            path = self.path

        structure = Structure.from_file(path)
        structure.remove_species(self.vac_as)
        Structure.to(structure,path,"poscar")
        return 
        # with open(path,'r',encoding="utf-8") as f:
        #     lines = f.readlines()
        #     num_lines = len(lines)
        #     #取出元素和数量
        #     ele_line = lines[5]
        #     num_ele_line = lines[6]
        #     el_pre = ele_line.strip().split(" ")
        #     num_pre = num_ele_line.strip().split(" ")
        #     el=[i.strip() for i in el_pre if i.strip()]
        #     num=[i.strip() for i in num_pre if i.strip()]  
            
        #     #处理元素序列
        #     el_index_lis=[]
        #     q = 0
        #     for i in num:
        #         t=int(i)
        #         ele = el[q]
        #         q += 1
        #         for j in range(t):
        #             el_index =  ele+str(j+1)
        #             el_index_lis.append(el_index)                                 #生成元素标签
        #     #print(el),print(num),print(el_index_lis)
        #     try:
        #         int(lines[9].rstrip()[-2:].strip())
        #     except ValueError:
        #         print('poscar末尾已有元素标记，将进行删除处理')                #需要优化，该判断仅针对cifPOSCAR
        #         for i in range(8,len(lines)):
        #             lines[i]=lines[i][:-3].rstrip()
        #     #元素行判断是否有空位并删除
        #     Vac_remove = ['V','Va', 'Vac','']                                #V为空位cif使用pymatgen生成的POSCAR并非元素V，当然该列表也可以用来删除元素
        #     for i in Vac_remove:
        #         for j in range(len(el)):
        #             if i == el[j]:
        #                 Vac_index = j
        #                 self.Vac_num = num[j]
        #                 ele_line = ele_line.replace(i,'')
        #                 lines[5] = ele_line
        #                 num_ele_line = num_ele_line.replace(self.Vac_num,'')
        #                 lines[6] = num_ele_line
        #     content=list(lines[:8])

        #     for i in range(8,num_lines):
        #         content.append(lines[i].rstrip()+' '+el_index_lis[i-8]+'\n')      #坐标末添加元素序列
        #         if i == len(el_index_lis)+7:       #
        #             break


        #     lines_to_remove=[]
        #     startpos = 8
        #     if Vac_index != 0:
        #         k = 0
        #         for i in range(Vac_index):
        #             k += int(num[i])
        #         startpos += k
        #     endpos = startpos + int(self.Vac_num)
        #     lines_to_remove = list(range(startpos, endpos))                #找出空位位置行

        #     for i in lines_to_remove:
        #         self.Vac_lines.append(content[i])                                  #储存空位行

        #     for i in sorted(lines_to_remove, reverse=True):                    #倒装删除空位坐标
        #         content.pop(i)

        #     with open(self.path,'w+') as copy:          #生成新的标注POSCAR
        #         copy.writelines(content)
        #     copy.close()
        #     int_numbers = [int(num) for num in num]
        #     totle = sum(int_numbers)-int(self.Vac_num)
        # f.close()

    

    def poscar_verse_trans(self,v_pos_path,con_path,new_path):
        '''
        将v_pos_path空位提取\n
        将con_path内容提取\n
        结合储存为new_path\n
        '''
        v_pos_path = v_pos_path                             #包含空位文件
        con_path = con_path                                 #不包含空位文件
        new_path=new_path                                   #输出文件

        structure = Structure.from_file(v_pos_path)
        con_stru = Structure.from_file(con_path)
        element_counts = structure.composition.as_dict()
        eles = []
        id = []
        total = 0
        for element, count in element_counts.items():
            count = int(count)
            total += count
            eles.extend([element] * int(count))
            id.extend(list(range(total - count + 1, total + 1)))
        atom_index = dict(zip(id, eles))

        v_lis=[]
        for k,v in atom_index.items():
            if v == self.vac_as:
                v_lis.append(k)
        coordinate = structure.frac_coords
        v_coordinate = []
        for i in v_lis:
            v_coordinate.append(coordinate[int(i)-1])
        for coor in v_coordinate:
            con_stru.append(self.vac_as,coor)
        Structure.to(con_stru,new_path,'poscar')