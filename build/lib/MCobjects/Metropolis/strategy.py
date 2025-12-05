'''
Descripttion: 
version: 
Author: sch
Date: 2022-03-29 19:06:20
LastEditors: sch
LastEditTime: 2022-03-31 08:58:54
'''
import math
import numpy as np 


class Exchange(object):
    # Boltzmann constant
    k = 8.617333262145E-5

    @classmethod
    def mark(cls, E_1:float, E_2:float, T:float) -> float:
        random_number = np.random.rand()
        
        delta_E = E_2 - E_1

        if ( E_1 > E_2):
            possibility=1
            return True,possibility        

        possibility = math.exp( - delta_E / (cls.k * T) )

        if ( possibility > random_number ):
            return True,possibility
        
        return False,possibility