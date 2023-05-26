# %%
# Data


import numpy as np
from quaternion import *

N = 16


q0 = quaternion(0,0,0,0)
q1 = quaternion(1,0,0,0)
qi = quaternion(0,1,0,0)
qj = quaternion(0,0,1,0)
qk = quaternion(0,0,0,1)
qq = quaternion(0.5,0.5,0.5,0.5)

map_quaternion_to_string = {
    q1 : "+",
    -q1: "-",
    qi : "i",
    -qi: "I",
    qj : "o",
    -qj: "O",
    qk : "u",
    -qk: "U",
    qq : "q",
    qq*-1: "Q",
    qq*qi : "ì",
    qq*-qi: "Ì",
    qq*qj : "ò",
    qq*-qj: "Ò",
    qq*qk : "ù",
    qq*-qk: "Ù",
}

Qplus = [i.copy() for i in map_quaternion_to_string.keys()]


# %%
# Class definition


class QS:

    def __init__(self, size):
        self.values = [q1.copy() for i in range(size)]
        self.size = size

    def set_values(self, values):
        self.values = values

    def set_value(self, value, index):
        self.values[index] = value



    def periodic_autocorrelation(self,t):

        sum_res = q0.copy()
        for i in range(self.size):
            sum_res += self.values[i]*(self.values[(i+t)%self.size]).conjugate()

        return sum_res

    def odd_periodic_autocorrelation(self,t):

        sum_res = q0.copy()
        for i in range(self.size):
            power = (-1)**((i+t)//self.size)
            sum_res += power*self.values[i]*(self.values[(i+t)%self.size]).conjugate()

        return sum_res



    def is_pqs(self):

        for t in range(1,self.size):
            if self.periodic_autocorrelation(t) != q0:
                return False

        return True

    def is_opqs(self):

        for t in range(1,self.size):
            if self.odd_periodic_autocorrelation(t) != q0:
                return False

        return True


    def sub_opqs(self):

        for start in range(1,(self.size-1)//2):
            size = self.size - 2*start
            pqs = QS(size)
            pqs.set_values(self.values[start:-start])
            if pqs.is_opqs():
                return pqs

        return None



    def __str__(self):
        res_string = "["
        for q in self.values:
            res_string += map_quaternion_to_string[q]

        return res_string + "]"





# %%
# Naive algorithm

count = 0

def find(size):
    global count
    count = 0
    pqs = QS(size)
    pqs.set_values([q1 for _ in range(size)])

    find_recursive(pqs, size, 1)

    print(f"total number: {count}")


def find_recursive(pqs, size, index):

    if index >= (size+1)//2:
        if pqs.is_opqs():
            global count
            count+=1
            print(pqs)
        return

    for index_value_to_test in range(N):
        pqs.set_value(Qplus[index_value_to_test], index)
        pqs.set_value(Qplus[index_value_to_test], size-1-index)

        find_recursive(pqs, size, index+1)













