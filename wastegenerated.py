from data import ProblemParameters
from utilities import *
from inventorytime import CommonTime
class UsefulWaste:
    def __init__(self, T1, T4,  params: ProblemParameters , ct: CommonTime ):
        self.params = params
        self.T5 = ct.solve_T5(T1,T4)
    def solids(self, T1, T4 ):
        return self.params.w_f*(self.params.P0*safe_exp(-self.params.delta_2*(T4-T1))+
                                self.params.R0*safe_exp(-self.params.delta_1*(self.T5-T4)))
    def cream(self, T1, T4 ):
        return self.params.w_z*(1-self.params.w_f)*(self.params.P0*safe_exp(-self.params.delta_2*(T4-T1))+
                                self.params.R0*safe_exp(-self.params.delta_1*(self.T5-T4)))
    def reusereturns(self):
        return np.sum(self.params.w_ri)*0.6

class UselessWaste:
    def __init__(self, T1, T4,  params: ProblemParameters, ct: CommonTime ):
        self.params = params
        self.T5 = ct.solve_T5(T1,T4)
    def spillage(self, T1, T4):
        return self.params.sigma*(self.params.P0*safe_exp(-self.params.delta_2*(T4-T1))+
                                self.params.R0*safe_exp(-self.params.delta_1*(self.T5-T4)))
    def unreuse(self,T1,T4):
        return np.sum(self.params.w_ri)*0.4