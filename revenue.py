from data import ProblemParameters
from utilities import *
from inventorytime import CommonTime
from wastegenerated import *
class CalculateRevenue:
    def __init__(self, T1, T4, demand, params:ProblemParameters, ct : CommonTime, uw:UsefulWaste):
        self.params = params
        self.demand = demand
        self.T1 = T1
        self.T4 = T4
        self.T5 = ct.solve_T5(self.T1,self.T4)
        self.solids = uw.solids(T1,T4)
        self.cream = uw.cream(T1,T4)
    def productrevenue(self):
        return self.params.sp * self.demand
    def solidsrevenue(self):
        return self.solids * self.params.r_f
    def creamrevenue(self):
        return self.cream * self.params.r_s
    def returnrevenue(self):
        return self.params.sp_re * self.params.omega * np.sum(self.params.w_ri)
    def total_revenue(self): 
        return self.productrevenue() + self.solidsrevenue() + self.creamrevenue() + self.returnrevenue()
                
