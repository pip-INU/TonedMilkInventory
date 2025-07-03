import numpy as np
from utilities import *
from data import ProblemParameters
from inventorytime import CommonTime
class RawInventory:
    def __init__(self, T1,T4, inventory, params: ProblemParameters,ct: CommonTime):
        self.params = params
        self.T1 = T1
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T4 = T4
    def dI_r0(self, t):
        return self.params.mu0 * t - self.params.P0*safe_exp(-self.params.delta_2*t)
    def I_r0(self, t):
        return (self.params.mu0 / 2) * (t**2 - self.T1**2) + \
               (self.params.P0 / self.params.delta_2) * (np.exp(-self.params.delta_2 * t) - np.exp(-self.params.delta_2 * self.T1))
    def integral_I_r0(self, a, b):
        term1 = (self.params.mu0 / 6) * (b**3 - a**3)
        term2 = -(self.params.mu0 * self.T1**2 / 2) * (b - a)
        term3 = (self.params.P0 / (self.params.delta_2**2)) * (
            safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b))
        term4 = (self.params.P0 / self.params.delta_2) * safe_exp(-self.params.delta_2 * self.T1) * (a - b)
        return term1 + term2 + term3 + term4
    def dI_r1(self, t):
        return - self.params.P0*safe_exp(-self.params.delta_2*t)
    def I_r1(self, t):
        return (self.params.P0 / self.params.delta_2) * (np.exp(-self.params.delta_2 * t) - np.exp(-self.params.delta_2 * self.T4))
    def integral_I_r1(self, a, b):
        term1 = (self.params.P0 / (self.params.delta_2**2)) * (
            safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b))
        term2 = (self.params.P0 / self.params.delta_2) * safe_exp(-self.params.delta_2 * self.T4) * (a - b)
        return term1 + term2
    def total_raw(self):
        return self.integral_I_r0(0,self.T3)+self.integral_I_r1(self.T3,self.T4)
    
class ReprocessableInventory:
    def __init__(self, T1,T4, inventory, params: ProblemParameters,ct: CommonTime):
        self.params = params
        self.T1 = T1
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T5 = ct.solve_T5(T1, T4)
        self.T4 = T4
    def dI_rw0(self, t):
        return self.params.H*self.params.P0*safe_exp(-self.params.delta_2*t)
    def I_rw0(self, t, T1):
        return (self.params.H * self.params.P0 / self.params.delta_2) * (np.exp(-self.params.delta_2 * T1) - np.exp(-self.params.delta_2 * t))
    def integral_I_rw0(self, a, b):
        A = (self.params.H * self.params.P0) / self.params.delta_2
        term1 = A * safe_exp(-self.params.delta_2 * self.T1) * (b - a)
        term2 = (A / self.params.delta_2) * (safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b))
        return term1 + term2
    def dI_rw1(self, t):
        return self.params.eeta*self.params.R0*safe_exp(-self.params.delta_3*t)
    def I_rw1(self, t):
        return (self.params.R0 / self.params.delta_3) * (np.exp(-self.params.delta_3 * t) - np.exp(-self.params.delta_3 * self.T5))
    def integral_I_rw1(self, a, b):
        term1 = (self.params.R0 / (self.params.delta_3**2)) * (
            safe_exp(-self.params.delta_3 * a) - safe_exp(-self.params.delta_3 * b))
        term2 = (self.params.R0 / self.params.delta_3) * safe_exp(-self.params.delta_3 * self.T5) * (a - b)
        return term1 + term2
    def total_reprocess(self):
        return self.integral_I_rw0(self.T1,self.T4) + self.integral_I_rw1(self.T4,self.T5)