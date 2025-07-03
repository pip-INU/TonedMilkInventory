import numpy as np
from utilities import *
from data import ProblemParameters
from inventorytime import CommonTime, BackorderTime
class ServiceableInventory:
    def __init__(self, T1, T4, params: ProblemParameters, ct: CommonTime, bt: BackorderTime, inventory = "backorder"):
        self.params = params
        self.T1 = T1
        self.T4 = T4
        self.T2=bt.solve_T2(T1)
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T5=ct.solve_T5(T1,T4)
    def dI_s0(self,t):
        return -self.params.D0*np.exp(self.params.delta_1 * t)
    def I_s0(self, t):
        return (self.params.D0 / self.params.delta_1) * (1 - np.exp(self.params.delta_1 * t))
    def integral_I_s0(self, a, b):
        term1 = (self.params.D0 / self.params.delta_1) * (b - a)
        term2 = (self.params.D0 / (self.params.delta_1**2)) * (safe_exp(self.params.delta_1 * a) - safe_exp(self.params.delta_1 * b))
        return term1 + term2
    def V(self,t):
        return self.params.K*self.params.P0*np.exp(-self.params.delta_2 * t)
    def integral_V(self, a, b):
        (self.params.K * self.params.P0 / self.params.delta_2) * (np.exp(-self.params.delta_2 * a) - np.exp(-self.params.delta_2 * b))
    def dI_s1(self,t):
        return self.params.V(t)-self.params.D0*np.exp(self.params.delta_1 * t)
    def I_s1(self, t, T1):
        return (((self.params.K * self.params.P0) / ( self.params.delta_2)) * (safe_exp(-self.params.delta_2 * self.T2) - safe_exp(-self.params.delta_2 * t)) + 
               (self.params.D0 / self.params.delta_1) * (safe_exp(self.params.delta_1 * self.T2) - safe_exp(self.params.delta_1 * t)))
    def integral_I_s1(self, a, b):
        A = (self.params.K * self.params.P0) / (self.params.delta_2)
        B = self.params.D0 / self.params.delta_1
        
        term1 = A * safe_exp(-self.params.delta_2 * self.T2) * (b - a)
        term2 = (A / self.params.delta_2) * (
            safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b)
        )
        term3 = B * safe_exp(self.params.delta_1 * self.T2) * (b - a)
        term4 = (B / self.params.delta_1) * (
            safe_exp(self.params.delta_1 * a) - safe_exp(self.params.delta_1 * b)
        )
        return term1 + term2 + term3 + term4
    def dI_s2(self,t):
        return self.params.V(t)-self.params.D0*np.exp(self.params.delta_1 * t)-self.params.theta*self.params.I_s2(t)
    def I_s2(self, t):
        return (((self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)) * 
               (np.exp(-(self.params.delta_2 - self.params.theta) * self.T2 - self.params.theta * t) - np.exp(-self.params.delta_2 * t)) + 
               (self.params.D0 / (self.params.delta_1 + self.params.theta)) * 
               (np.exp((self.params.delta_1 + self.params.theta) * self.T2 - self.params.theta * t) - np.exp(self.params.delta_1 * t)))
    def integral_I_s2(self, a, b):
        A = (self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)
        B = self.params.D0 / (self.params.delta_1 + self.params.theta)
        term1 = A * safe_exp(-(self.params.delta_2 - self.params.theta) * self.T2) * (-1/self.params.theta) * (safe_exp(-self.params.theta * b)
                                                                                                                - safe_exp(-self.params.theta * a))
        term2 = -A * (-1/self.params.delta_2) * (
            safe_exp(-self.params.delta_2 * b) - safe_exp(-self.params.delta_2 * a))
        term3 = B * safe_exp((self.params.delta_1 + self.params.theta) * self.T2) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term4 = -B * (1/self.params.delta_1) * (
            safe_exp(self.params.delta_1 * b) - safe_exp(self.params.delta_1 * a))
        return term1 + term2 + term3 + term4
    def dI_s3(self,t):
        return self.params.V(t)-self.params.theta*self.params.I_s2(t)
    def I_s3(self, t,):
        return (((self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)) *
               (safe_exp(-(self.params.delta_2 - self.params.theta) * self.T2 - self.params.theta * t) - safe_exp(-self.params.delta_2 * t)) +
               (self.params.D0 / (self.params.delta_1 + self.params.theta)) * 
               (safe_exp((self.params.delta_1 + self.params.theta) * self.T2 - self.params.theta * t) - 
                                        safe_exp((self.params.delta_1 + self.params.theta) * self.T3 - self.params.delta_1 * t)))
    def integral_I_s3(self, a, b):
        A = (self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)
        B = self.params.D0 / (self.params.delta_1 + self.params.theta)
        term1 = A * safe_exp(-(self.params.delta_2 - self.params.theta) * self.T2) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term2 = -A * (-1/self.params.delta_2) * (
            safe_exp(-self.params.delta_2 * b) - safe_exp(-self.params.delta_2 * a))
        term3 = B * safe_exp((self.params.delta_1 + self.params.theta) * self.T2) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term4 = -B * safe_exp((self.params.delta_1 + self.params.theta) * self.T3 - self.params.delta_1 * a) * (
            1/self.params.delta_1) * (safe_exp(self.params.delta_1 * b) - safe_exp(self.params.delta_1 * a))
        return term1 + term2 + term3 + term4
    def dI_s4(self,t):
        return self.params.eeta* self.params.R0*np.exp(-self.params.delta3*t) - self.params.theta*self.params.I_s4(t)
    def I_s4(self, t):
        return (((self.params.eeta * self.params.R0) /(self.params.delta_3 - self.params.theta)) * 
               (safe_exp(-(self.params.delta_3 - self.params.theta) * self.T4 - self.params.theta * t) - safe_exp(-self.params.delta_3 * t)) + 
               ((self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)) * 
               (safe_exp(-(self.params.delta_2 - self.params.theta) * self.T2 - self.params.theta * t) - 
                        safe_exp(-(self.params.delta_2 - self.params.theta) * self.T4 - self.params.theta * t)) + 
               (self.params.D0 / (self.params.delta_1 + self.params.theta)) * 
               (safe_exp((self.params.delta_1 + self.params.theta) * self.T2 - self.params.theta * t) - 
                        safe_exp((self.params.delta_1 + self.params.theta) * self.T3 - self.params.delta_1 * t)))
    def integral_I_s4(self, a, b):
        A = (self.params.eeta * self.params.R0) / (self.params.delta_3 - self.params.theta)
        B = (self.params.K * self.params.P0) / (self.params.delta_2 - self.params.theta)
        C = self.params.D0 / (self.params.delta_1 + self.params.theta)
        term1 = A * safe_exp(-(self.params.delta_3 - self.params.theta) * self.T4) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term2 = -A * (-1/self.params.delta_3) * (
            safe_exp(-self.params.delta_3 * b) - safe_exp(-self.params.delta_3 * a))
        term3 = B * safe_exp(-(self.params.delta_2 - self.params.theta) * self.T2) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term4 = -B * safe_exp(-(self.params.delta_2 - self.params.theta) * self.T4) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term5 = C * safe_exp((self.params.delta_1 + self.params.theta) * self.T2) * (
            -1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term6 = -C * safe_exp((self.params.delta_1 + self.params.theta) * self.T3 - self.params.delta_1 * a) * (
            1/self.params.delta_1) * (safe_exp(self.params.delta_1 * b) - safe_exp(self.params.delta_1 * a))
        return term1 + term2 + term3 + term4 + term5 + term6
    def dI_s5(self,t):
        return -self.params.D0*np.exp(self.params.delta1*t)-self.params.theta*self.params.I_s5(t)
    def I_s5(self, t):
        return ((self.params.D0 / (self.params.delta_1 + self.params.theta)) * (safe_exp((self.params.delta_1 + self.params.theta) * self.params.T - self.params.theta * t) 
                                                            - safe_exp(self.params.delta_1 * t)))
    def integral_I_s5(self, a, b):
        term1 = (self.params.D0 / (self.params.delta_1 + self.params.theta)) * safe_exp(
            (self.params.delta_1 + self.params.theta) * self.params.T) * (-1/self.params.theta) * (safe_exp(-self.params.theta * b) - safe_exp(-self.params.theta * a))
        term2 = -(self.params.D0 / (self.params.delta_1 * (self.params.delta_1 + self.params.theta))) * (
            safe_exp(self.params.delta_1 * b) - safe_exp(self.params.delta_1 * a))
        return term1 + term2
    def total_serviceable(self):
        return (self.integral_I_s1(self.T1,self.T2)+self.integral_I_s2(self.T2,self.T3)+self.integral_I_s3(self.T3,self.T4)+
                self.integral_I_s4(self.T4,self.T5)+self.integral_I_s5(self.T5,self.params.T))
    
class WasteInventory:
    def __init__(self, T1, T4, params: ProblemParameters, ct: CommonTime, bt: BackorderTime, inventory = "backorder"):
        self.params = params
        self.T1 = T1
        self.T4 = T4
        self.T2=bt.solve_T2(T1)
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T5=ct.solve_T5(T1,T4)
        self.w   =  params.w_z * (1 - params.w_f + params.w_f)
    def I_w0(self, t):
        return ((self.w * self.params.P0) / self.params.delta_2) * (np.exp(-self.params.delta_2 * self.T1) - np.exp(-self.params.delta_2 * t)) + \
               np.sum(self.params.w_ri) * (t - self.T1)
    def integral_I_w0(self, a, b):
        A = (self.w * self.params.P0) / self.params.delta_2
        term1 = A * safe_exp(-self.params.delta_2 * self.T1) * (b - a)
        term2 = (A / self.params.delta_2) * (safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b))
        term3 = (np.sum(self.params.w_ri) / 2) * (b**2 - a**2)
        term4 = -np.sum(self.params.w_ri) * self.T1 * (b - a)
        return term1 + term2 + term3 + term4
    def I_w1(self, t):
        return ((self.w * self.params.P0) / self.params.delta_2) * (np.exp(-self.params.delta_2 * self.T1) - np.exp(-self.params.delta_2 * t)) + \
               np.sum(self.params.w_ri) * (self.T3 - self.T1)
    def integral_I_w1(self, a, b):
        A = (self.w * self.params.P0) / self.params.delta_2
        term1 = A * safe_exp(-self.params.delta_2 * self.T1) * (b - a)
        term2 = (A / self.params.delta_2) * (safe_exp(-self.params.delta_2 * a) - safe_exp(-self.params.delta_2 * b))
        term3 = np.sum(self.params.w_ri) * (self.T3 - self.T1) * (b - a)
        return term1 + term2 + term3
    def I_w2(self, t):
        return (((1 - self.params.eeta) * self.params.R0) / (self.params.delta_3)) * \
               (np.exp(-self.params.delta_3 * self.T4) - np.exp(-self.params.delta_3 * t)) + \
               ((self.w * self.params.P0) / self.params.delta_2) * \
               (np.exp(-self.params.delta_2 * self.T1) - np.exp(-self.params.delta_2 * self.T4)) + \
               np.sum(self.params.w_ri) * (self.T3 - self.T1)
    def integral_I_w2(self, a, b):
        A = ((1 - self.params.eeta) * self.params.R0) / (self.params.delta_3)
        B = (self.w * self.params.P0) / self.params.delta_2
        term1 = A * safe_exp(-self.params.delta_3 * self.T4) * (b - a)
        term2 = (A / self.params.delta_3) * (safe_exp(-self.params.delta_3 * a) - safe_exp(-self.params.delta_3 * b))
        term3 = B * (safe_exp(-self.params.delta_2 * self.T1) - safe_exp(-self.params.delta_2 * self.T4)) * (b - a)
        term4 = np.sum(self.params.w_ri) * (self.T3 - self.T1) * (b - a)
        return term1 + term2 + term3 + term4
    def total_waste(self):
        waste_w0 = self.integral_I_w0(self.T1, self.T3)
        waste_w1 = self.integral_I_w1(self.T3, self.T4)
        waste_w2 = self.integral_I_w2(self.T4, self.T5)
        total_waste = waste_w0 + waste_w1 + waste_w2
        return total_waste