from utilities import *
from data import ProblemParameters 
from inventorytime import CommonTime
class ServiceableInventory:
    def __init__(self, T1, T4, params: ProblemParameters, ct: CommonTime, inventory = "buffer"):
        self.params = params
        self.T1 = T1
        self.T4 = T4
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T5=ct.solve_T5(T1,T4)
        
    def I_s0(self, t):
        return self.S0 * safe_exp(-self.theta * t) + \
                (self.D0 / (self.theta + self.delta1)) * (safe_exp(-self.theta * t) - safe_exp(self.delta1 * t))

    def integral_I_s0(self, T1):
        return self.S0 * (1 - safe_exp(-self.theta * T1)) / self.theta + \
                (self.D0 / (self.theta + self.delta1)) * ((1 - safe_exp(-self.theta * T1)) / self.theta - (safe_exp(self.delta1 * T1) - 1) / self.delta1)

    def I_s1(self, t, T1):
        return self.S0 * safe_exp(-self.theta * t) + \
                (self.K * self.P0) / (self.lam * (self.delta2 - self.theta)) * (safe_exp(-(self.delta2 - self.theta) * T1 - self.theta * t) - safe_exp(-self.delta2 * t)) + \
                (self.D0 / (self.theta + self.delta1)) * (safe_exp(-self.theta * t) - safe_exp(self.delta1 * t))

    def I_s2(self, t, T1, T3):
        return self.S0 * safe_exp(-self.theta * t) + \
                (self.K * self.P0) / (self.lam * (self.delta2 - self.theta)) * (safe_exp(-(self.delta2 - self.theta) * T1 - self.theta * t) - safe_exp(-self.delta2 * t)) + \
                (self.D0 / (self.delta1 + self.theta)) * (safe_exp(-self.theta * t) - safe_exp(self.delta1 * T1 + self.theta * (T3 - t)))

    def I_s3(self, t, T1, T3, T4):
        return self.S0 * safe_exp(-self.theta * t) + \
                (self.eta * self.self.params.R0) / (self.lam * (self.delta3 - self.theta)) * (safe_exp(-(self.delta3 - self.theta) * T4 - self.theta * t) - safe_exp(-self.delta3 * t)) + \
                (self.K * self.P0) / (self.lam * (self.delta2 - self.theta)) * (safe_exp(-(self.delta2 - self.theta) * T1 - self.theta * t) - safe_exp(-(self.delta2 - self.theta) * T4 - self.theta * t)) + \
                (self.D0 / (self.delta1 + self.theta)) * (safe_exp(self.theta * (T4 - t)) - safe_exp(self.delta1 * T1 + self.theta * (T4 - T3 - t)))

    def I_s4(self, t, T4, T5, T3, T1):
        A = self.S0 * safe_exp(-self.theta * T5) + \
                (self.eta * self.self.params.R0) / (self.lam * (self.delta3 - self.theta)) * (safe_exp(-(self.delta3 - self.theta) * T4 - self.theta * T5) - safe_exp(-self.delta3 * T5)) + \
                (self.K * self.P0) / (self.lam * (self.delta2 - self.theta)) * (safe_exp(-(self.delta2 - self.theta) * T1 - self.theta * T5) - safe_exp(-(self.delta2 - self.theta) * T4 - self.theta * T5)) + \
                (self.D0 / (self.delta1 + self.theta)) * safe_exp(self.theta * (T4 - T5))
        B = self.S0 * safe_exp(self.theta * (self.T - T5)) + \
                (self.D0 / (self.delta1 + self.theta)) * (safe_exp((self.delta1 + self.theta) * self.T - self.theta * T5) - safe_exp(self.delta1 * T5))
        condition = self.theta * A > self.xi
        if condition:
            return self.S0 * safe_exp(self.theta * (self.T - t)) + \
                    (self.D0 / (self.delta1 + self.theta)) * (safe_exp((self.delta1 + self.theta) * self.T - self.theta * t) - safe_exp(self.delta1 * t))
        else:
            return self.S0 * safe_exp(2 * self.theta * (self.T - t)) + \
                    (self.D0 / (self.delta1 + 2 * self.theta)) * (safe_exp((self.delta1 + 2 * self.theta) * self.T - 2 * self.theta * t) - safe_exp(self.delta1 * t))
        
class WasteInventory:
    def __init__(self, T1, T4, params: ProblemParameters, ct:CommonTime, inventory = "buffer"):
        self.params = params
        self.w   =  params.w_z * (1 - params.w_f + params.w_f)
        self.T1=T1
        self.T3=ct.solve_T3(inventory,T1,T4)
        self.T4=T4
    def dI_w0(self):
        return np.sum(self.params.w_ri)
    def I_w0(self,t):
        return np.sum(self.params.w_ri)*t
    def dI_w1(self,t):
        return self.params.w * self.params.P0*np.exp(-self.params.delta_2*t)+np.sum(self.params.w_ri)
    def I_w1(self,t):
        return ((self.w*self.params.P0)/self.params.delta_2)*(np.exp(-self.params.delta_2*self.T1)-np.exp(-self.params.delta_2*t)) +np.sum(self.params.w_ri)*(t-self.T1)
    def dI_w2(self,t):
        return self.params.w * self.params.P0*np.exp(-self.params.delta_2*t)
    def I_w2(self,t):
        return  ((self.w*self.params.P)/self.params.delta_2)*(np.exp(-self.params.delta_2*self.T1)-np.exp(-self.params.delta_2*t)) +np.sum(self.params.w_ri)*(self.T3-self.T1)
    def dI_w3(self,t):
        return ((1-self.params.eeta)*self.params.R0*np.exp(-self.params.delta_3*t))
    def I_w3(self,t):
        return ((((1-self.params.eeta)*self.params.R0)/self.params.delta_3)*(np.exp(-self.params.delta_3*self.T4)-np.exp(-self.params.delta_3*t)) 
                    + ((self.w*self.params.P)/self.params.delta_2)*(np.exp(-self.params.delta_2*self.T1)-np.exp(-self.params.delta_2*t)) +np.sum(self.params.w_ri)*(self.T3-self.T1))