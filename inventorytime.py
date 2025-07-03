import numpy as np
from utilities import *
from data import ProblemParameters
class CommonTime:
    def __init__(self, params: ProblemParameters):
        self.params = params
    def solve_T5(self, T1, T4):
        term1 = (self.params.H * self.params.P0 / self.params.delta_2) * (safe_exp(-self.params.delta_2 * T4) - safe_exp(-self.params.delta_2 * T1))
        term2 = (self.params.delta_3 / self.params.R0) * term1 + safe_exp(-self.params.delta_3 * T4)
        return -safe_log(term2) / self.params.delta_3
    def solve_T3(self, T1, T4, inventory):
        if inventory == "backorder":
            T5 = self.solve_T5(T1, T4)
            bt = BackorderTime()
            T2 = bt.solve_T2(T1)

            part1 = (self.params.eeta * self.params.R0) / (self.params.delta_3 - self.params.theta)
            part2 = self.params.D0 / (self.params.delta_1 + self.params.theta)
            term1 = (part1 / part2) * (safe_exp(-(self.params.delta_3 - self.params.theta) * T4 - self.params.theta * T5) - safe_exp(-self.params.delta_3 * T5))

            part3 = (self.params.K * self.params.P0) /  (self.params.delta_2 - self.params.theta)
            term2 = (part3 / part2) * (safe_exp(-(self.params.delta_2 - self.params.theta) * T2 - self.params.theta * T5) - 
                                        safe_exp(-(self.params.delta_2 - self.params.theta) * T4 - self.params.theta * T5))

            term3 = (safe_exp((self.params.delta_1 + self.params.theta) * T2 - self.params.theta * T5) - 
                        safe_exp((self.params.delta_1 + self.params.theta) * self.params.T - self.params.theta * T5) +
                        safe_exp(self.params.delta_1 * T5))

            log_arg = term1 + term2 + term3
            return (self.params.theta * T5 + safe_log(log_arg)) / (self.params.delta_1 + self.params.theta)
    

class BackorderTime:  
    def __init__(self, params: ProblemParameters):
        self.params = params
    def solve_T2(self, T1):
        numerator   = self.params.K * self.params.P0 * (safe_exp(-self.params.delta_2 * T1) - 1)
        denominator = self.params.delta_2 * (self.params.D0 - self.params.K * self.params.P0)
        if abs(denominator) < 1e-10:
            raise ValueError("Denominator is zero — no feasible T2 solution")
        return numerator / denominator
        
    
    