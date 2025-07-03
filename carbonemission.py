import numpy as np
from data import ProblemParameters
from utilities import *
from inventorytime import CommonTime
import wastegenerated
import numpy as np
class CalculateCarbonEmission:
    def __init__(self, T1, T4, inventory, params: ProblemParameters, ct:CommonTime, wg: wastegenerated):
        self.params = params 
        self.T1=T1
        self.T3 = ct.solve_T3(T1,T4, inventory)
        self.T4=T4
        self.T5=ct.solve_T5(T1,T4)
        self.total_waste= wg.UsefulWaste
        self.unuse_waste_generate=  wg.UselessWaste()
        self.carbon_cap = self.params.emission_cap
    def processemission(self):
        return(self.params.e_p * (self.params.P0 * safe_exp(-self.params.delta_2 * (self.T4 - self.T1)) + 
                                self.params.e_rp * self.params.R0 * safe_exp(-self.params.delta_3 * (self.T5 - self.T4))))
    def coldstorageemission(self):
        return self.params.e_c * ((self.T3-self.T1) + (self.params.T-self.T1) + (self.T4-self.T1) + (self.params.T-self.T1))
    def wasteprocessemission(self):
        return self.params.e_w * self.total_waste 
    def dumpwasteemission(self):
        return self.params.e_dump * max(0, self.unuse_waste_generate)
    def calculate_total_emission(self):
        return (self.processemission() + self.coldstorageemission() + 
                    self.wasteprocessemission() + self.dumpwasteemission())
    def emissiondifference(self):
        return self.carbon_cap - self.calculate_total_emission()
    def carbon_credits(self):
        return self.emissiondifference() * 1 # 1 ton carbon = 1 credit
    