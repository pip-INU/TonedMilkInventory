from data import ProblemParameters
from utilities import *
import inventorymodelbackorder
import inventorymodelbuffer 
import inventorytime
import commoninventory
import delay
import carbonemission
def inventoryobject(T1,T4):
    RI1 = commoninventory.RawInventory(T1,T4, inventory="backorder")
    SI1 = inventorymodelbackorder.ServiceableInventory(T1,T4)
    RWI1 = commoninventory.ReprocessableInventory(T1,T4, inventory="backorder")
    WI1 = inventorymodelbackorder.WasteInventory(T1,T4)
    RI2 = commoninventory.RawInventory(T1,T4, inventory="buffer")
    SI2 = inventorymodelbuffer.ServiceableInventory(T1,T4)
    RWI2 = commoninventory.ReprocessableInventory(T1,T4, inventory="buffer")
    WI2 = inventorymodelbuffer.WasteInventory(T1,T4)
    return RI1,SI1, RWI1, WI1, RI2, SI2, RWI2, WI2
def time(T1,T4, inventory):
    ct = inventorytime.CommonTime()
    T3 = ct.solve_T3(inventory,T1,T4)
    T5 = ct.solve_T5(T1,T4)
    return T3,T5
class CalculateBackorderTotalCost:
    def __init__(self, T1,T4, x, params: ProblemParameters):
        self.params = params
        self.T1=T1
        self.T4=T4
        bt = inventorytime.BackorderTime()
        self.T2 = bt.solve_T2(T1)
        self.T3, self.T5 = time(T1,T4, inventory="backorder")
        self.totalprocured = self.params.mu0 * (self.T3-T1)
        self.totaltoreprocess = inventoryobject(T1,T4)[2].integral_I_rw0(self.T1,self.T4)
        self.totalraw = inventoryobject(T1,T4)[0].total_raw()
        self.totalservice = inventoryobject(T1,T4)[1].total_serviceable()
        self.totalreprocess = inventoryobject(T1,T4)[2].total_reprocess()
        self.totalwaste = inventoryobject(T1,T4)[3].total_waste()
        self.totaldelay = np.sum(delay.generate_random_delays(self.params.n))
        self.totalreturn =   np.sum(self.params.w_ri) 
        self.backorder = inventoryobject(T1,T4)[1].integral_I_s0(0,T1)
        self.x = x
        self.carboncredit = carbonemission.CalculateCarbonEmission(T1,T4, inventory= "backorder").carbon_credits()
    def procurementCost(self):
        return self.params.CP * self.totalprocured
    def processingCost(self):
        return self.params.c_p * self.totalprocured
    def reprocessingCost(self):
        return self.params.c_r * self.totaltoreprocess
    def rawHoldingCost(self):
        return self.params.c_raw * self.totalraw
    def serviceableHoldingCost(self):
        return self.params.c_s * self.totalservice
    def reprocessHoldingCost(self):
        return self.params.c_rw *self.totalreprocess
    def wasteHoldingCost(self):
        return self.params.c_w * self.totalwaste
    def fixedCost(self):
        processing_setup_cost = self.params.setup_p
        setup_to_value_add_waste = self.params.setup_w
        cold_storage_setup = self.params.c_s * 4 # 4 inventories are end to end cold storage
        return processing_setup_cost + setup_to_value_add_waste + cold_storage_setup
    def packagingCost(self):
        total_packets = self.params.P0*safe_exp(-self.delta_2*(self.T4-self.T1))//self.params.lmbd
        return self.params.c_pack * total_packets
    def transportCost(self):
        return self.params.c_t*self.totalservice*self.x
    def returnpolicycost(self):
        return 0.7*self.params.sp * self.totalreturn #brought back at discounted price
    def penalty(self):
        penalty_backorder = self.params.p_b * self.backorder
        penalty_delay = self.params.p_d * self.totaldelay
        return penalty_backorder + penalty_delay
    def carbonbalance(self):
        return self.carboncredit * self.params.credit_price
    def TotalCost(self):
        return (self.procurementCost() + self.processingCost() + self.reprocessingCost() +
                self.rawHoldingCost() + self.serviceableHoldingCost() + self.reprocessHoldingCost() +
                self.wasteHoldingCost() + self.fixedCost() + self.packagingCost() +
                self.transportCost() + self.returnpolicycost() + self.penalty() + self.carbonbalance())

    
    