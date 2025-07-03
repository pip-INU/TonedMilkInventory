from utilities import *
from data import ProblemParameters as pp
from inventorymodelbackorder import ServiceableInventory
import inventorytime
def time(T1,T4, inventory):
    ct = inventorytime.CommonTime()
    bt = inventorytime.BackorderTime()
    T2 = bt.solve_T2(T1)
    T3 = ct.solve_T3(inventory,T1,T4)
    T5 = ct.solve_T5(T1,T4)
    return T2,T3,T5
def timeconstraints(T1, T4):
        T2,T3,T5 = time(T1,T4, inventory="backorder")
    
        if 0<T1<T2<T3<T4<T5<24:
            return True
        else :
            return False
def quantityconstraint(T1,T4, x):
        s = ServiceableInventory(T1,T4)
        if pp.minQ <= s.total_serviceable()/x <= pp.maxQ:
            return True
        else:
            return False 
