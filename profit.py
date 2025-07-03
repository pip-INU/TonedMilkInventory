from totalcost import CalculateBackorderTotalCost
from revenue import CalculateRevenue
class CalculateProfit:
    def __init__(self, T1, T4, x, CP:CalculateBackorderTotalCost, SP:CalculateRevenue):
        self.T1 = T1
        self.T4 = T4
        self.x = x
        self.CP = CP
        self.SP = SP
    def calculate_profit(self):
        return self.SP.total_revenue() - self.CP.TotalCost()