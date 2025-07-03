from profit import CalculateProfit
from tqdm import tqdm
if __name__ == "__main__":
    results=[]
    for T1 in tqdm(range(0, 25)):
        for T4 in range(0, 25):
            for x in range(4, 30):
                optimizer = CalculateProfit(T1,T4,x)
                if optimizer.constraints(T1,T4) and optimizer.quantityconstraint(T1,T4,x):
                        res = optimizer.calculate_profit()
                        if res:
                            results.append(res)