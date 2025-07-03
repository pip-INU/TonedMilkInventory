from profit import CalculateProfit
from tqdm import tqdm
import pandas as pd
import plots
import constraints
if __name__ == "__main__":
    results=[]
    for T1 in tqdm(range(0, 25)):
        for T4 in range(0, 25):
            for x in range(4, 30):
                optimizer = CalculateProfit(T1,T4,x)
                if constraints.timeconstraints(T1,T4) and constraints.quantityconstraint(T1,T4,x):
                        res = optimizer.calculate_profit()
                        if res:
                            results.append(res)
    
    df = pd.DataFrame(results)
    max_profit_row = df.loc[df['TotalProfit'].idxmax()]

    # Extract corresponding values
    best_x  = max_profit_row['NumberOfDispatches']
    best_T1 = max_profit_row['ProcessingStartTime']
    best_T4 = max_profit_row['ProcessingEndTime']

    print(f"Maximum Profit: Rs.{max_profit_row['TotalProfit']:.2f} at T1={best_T1}, T4={best_T4}, x={best_x}")

    # Now plot surfaces with those values
    plots.plot_T1_T4(x_fixed=best_x)
    plots.plot_T4_x(T1_fixed=best_T1)
    plots.plot_T1_x(T4_fixed=best_T4)
