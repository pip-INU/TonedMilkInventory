def constraints(self, T1, T4):
        T2 = self.solve_T2(T1)
        T3 = self.solve_T3(T1, T4)
        T5 = self.solve_T5(T1, T4)
    
        if 0<T1<T2<T3<T4<T5<self.params['T']:
            return True
        else :
            return False