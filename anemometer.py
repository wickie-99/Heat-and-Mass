import sympy as sp
import numpy as np

class fluid():
    def __init__(self, Ts):
        self.Ts = Ts
        self.Tf = 1 
        self.density = 1
        self.viscosity = 1
        self.conductivity = 1
        self.Pr = 1
        self.Pr_s = 1

    def update_properties(self):
        self.Tf = (self.Ts+300)/2
        self.density = 1.1614 + (self.Tf-300)/(350-300)*(0.995-1.1614)
        self.viscosity = (184.6 + (self.Tf-300)/(350-300)*(208.2-184.6))*10**-7
        self.conductivity = (26.3 + (self.Tf-300)/(350-300)*(30.0-26.3))*10**-3
        self.Pr = (0.707 + (self.Tf-300)/(350-300)*(0.70-0.707))
        self.Pr_s = 0.707
        
class Zukauskas():
    def __init__(self, Re):
        self.C = 1 
        self.m = 1
        self.n = 1
        self.Re = Re

    def update_Zukauskas(self):
        if (self.Re<=40):       
            self.C = 0.75 
            self.m = 0.4
            self.n = 0.33
        elif (self.Re<=1000):       
            self.C = 0.51 
            self.m = 0.5
            self.n = 0.33
        elif (self.Re<=200000):       
            self.C = 0.26
            self.m = 0.6
            self.n = 0.33
        else:       
            self.C = 0.076
            self.m = 0.7
            self.n = 0.33

I, k,C,m,n,Pr,Prs,A,s,l,Ts, Too = sp.symbols("I, k,C,m,n,Pr,Pr_s,A,s,l,T_s, T_oo")
Re, rho, d, mu, v = sp.symbols("Re, rho, d, mu, v")
dI, dv = sp.symbols("dI, dv")
Q1, Q2 = sp.symbols("Q_1, Q_2")

eq_1 = sp.Eq(I,2*v/m*dI/dv)
eq_2 = sp.Eq(Re,rho*v*d/mu)
eq_3 = sp.Eq(I**2, (k/d)*C*(Re**m)*(Pr**n)*((Pr/Prs)**(1/4))*(((sp.pi*d*l)*(sp.pi*d**2/4))/(s*l))*(Ts-Too))

values_1 = {v:50,I:0.2,dI:10*10**-6,dv:0.01}       
m_min = sp.solve(eq_1.subs(values_1), m)
print(m_min)
Re_min = 40


min_diameters = []

for Temp in range(312,327):
    air = fluid(Temp)
    air.update_properties()
    values_2 = {Re:40, rho:air.density,v:50, mu:air.viscosity}
    d_min = sp.solve(eq_2.subs(values_2),d)
    min_diameters.append((Temp, d_min[0]*10**6))
    
print(min_diameters)
        
Q1 = I**2*s*l/(np.pi*d**4/4)
Q2 = (k/d)*C*(Re**m)*(Pr**n)*((Pr/Prs)**(1/4))*(np.pi*d*l)*(Ts-Too)

curernt = 0.2
velocity = 50

values = {I:curernt, Too:300,s:10.6*10**-8,l:0.005, v:velocity}

diameter = 13
Temperature = 313
Error = 1*10**10
           
temp_accuracy = 100
dia_accuracy = 100

for temp in np.linspace(313, 327, temp_accuracy):
    air = fluid(temp)
    for dia in np.linspace(13,50, dia_accuracy):
        print(temp, dia, Error)
        Re_num = air.density*velocity*dia/air.viscosity
        zuka = Zukauskas(Re_num)
        Error_temp = abs(Q1.subs(values).subs({d:dia}) - (Q2.subs(values).subs({k:air.conductivity,rho:air.density,mu:air.viscosity,Pr:air.Pr, Prs:air.Pr_s , C:zuka.C,m:zuka.m,n:zuka.n, Ts:air.Ts, Re:Re_num})))
        if Error_temp<Error:
            Error = Error_temp
            diameter = dia
            Temperature = temp
            print(temp)

print(diameter, Temperature, Error)
