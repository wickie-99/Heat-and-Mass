import sympy as sp
import numpy as np

def property(T):
        Tf = (T+300)/2
        density = 1.1614 + (Tf-300)/(350-300)*(0.995-1.1614)
        viscosity = (184.6 + (Tf-300)/(350-300)*(208.2-184.6))*10**-7
        conductivity = (26.3 + (Tf-300)/(350-300)*(30.0-26.3))*10**-3
        Pr = (0.707 + (Tf-300)/(350-300)*(0.70-0.707))
        Pr_s = 0.707
        return(density,viscosity,conductivity,Pr, Pr_s)

def zuka(Re_num):
    Re = Re_num
    if (Re<=40):       
        C = 0.75 
        m = 0.4
        n = 0.33
    elif (Re<=1000):       
        C = 0.51 
        m = 0.5
        n = 0.33
    elif (Re<=200000):       
        C = 0.26
        m = 0.6
        n = 0.33
    else:       
        C = 0.076
        m = 0.7
        n = 0.33
    return (C,m,n)

I, k,C,m,n,Pr,Prs,A,s,l,Ts, Too = sp.symbols("I, k,C,m,n,Pr,Pr_s,A,s,l,T_s, T_oo")
Re, rho, d, mu, v = sp.symbols("Re, rho, d, mu, v")
dI, dv = sp.symbols("dI, dv")
Q1, Q2 = sp.symbols("Q_1, Q_2")
Err = sp.Symbol("Err")

eq_sensitivity = sp.Eq(I,2*v/m*dI/dv)
eq_Re = sp.Eq(Re,rho*v*d/mu)
eq_Q1 = sp.Eq(Q1, I**2*s*l/(sp.pi/4*d**2))
eq_Q2 = sp.Eq(I**2, (k/d)*C*(Re**m)*(Pr**n)*((Pr/Prs)**(1/4))*(sp.pi*d*l)*(Ts-Too))

l_num = 0.01
s_num = 10.6 * 10**-8
Too_num = 300


v_max = 50
I_max = 0.2

values_1 = {v:v_max,I:I_max,dI:10*10**-6,dv:0.01}       
m_min = round(sp.solve(eq_sensitivity.subs(values_1), m)[0],1)

m_Re = 0
Re_min = 0

while m_Re < m_min:
    z = zuka(Re_min)[1]
    m_Re = z
    Re_min += 1
Re_min -= 1

print(Re_min)
min_diameters = []

for Temp in range(312,327):
    prop = property(Temp)
    values = {Re:Re_min, rho:prop[0],v:v_max, mu:prop[1]}
    d_min = sp.solve(eq_Re.subs(values),d)
    min_diameters.append(d_min[0]*10**6)
    
d_min = int(round(max(min_diameters)))
print(d_min)
temp_acc = 100
dia_acc = 100

values_1 = {v:v_max, I:I_max, l:l_num,s:s_num, Too:Too_num, sp.pi:np.pi}

Error = 10**10
diameter = 0
temperature = 0

for temp in np.linspace(313,327,temp_acc):
    prop = property(temp)
    values_2 = {k:prop[2], mu:prop[1], rho: prop[0], Pr:prop[3], Prs:prop[4]}
    for dia in np.linspace(d_min,50,dia_acc):
        dia = dia * 10**-6
        Re_num = prop[0]*v_max*dia/prop[1]
        z = zuka(Re_num)

        values_3 = {C:z[0], m:z[1], n:z[2], Re:Re_num, Ts:temp, d:dia}
        Error_temp = abs((eq_Q1.rhs - eq_Q2.rhs).subs(values_1).subs(values_2).subs(values_3))
        if Error_temp< Error:
            diameter = dia
            temperature = temp
            Error = Error_temp

print("diameter = ", end="")
print(diameter*10**6) 
print("Temperature = ", end="")
print(temperature)
print("Error = ", end="")
print(Error)