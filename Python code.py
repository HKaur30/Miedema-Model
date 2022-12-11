from sympy import *
import numpy as np
import sympy as sym
sym.init_printing()
import re
import sympy
def file_read(fname, EL_1, EL_2):
E1data_array = []
E2data_array = []
with open(fname) as f:
flag=0
index=0
#Content_list is the list that contains the read lines.
for line in f:
index=index+1
element_1 = re.match(EL_1,line)
if element_1:
flag=1
break
if flag == 0:
print ("Element doesnot exist")
else:
E1data_array=line.split()
for i in range(1, len(E1data_array)):
E1data_array[i-1] = float(E1data_array[i])
with open(fname) as f:
flag=0
index=0
#Content_list is the list that contains the read lines.
for Line in f:
index=index+1
element_2= re.match(EL_2,Line)
if element_2:
flag=1
break
if flag == 0:
print ("Element doesnot exist")
else:
E2data_array=Line.split()
for i in range(1, len(E2data_array)):
E2data_array[i-1] = float(E2data_array[i])
return E1data_array,E2data_array
def P(E1data,E2data):
P_Value = [14.2, 12.35, 10.7]
if (E1data[5] + E2data[5]) == 2:
# Both elementA and elementB are Transition Metals.
return P_Value[0]
elif (E1data[5] + E2data[5]) == 1:
# Only one of elementA and elementB are Transition Metals.
return P_Value[1]
else:
# Neither elementA nor elementB are Transition Metals.
return P_Value[2]
print(P_Value )
def RtoP(E1,E2):
tmrange = []
tmrange.extend(list(range(20, 30)))
tmrange.extend(list(range(38, 48)))
tmrange.extend(list(range(56, 58)))
tmrange.extend(list(range(72, 80)))
tmrange.extend([90, 92, 94])
nontmrange = []
nontmrange.extend(list(range(3, 8)))
nontmrange.extend(list(range(11, 16)))
nontmrange.extend([19])
nontmrange.extend(list(range(30, 34)))
nontmrange.extend([37])
nontmrange.extend(list(range(48, 52)))
nontmrange.extend([55])
nontmrange.extend(list(range(80, 84)))
if (E1[3] in tmrange) and (E2[3] in nontmrange):
RtoP = E1[6] * E2[6]
elif (E1[3] in nontmrange) and (E2[3] in tmrange):
RtoP = E1[6] * E2[6]
else:
RtoP = 0.0
return RtoP
def gamma(E1,E2,P_val,RtoP_val):
QtoP = 9.4 # Constant from Miedema's Model.
phi = [E1[0], E2[0]]
rho = [E1[1], E2[1]]
d_phi = phi[0] - phi[1]
d_rho = rho[0] - rho[1]
m_rho = (1/2)*(1/(rho[0] + rho[1]))
gamma = (P_val * (QtoP * (d_rho ** 2 )- (d_phi ** 2) - RtoP_val) /
m_rho)
return int(round(gamma))
def a_A(E1,E1d,E2d):
return pick_a(E1,E1d,E2d)
def a_B(E2,E1d,E2d):
return pick_a(E2,E1d,E2d)
def pick_a(elt,E1d,E2d):
"""Choose a value of a based on the valence of element A."""
possible_a = [0.14, 0.1, 0.07, 0.04]
if elt == E1:
params = E1d
else:
params = E2d
if params[4] == 1:
return possible_a[0]
elif params[4] == 2:
return possible_a[1]
elif params[4] == 3:
return possible_a[2]
# elif elementA in ["Ag","Au","Ir","Os","Pd","Pt","Rh","Ru"]:
elif elt in ["Ag", "Au", "Cu"]:
return possible_a[2]
else:
return possible_a[3]
def H_form_ord(E1,E2,Cs,gamma_val,E1_A,E1_B):
"""Calculate the enthalpy of formation for an ordered compound of
elements A
and B with a composition xB of element B."""
vol0_A = E1[2]
vol0_B = E2[2]
phi = [E1[0], E2[0]]
#htrans = [E1[7], E2[7]]
# Determine volume scale parameter a.
# Calculate surface concentrations using original volumes.
c_S_A = Cs * vol0_A / (Cs * vol0_A + (1-Cs) * vol0_B)
c_S_B =(1-Cs) * vol0_B / (Cs * vol0_A + (1-Cs) * vol0_B)
# Calculate surface fractions for ordered compounds using original
volumes.
f_BA = c_S_B * (1 + 8 * ((c_S_A * c_S_B) ** 2))
f_AB = c_S_A * (1 + 8 *( (c_S_A * c_S_B) ** 2))
# Calculate new volumes using surface fractions (which use
original
# volumes).
vol_A = vol0_A * (1 + E1_A * f_BA * (phi[0] - phi[1]))
vol_B = vol0_B * (1 + E1_B * f_AB * (phi[1] - phi[0]))
#Recalculate surface concentrations using new volumes.
c_S_A = (1 - Cs) * vol_A / ((1 - Cs) * vol_A + Cs * vol_B)
c_S_B = Cs * vol_B / ((1 - Cs) * vol_A +Cs * vol_B)
# Recalculate surface fractions for ordered compounds using new
volumes.
f_BA = c_S_B * (1 + 8 * ((c_S_A * c_S_B) ** 2))
f_AB = c_S_A * (1 + 8 * ((c_S_A * c_S_B) ** 2))
#D_htrans = Cs * htrans[1] + (1 - Cs) * htrans[0]
Sc = 1 – ((c_S_B (vol_A - vol_B ))/ (((c_S_A * vol_A)+ (c_S_B *
vol_B))^2)
H_ord = (
gamma_val
*vol_A
*Cs
*f_BA
)*Sc
return round(H_ord ,2)
if __name__=="__main__":
El_1=input('Enter Solute Element ')
El_2=input('Enter SolventElement ')
CS=float(input('Enter concentration of solute '))
E1Data, E2Data = file_read("C:\\Users\\User\\Desktop\\DATASET.txt", El_1,
El_2)
P_Val =P(E1Data,E2Data)
RtoP_Val=RtoP(E1Data,E2Data)
gamma_Val= gamma(E1Data,E2Data,P_Val,RtoP_Val)
E1_a = a_A(El_1,E1Data,E2Data)
E1_b = a_B(El_1,E1Data,E2Data)
H_Val=H_form_ord(E1Data,E2Data,CS,gamma_Val,E1_a,E1_b)
print (H_Val)