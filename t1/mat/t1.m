close all
clear all
format long

% Circuit T1 - Theoretical analysis

Z = 0.000000000000
O = 1.000000000000
R1 = 1.04001336091e+03
R2 = 2.04372276851e+03
R3 = 3.11359737601e+03
R4 = 4.17085404861e+03
R5 = 3.02859283303e+03
R6 = 2.07054576700e+03
R7 = 1.01835949725e+03
Va = 5.20102702949 
Id = 1.00460501759e-03
Kb = 7.19043597753e-03
Kc = 8.06397385506e+03


% Mesh analysis

A = [R4+R3+R1,-R3,-R4;-Kb*R3,Kb*R3-O,Z;R4,Z,-R6-R7+Kc-R4]
B = [-Va; Z; Z]

C = A\B

IA = C(1,1)
IB = C(2,1)
IC = C(3,1)

V1 = R1*IA
V2 = V1+R2*IB
V3 = -Va
V4 = V1-R3*(-IA+IB)
V5 = V4-R5*(IB-Id)
V6 = V3-R6*IC
V8 = V6-R7*IC

Vn = [R1*IA, V1+R2*IB, -Va, V1-R3*(-IA+IB), V4-R5*(IB-Id), V3-R6*IC, V6-R7*IC]
NodeN = [1:6,8]
mT = table( NodeN, Vn)

% Node analysis

D = [(-O/R2)-(O/R3)-(O/R1),O/R2,Z,O/R3,Z,Z,Z;Kb+(O/R2),-O/R2,Z,-Kb,Z,Z,Z;Z,Z,O,Z,Z,Z,Z;O/R3,Z,O/R4,(-O/R4)+(-O/R3)-(O/R5),O/R5,O/R7,-O/R7;Kb,Z,Z,(-O/R5)-(Kb),O/R5,Z,Z;Z,Z,O/R6,Z,Z,(-O/R6)-(O/R7),O/R7;Z,Z,Kc/R6,-O,Z,-Kc/R6,O]
E = [Z;Z;-Va;Id;Id;Z;Z]

F = D\E
