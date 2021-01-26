close all;
clear all;
clc;

n = Nondominated();
n.T= 1;
n.A = [1,1;1,0];
n.B = [0;1];
n.x0 = [1;0];

par = n.paretoT();
n.plotPareto(par)

WF = n.searchBack();

[t,x] = n.timeDepODE(WF(1), WF(2));
figure();
plot(t,x)

