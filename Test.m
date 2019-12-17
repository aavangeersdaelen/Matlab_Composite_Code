%Desciption: Test code for comp and laminate objects 
clear;clc;
load('UniCarbonEpoxy.mat');
nasa=Comp('Nasa',20.01e7,1.301e6,1.001e6,.3,.005);
nasa2=Comp('Nasa2',30.01e6,1.301e6,1.001e6,.3,.005);
test=Laminate('Test',[nasa,nasa2,nasa],[0,90,0],0);
printprop(test)
N=[10 30 40]';
M=[30 40 20]';
test=strainNM(test,N,M);

