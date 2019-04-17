function main_JPDAF

clc;clear all;close all;

n=20;                                                                             %采样次数
T=1;                                                                               %T为采样间隔
MC_number=1;                                                             %Monte Carlo仿真次数
c=2;                                                                                %目标个数
target_position=[1500 300 500 400; 
                             500 400 1500 300];                              %目标的起始位置和速度(m,m/s)                   
JPDAF(target_position,n,T,MC_number,c);                                            
%调用子函数PDAF 