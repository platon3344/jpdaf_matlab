function main_JPDAF

clc;clear all;close all;

n=20;                                                                             %��������
T=1;                                                                               %TΪ�������
MC_number=1;                                                             %Monte Carlo�������
c=2;                                                                                %Ŀ�����
target_position=[1500 300 500 400; 
                             500 400 1500 300];                              %Ŀ�����ʼλ�ú��ٶ�(m,m/s)                   
JPDAF(target_position,n,T,MC_number,c);                                            
%�����Ӻ���PDAF 