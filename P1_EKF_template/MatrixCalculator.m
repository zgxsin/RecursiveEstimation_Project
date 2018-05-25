close all 
clear all
clc

syms q_1 q_2 q_3 q_4 q_5 q_6 a_1 a_2 dragCoefficient rudderCoefficient
syms p_11 p_12 p_13 p_14 p_15 p_16 ...
     p_21 p_22 p_23 p_24 p_25 p_26 ...
     p_31 p_32 p_33 p_34 p_35 p_36 ...
     p_41 p_42 p_43 p_44 p_45 p_46 ...
     p_51 p_52 p_53 p_54 p_55 p_56 ... 
     p_61 p_62 p_63 p_64 p_65 p_66
 syms v_d v_r v_b

a33 = -2*cos(q_5)*dragCoefficient*(1+0)*q_3;
a34 = -2*cos(q_5)*dragCoefficient*(1+0)*q_4;
a35 = -sin(q_5)*(tanh(a_1)-dragCoefficient*(q_3^2+q_4^2)*(1+0));
a43 = -2*sin(q_5)*dragCoefficient*(1+0)*q_3;
a44 = -2*sin(q_5)*dragCoefficient*(1+0)*q_4;
a45 = cos(q_5)*(tanh(a_1)-dragCoefficient*(q_3^2+q_4^2)*(1+0));

A = [0,0,1,0,0,0;
     0,0,0,1,0,0;
     0,0,a33,a34,a35,0;
     0,0,a43,a44,a45,0;
     0,0,0,0,0,0;
     0,0,0,0,0,0];

P =  [ p_11 p_12 p_13 p_14 p_15 p_16 ;
     p_21 p_22 p_23 p_24 p_25 p_26 ;
     p_31 p_32 p_33 p_34 p_35 p_36 ;
     p_41 p_42 p_43 p_44 p_45 p_46 ;
     p_51 p_52 p_53 p_54 p_55 p_56 ; 
     p_61 p_62 p_63 p_64 p_65 p_66 ];
 
Q = [v_d 0 0; 
     0 v_r 0;
     0 0 v_b];

l31 = cos(q_5)*(-dragCoefficient*(q_3^2+q_4^2));
l41 = sin(q_5)*(-dragCoefficient*(q_3^2+q_4^2));
l52 = rudderCoefficient*a_2;
 
L = [0,0,0;
     0,0,0;
     l31,0,0;
     l41,0,0;
     0,l52,0;
     0,0,1];
 
 A*P+P*A+L*Q*transpose(L)
  