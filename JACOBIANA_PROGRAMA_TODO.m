clc
close all
clear all

%Defino variables
syms q1 d2 d3 q4 l1 l4
PI=sym('pi'); 

%Matrices a crear  %matrices de denavit 
T0=DenavitH(q1,l1,0,0);
T1=DenavitH(0,d2,0,-PI/2); %Datos de tablita
T2=DenavitH(0,d3,0,0);
T3=DenavitH(q4,l4,0,0);

%HTM
MultDH=simplify(T0*T1*T2*T3) %multiplicar las de denavit 

%Matriz lineal jacobiana de velocidad
X=MultDH(1:3,4); %las de posicion 
Joint=[q1; d2; d3; q4;]%moviimientos de eje z

JV = simplify(jacobian(X,Joint))%derivadas de las de posicion 

%Matriz angular
Y=MultDH(1,1:3); %rotacion 
Y1=MultDH(2,1:3);
Y2=MultDH(3,1:3);
R=[Y;Y1;Y2]; %matriz 3x3 

%Matriz rotacion
r1=R(1,1);
r2=R(1,2);
r3=R(1,3);      %guardar cada posicion en una variable 
r4=R(2,1);
r5=R(2,2);
r6=R(2,3);
r7=R(3,1);
r8=R(3,2);
r9=R(3,3);

Rt= [r1 r4 r7;r2 r5 r8;r3 r6 r9];  %transpuesta

RQ1=diff(R,q1);
RQ2=diff(R,d2);   
RQ3=diff(R,d3);  %diff sirve para derivar 1er argumento q voy a derivar 2do contra q lo voy a derivar
RQ4=diff(R,q4);

Omega1=simplify(RQ1*Rt)
Omega2=simplify(RQ2*Rt)   %J=matriz omega 
Omega3=simplify(RQ3*Rt)   %RQ ya es la R derivada 
Omega4=simplify(RQ4*Rt)


Wz1=Omega1(2,1) %q1
Wz2=Omega2(2,1)
Wz3=Omega3(2,1)
Wz4=Omega4(2,1) 

Wy4=Omega4(1,3)
Wy1=Omega1(1,3)    
Wy2=Omega2(1,3)
Wy3=Omega3(1,3)

Wx4=Omega4(3,2)
Wx1=Omega1(3,2)
Wx2=Omega2(3,2)
Wx3=Omega3(3,2)

Wz= Wz1+Wz2+Wz3+Wz4
Wx= Wx4+Wx3+Wx2+Wx1 
Wy= Wy4+Wy3+Wy2+Wy1

JW=[ Wx1 Wx2 Wx3 Wx4 ;   % si Wx# existe, entonces el valor de (1,#) será Wx
     Wy1 Wy2 Wy3 Wy4 ;  % si Wy# existe, entonces el valor de (2,#) será Wy  
     Wz1 Wz2 Wz3 Wz4]    % si Wz# existe, entonces el valor de (3,#) será Wz
 
 %Jacobiana

JQ=[JV;JW]

%INVERSA
Jtranspuesta= transpose(JQ)
Jpseudo=simplify((((Jtranspuesta)*(JQ))^-1)*(Jtranspuesta),'Steps',2000) %n menor q 6
%si n es = a 6 se invierte JQ y si n es mayor q 6 otra fomrula 