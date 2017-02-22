clear all
close all
clc

%parte a
A = 0:0.5:10; 
length(A) %longitud del vector 

%parte b
B=3*A;
C=A+B;
F=3+A;
D=A.^2; %elevamos elemento a elemento al cuadrado
E=A*B'; %arreglamos problemas dimensinoales al trasponer B

%parte c
figure; %abre la ventana para cargar el grafico
plot (A,D+F);

%parte d
figure; %abre una nueva ventana para no superponer graficas
plot(A,0.1*sin(sqrt(20/1.5)*A))
xlabel ('tiempo(s)') %nombramos a los ejes
ylabel('posicion(m)')
title('Parte D')

%parte e
figure;
G=0:.1:10; %definimos un nuevo vector con espaciado 0.1
plot(G,0.1*sin(sqrt(20/1.5)*G))%en teoria el espaciado optimo seria un infinitesimal
xlabel ('tiempo(s)')
ylabel('posicion(m)')
title('Parte E')

%parte f
figure;
plotyy(G,0.1*sin(sqrt(20/1.5)*G), G,sqrt(20/1.5)*0.1*cos(sqrt(20/1.5)*G))%graficas de posicion y tiempo superpuestas
pp=plotyy(G,0.1*sin(sqrt(20/1.5)*G), G,sqrt(20/1.5)*0.1*cos(sqrt(20/1.5)*G));
axes(pp(1));
ylabel('posicion(m)');
axes(pp(2));
ylabel('velocidad(m/s)');
xlabel('tiempo(s)');
title('Parte F')

%parte g
figure;
subplot(2,1,1); %seleccionamos la primer subfigura
plot(G,0.1*sin(sqrt(20/1.5)*G)) 
xlabel('tiempo(s)');
ylabel('posicion(m)');
title('Parte G')

subplot(2,1,2); %seleccionamos la segunda subfigura
plot(G,sqrt(20/1.5)*0.1*cos(sqrt(20/1.5)*G)) 
xlabel('tiempo(s)');
ylabel('velocidad(m/s)');


