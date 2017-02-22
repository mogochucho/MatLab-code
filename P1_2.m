clear all;
close all;
clc;

%Carga de datos
perimetros = load('Perimetro.txt');
R = perimetros(:,1); %radios originales en cm
Se = perimetros(:,2)*10; %perimetros originales en dm

areas = load('Areas.txt');
Ae = areas(:,2)/100; %areas en mm^2

%Valores teorcos
S = 2*pi*R;
A = pi*R.^2;

figure;
plotyy(R, S, R, A);

%Comparamos perimetros teoricos y experimentales
figure;
plot(R,S, R,Se, 'r.');

%Graficamos con error 0.1
figure;
hold on;
errorbar(R, Se, [.1;.1;.1;.1],'r.'); %el resultado teorico no esta en el intervalo
plot(R,S);
hold off;

%Comparamos areas teoricos y experimentales
figure;
plot(R,A, R,Ae, 'r.');

%Graficamos con error 
figure;
hold on;
errorbar(R, Ae, Ae*.1,'r.'); %el resultado teorico no esta en el intervalo
plot(R,A);
hold off;

%Comparamos valores teoricos y experimentales en una misma ventana
figure;
subplot(2,1,1);
plot(R,S, R,Se, 'r.');
subplot(2,1,2);
plot(R,A, R,Ae, 'r.');
