clc;
close all;
clear all;

datos = load('Cronometro_fotosensor.txt');
VTcron = datos(:,1)/20; %Vector de per�odos tomados con el cron�metro.
VTfoto = datos(:,2)/20; %Vector de per�odos tomados con el fotosensor.
Vgcron = .58 * (2*pi*VTcron.^(-1)).^2; %Vector de la gravedad calculada a partir de los datos del cron�metro.
Vgfoto = .58 * (2*pi*VTfoto.^(-1)).^2; %Vector de la gravedad calculada a partir de los datos del fotosensor.
Tcron = mean(VTcron); %Mejor estimaci�n del per�odo basada en los datos del cron�metro.
Tfoto = mean(VTfoto); %Mejor estimaci�n del per�odo basada en los datos del fotosensor.
Tcron = mean(VTcron); %Mejor estimaci�n del per�odo basada en los datos del cron�metro.
Tfoto = mean(VTfoto); %Mejor estimaci�n del per�odo basada en los datos del fotosensor.
gcron = mean(Vgcron); %Mejor estimaci�n de la gravedad basada en los datos del cron�metro.
gfoto = mean(Vgfoto); %Mejor estimaci�n de la gravedad basada en los datos del fotosensor.
uTcron = std(VTcron); %Desviaci�n estandar de los valores tomados con el cron�metro.
uTfoto = std(VTfoto); %Desviaci�n estandar de los valores tomados con el fotosensor.
ugcron = std(Vgcron); %Desviaci�n estandar de los valores calculados de la gravedad basados en el cron�metro.
ugfoto = std(Vgfoto); %Desviaci�n estandar de los valores calculados de la gravedad basados en el fotosensor.

v = zeros(size(datos));
v(:,1) = Vgcron;
v(:,2) = Vgfoto;

% Para quitar todas las medidas que est�n a m�s de 2 veces la desviaci�n estandar.
for i = 1:length(VTcron)
    if abs(VTcron(i) - Tcron) > 2.5 * uTcron
        VTcron(i) = 0;
    end
    if abs(Vgcron(i) - gcron) > 2.5 * ugcron
        Vgcron(i) = 0;
    end
    if abs(VTfoto(i) - Tfoto) > 2.5 * uTfoto
        VTfoto(i) = 0;
    end
    if abs(Vgfoto(i) - gfoto) > 2.5 * ugfoto
        Vgfoto(i) = 0;
    end
end
VTcron(VTcron==0) = [];
VTfoto(VTfoto==0) = [];
Vgcron(Vgcron==0) = [];
Vgfoto(Vgfoto==0) = [];

% Recalculando los estimadores y la desviaci�n estandar.
Tcron = mean(VTcron);
gcron = mean(Vgcron);
Tfoto = mean(VTfoto);
gfoto = mean(Vgfoto);
uTcron = std(VTcron);
ugcron = std(Vgcron);
uTfoto = std(VTfoto);
ugfoto = std(Vgfoto);

x = [Tcron - 3*uTcron:.001:Tcron + 3*uTcron];
y = (length(VTcron) * (max(VTcron)-min(VTcron))/20) * exp(-(x - Tcron).^2/(2*uTcron^2))/sqrt(2*pi*uTcron^2);
%Valores de la campana de gauss multiplicados por el �rea total del histograma.

figure;
hold on;
hist(VTcron, 20);
plot(x,y,'r');
hold off;

x = [Tfoto - 3*uTfoto:.01:Tfoto + 3*uTfoto];
y = (length(VTfoto) * (max(VTfoto)-min(VTfoto))/20) * exp(-(x - Tfoto).^2/(2*uTfoto^2))/sqrt(2*pi*uTfoto^2);
%Valores de la campana de gauss multiplicados por el �rea total del histograma.

figure;
hold on;
hist(VTfoto, 20);
plot(x,y,'r');
hold off;

x = [gcron - 3*ugcron:.01:gcron + 3*ugcron];
y = (length(Vgcron) * (max(Vgcron)-min(Vgcron))/20) * exp(-(x - gcron).^2/(2*ugcron^2))/sqrt(2*pi*ugcron^2);
%Valores de la campana de gauss multiplicados por el �rea total del histograma.

figure;
hold on;
hist(Vgcron, 20);
plot(x,y,'r');
hold off;

x = [gfoto - 3*ugfoto:.01:gfoto + 3*ugfoto];
y = (length(Vgfoto) * (max(Vgfoto)-min(Vgfoto))/20) * exp(-(x - gfoto).^2/(2*ugfoto^2))/sqrt(2*pi*ugfoto^2);
%Valores de la campana de gauss multiplicados por el �rea total del histograma.

figure;
hold on;
hist(Vgfoto, 20);
plot(x,y,'r');
hold off;

'Medida de g seg�n la aplicaci�n de m�todos estadisticos a el vector de g`s calculado:'

gcron
ugcron
gfoto
ugfoto

'Medida de g seg�n la aplicaci�n de m�todos estadisticos a las medidas originales:'

gcron = .58 * (2*pi/Tcron)^2 % Se muestra la estimaci�n de g seg�n las medidas realizadas con el cron�metro
ugcron = sqrt((4*pi^2*.01/Tcron^2)^2+(.58*8*pi^2*uTcron/Tcron^3)^2) %y su correspondiente incertidumbre.

gfoto = .58 * (2*pi/Tfoto)^2 % Se muestra la estimaci�n de g seg�n las medidas realizadas con el fotosensor
ugfoto = sqrt((4*pi^2*.01/Tfoto^2)^2+(.58*8*pi^2*uTfoto/Tfoto^3)^2) %y su correspondiente incertidumbre.

