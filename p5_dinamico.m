clear all;
close all;
clc

cant = 50; % Cantidad de archivos
a_n = zeros(50, 2);
for i = 1:cant
    datos = load(['practica_final2/Datos_' int2str(i) '.txt']);
    %datos = load(['practica_final/Datos_' int2str(i) '.txt']);
    datos_t_p_v = [0 0 0];
    aux_v = 0;
    j = 1;
    h = round(length(datos(:,1))/10);
    while j+h < round(length(datos(:,1))*4/5) && not(isnan(datos(j+h,4)))% && datos(j+1,4) > aux_v
        datos_t_p_v(j, 1) = datos(j+1,1);
        datos_t_p_v(j, 2) = datos(j+1,3);
        datos_t_p_v(j, 3) = datos(j+1,4);
        j = j + 1;
        aux_v = datos(j,4);
    end
    A_B_C = polyfit(datos_t_p_v(:, 1), datos_t_p_v(:, 2), 2);
    M_N = polyfit(datos_t_p_v(:, 1), datos_t_p_v(:, 3), 1);
    a_n(i,1) = A_B_C(1)*2;
    a_n(i,2) = M_N(1);
end

a_n_1 = [0];
a_n_2 = [0];
ua = [std(a_n(:,1)), std(a_n(:,2))];
a = [mean(a_n(:,1)), mean(a_n(:,2))];
j = 1;
k = 1;
for i = 1:length(a_n(:,1))
    if abs(a_n(i,1) - a(1)) < 2.5*ua(1)
        a_n_1(j) = a_n(i,1);
        j = j + 1;
    end
    if abs(a_n(i,2) - a(2)) < 2.5*ua(2)
        a_n_2(k) = a_n(i,2);
        k = k + 1;
    end
end
ua = [std(a_n_1), std(a_n_2)];
a = [mean(a_n_1), mean(a_n_2)];

x = a(1)-2.5*ua(1):ua(1)/20:a(1)+2.5*ua(1);
y = (length(a_n_1)*(max(a_n_1) - min(a_n_1))/40)*exp(-(x-a(1)).^2/(2*ua(1)^2))/sqrt(2*pi*ua(1)^2);

figure;
hold on;
hist(a_n_1,40);
plot(x,y,'k');
hold off;
grid on
title({'Histograma de aceleraciones según la posición en función del tiempo';'para la lija como superficie de rozamiento'});
%title({'Histograma de aceleraciones según la posición en función del tiempo';'para la madera como superficie de rozamiento'});
ylabel('instancias');
xlabel('aceleración (m/s²)');

x = a(2)-2.5*ua(2):ua(2)/20:a(2)+2.5*ua(2);
y = (length(a_n_2)*(max(a_n_2) - min(a_n_2))/40)*exp(-(x-a(2)).^2/(2*ua(2)^2))/sqrt(2*pi*ua(2)^2);

figure; 
hold on;
hist(a_n_2,40);
plot(x,y,'k');
hold off;
grid on
title({'Histograma de aceleraciones según la velocidad en función del tiempo';'para la lija como superficie de rozamiento'});
%title({'Histograma de aceleraciones según la velocidad en función del tiempo';'para la madera como superficie de rozamiento'});
ylabel('instancias');
xlabel('aceleración (m/s²)');

m1 = .23345;
m2 = .17864;
um = .00001;
g = 9.8;
ug = 0.1;
mu = (m2 - a.*(m2+m1)/g)/m1;
umu_a = ((m2+m1).*ua/(g*m1)).^2;
umu_g = ((m2+m1).*a*ug/(g^2*m1)).^2;
umu_m1 = (- a/g*m1 + (- m2 + a.*(m2+m1)/g)/m1^2).^2 .* um^2;
umu_m2 = (1/m1 - a/g*m1).^2 .* um^2;
umu = sqrt(umu_a + umu_g + umu_m1 + umu_m2);

showmedida('mu1', mu(1), umu(1), '');
showmedida('mu2', mu(2), umu(2), '');
