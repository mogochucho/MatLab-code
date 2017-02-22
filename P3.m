clear all;
close all;
clc;
%Carga de datos
boyle_V_P = load('boyle.txt');
boyle1_V_P = load('boyle1.txt');
boyle3_V_P = load('boyle3.txt');

boyle_r1 = corrcoef(boyle1_V_P(:,1).^(-1), boyle1_V_P(:,2));
boyle_r = corrcoef(boyle_V_P(:,1).^(-1), boyle_V_P(:,2));
boyle_r3 = corrcoef(boyle3_V_P(:,1).^(-1), boyle3_V_P(:,2));

fprintf('Los 3 coeficientes de correlación para los datos de Boyle son: %f, %f y %f\n', boyle_r1(2,1), boyle_r(2,1), boyle_r3(2,1));

boyle1_V_P(:,1) = boyle1_V_P(:,1) + .4; %Sumandole el volumen de la válvula.
V = [min(boyle_V_P(:,1)):0.1:max(boyle_V_P(:,1))];

boyle_a_b = polyfit(boyle_V_P(:,1).^(-1), boyle_V_P(:,2),1); %Valores de los coeficientes del mejor ajuste por minimos cuadrados P(1/V)
figure;
plot(boyle_V_P(:,1).^(-1), boyle_V_P(:,2),'.', boyle_V_P(:,1).^(-1), boyle_a_b(1).*boyle_V_P(:,1).^(-1)+boyle_a_b(2)); %Grafica P(1/V) del ajuste y de los datos experimentales
figure;
plot(boyle_V_P(:,1), boyle_V_P(:,2),'.', V, boyle_a_b(1).*V.^(-1)+boyle_a_b(2)); %Grafica P(V) del ajuste y de los datos experimentales

boyle_ins_a = abs(boyle_a_b(1))*sqrt((1/boyle_r(1,2)^2)-1)/sqrt(length(boyle_V_P(:,1))-2);
if(boyle_ins_a < 2.5*10^floor(log10(boyle_ins_a)))
    fprintf( 'a = (%.1f \\pm %.1f) \\times 10^%d\n', round(boyle_a_b(1)/10^floor(log10(boyle_ins_a)),1), round(boyle_ins_a/10^floor(log10(boyle_ins_a)),1), floor(log10(boyle_ins_a)));
else
    fprintf( 'a = (%.0f \\pm %.0f) \\times 10^%d\n', round(boyle_a_b(1)/10^floor(log10(boyle_ins_a))), round(boyle_ins_a/10^floor(log10(boyle_ins_a))), floor(log10(boyle_ins_a)));
end
boyle_ins_n = sqrt((boyle_ins_a*10^(-3)/(8.314*296.15))^2+(boyle_a_b(1)*10^(-3)/(8.314*296.15^2))^2);
if(boyle_ins_n < 2.5*10^floor(log10(boyle_ins_n)))
    fprintf( 'n = (%.1f \\pm %.1f) \\times 10^%d\n', round((boyle_a_b(1)*10^(-3)/(8.314*296.15))/10^floor(log10(boyle_ins_n)),1), round(boyle_ins_n/10^floor(log10(boyle_ins_n)),1), floor(log10(boyle_ins_n)));
else
    fprintf( 'n = (%.0f \\pm %.0f) \\times 10^%d\n', round((boyle_a_b(1)*10^(-3)/(8.314*296.15))/10^floor(log10(boyle_ins_n))), round(boyle_ins_n/10^floor(log10(boyle_ins_n))), floor(log10(boyle_ins_n)));
end

%Carga de datos
gay_lussac_T_R_P = load('gay_lussac.txt');

gay_lussac_a_b = polyfit(gay_lussac_T_R_P(:,1), gay_lussac_T_R_P(:,3),1); %Valores de los coeficientes del mejor ajuste por minimos cuadrados P(T)
T = [-gay_lussac_a_b(2)/gay_lussac_a_b(1) max(gay_lussac_T_R_P(:,1))];

gay_lussac_r = corrcoef(gay_lussac_T_R_P(:,1), gay_lussac_T_R_P(:,3));
fprintf('El coeficiente de correlación para para los datos de la Gay-Lussac es %f\n', gay_lussac_r(1,2));

figure;
plot(gay_lussac_T_R_P(:,1)+273, gay_lussac_T_R_P(:,3),'.', gay_lussac_T_R_P(:,1)+273, gay_lussac_a_b(1).*gay_lussac_T_R_P(:,1)+gay_lussac_a_b(2)); %Grafica P(T) del ajuste y de los datos experimentales
figure;
plot(T, gay_lussac_a_b(1).*T+gay_lussac_a_b(2)); %Grafica P(T) del ajuste

gay_lussac_r = corrcoef(gay_lussac_T_R_P(:,1), gay_lussac_T_R_P(:,3));
gay_lussac_ins_a = abs(gay_lussac_a_b(1))*sqrt((1/gay_lussac_r(1,2)^2)-1)/sqrt(length(gay_lussac_T_R_P(:,1))-2);
gay_lussac_ins_b = gay_lussac_ins_a*sqrt(sum(gay_lussac_T_R_P(:,1).^2)/length(gay_lussac_T_R_P(:,1)));
gay_lussac_ins_T0 = sqrt((gay_lussac_ins_b/gay_lussac_a_b(1))^2+(gay_lussac_ins_a*gay_lussac_a_b(2)/gay_lussac_a_b(1)^2)^2);
if(gay_lussac_ins_T0 < 2.5*10^floor(log10(gay_lussac_ins_T0)))
    fprintf( 'T(P=0) = (%.1f \\pm %.1f) \\times 10^%d\n', round((-gay_lussac_a_b(2)/gay_lussac_a_b(1))/10^floor(log10(gay_lussac_ins_T0)),1), round(gay_lussac_ins_T0/10^floor(log10(gay_lussac_ins_T0)),1), floor(log10(gay_lussac_ins_T0)));
else
    fprintf( 'T(P=0) = (%.0f \\pm %.0f) \\times 10^%d\n', round((-gay_lussac_a_b(2)/gay_lussac_a_b(1))/10^floor(log10(gay_lussac_ins_T0))), round(gay_lussac_ins_T0/10^floor(log10(gay_lussac_ins_T0))), floor(log10(gay_lussac_ins_T0)));
end
if(gay_lussac_ins_a < 2.5*10^floor(log10(gay_lussac_ins_a)))
    fprintf( 'gay_lussac_a = (%.1f \\pm %.1f) \\times 10^%d\n', round(gay_lussac_a_b(1)/10^floor(log10(gay_lussac_ins_a)),1), round(gay_lussac_ins_a/10^floor(log10(gay_lussac_ins_a)),1), floor(log10(gay_lussac_ins_a)));
else
    fprintf( 'gay_lussac_a = (%.0f \\pm %.0f) \\times 10^%d\n', round(gay_lussac_a_b(1)/10^floor(log10(gay_lussac_ins_a))), round(gay_lussac_ins_a/10^floor(log10(gay_lussac_ins_a))), floor(log10(gay_lussac_ins_a)));
end
if(gay_lussac_ins_b < 2.5*10^floor(log10(gay_lussac_ins_b)))
    fprintf( 'gay_lussac_b = (%.1f \\pm %.1f) \\times 10^%d\n', round(gay_lussac_a_b(2)/10^floor(log10(gay_lussac_ins_b)),1), round(gay_lussac_ins_b/10^floor(log10(gay_lussac_ins_b)),1), floor(log10(gay_lussac_ins_b)));
else
    fprintf( 'gay_lussac_b = (%.0f \\pm %.0f) \\times 10^%d\n', round(gay_lussac_a_b(2)/10^floor(log10(gay_lussac_ins_b))), round(gay_lussac_ins_b/10^floor(log10(gay_lussac_ins_b))), floor(log10(gay_lussac_ins_b)));
end
gay_lussac_ins_n = sqrt((gay_lussac_ins_a*10.4*10^(-3)/8.314)^2+(gay_lussac_a_b(1)*10^(-3)/8.314)^2);
if(gay_lussac_ins_n < 2.5*10^floor(log10(gay_lussac_ins_n)))
    fprintf( 'gay_lussac_n = (%.1f \\pm %.1f) \\times 10^%d\n', round((gay_lussac_a_b(1)*10.4*10^(-3)/8.314)/10^floor(log10(gay_lussac_ins_n)),1), round(gay_lussac_ins_n/10^floor(log10(gay_lussac_ins_n)),1), floor(log10(gay_lussac_ins_n)));
else
    fprintf( 'gay_lussac_n = (%.0f \\pm %.0f) \\times 10^%d\n', round((gay_lussac_a_b(1)*10.4*10^(-3)/8.314)/10^floor(log10(gay_lussac_ins_n))), round(gay_lussac_ins_n/10^floor(log10(gay_lussac_ins_n))), floor(log10(gay_lussac_ins_n)));
end

calib_T = gay_lussac_T_R_P(:,1) + 273.15;
calib_x = calib_T.^(-1) - calib_T(1).^(-1);
calib_y = log(gay_lussac_T_R_P(:,2));
calib_a_b = polyfit(calib_x, calib_y,1); %Valores de los coeficientes del mejor ajuste por minimos cuadrados P(1/V)
figure;
plot(calib_x, calib_y, '.', calib_x, calib_a_b(1)*calib_x + calib_a_b(2));

calib_r = corrcoef(calib_x, calib_y);
fprintf('El coeficiente de correlación para los datos de la calibración es %f\n', calib_r(1,2));

calib_ins_a = abs(calib_a_b(1))*sqrt((1/calib_r(1,2)^2)-1)/sqrt(length(calib_x)-2);
calib_ins_b = calib_ins_a*sqrt(sum(calib_x.^2)/length(calib_x));
if(calib_ins_a < 2.5*10^floor(log10(calib_ins_a)))
    fprintf( 'B = (%.1f \\pm %.1f) \\times 10^%d\n', round(calib_a_b(1)/10^floor(log10(calib_ins_a)),1), round(calib_ins_a/10^floor(log10(calib_ins_a)),1), floor(log10(calib_ins_a)));
else
    fprintf( 'B = (%.0f \\pm %.0f) \\times 10^%d\n', round(calib_a_b(1)/10^floor(log10(calib_ins_a))), round(calib_ins_a/10^floor(log10(calib_ins_a))), floor(log10(calib_ins_a)));
end

