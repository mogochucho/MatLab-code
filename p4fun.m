function p4fun( medio , mlocal, umlocal )

fprintf(['Todos los cálculos y medidas son para el caso del ' medio '\n\n']);
pause;

opt = input('Elija una opción:\n1- Utilizar la envolvente guardada.\n2- Construir una nueva envolvente.\n- ');
if (opt == 1)
    datos = load(['datos' medio '.txt']);
    t_x = load(['envolvente' medio '.txt']);
    r = corrcoef(t_x(1,:),log(t_x(2,:)));
    r_best = r(1,2);
elseif (opt == 2)
    datos = cat(3, load(['datos1' medio '.txt']), load(['datos2' medio '.txt']), load(['datos3' medio '.txt']));
    eq = load(['distancia_equilibrio_' medio '.txt']);
    eq = mean(eq(:,2));
    t_x = zeros(2,20);
    r = zeros(2,2);
    for i = 1:3
        datos(:,2,i) = datos(:,2,i) - eq;
        figure;
        plot(datos(:,1,i),datos(:,2,i));
        [t_x(1,:,i), t_x(2,:,i)] = ginput(20);
        r(:,:,i) = corrcoef(t_x(1,:,i),log(t_x(2,:,i)));
        fprintf('El coeficiente de correlación %d es %f\n', i, r(1,2,i));
    end
    opt = input('Elija una opción:\n1- Comparar con los datos guardados.\n2- Guardar el mejor.\n3- Descartar.\n- ');
    if (opt == 3)
        t_x = load(['envolvente' medio '.txt']);
        datos = load(['datos' medio '.txt']);
        r = corrcoef(t_x(1,:),log(t_x(2,:)));
        r_best = r(1,2);
        fprintf('Se utilizaran los datos guardados\n');
    else
        r_best = r(1,2,1);
        j = 1;
        for i = 2:3
            if abs(r(1,2,i)) > abs(r_best)
                r_best = r(1,2,i);
                j = i;
            end
        end
        datos = datos(:,:,j);
        t_x = t_x(:,:,j);
        if (opt == 2)
            save(['datos' medio '.txt'], 'datos', '-ascii');
            save(['envolvente' medio '.txt'], 't_x', '-ascii');
            fprintf('Se guardaron los nuevos datos\n');
        elseif (opt == 1)
            t_x_old = load(['envolvente' medio '.txt']);
            r = corrcoef(t_x_old(1,:),log(t_x_old(2,:)));
            if abs(r_best) > abs(r(1,2))
                save(['datos' medio '.txt'], 'datos', '-ascii');
                save(['envolvente' medio '.txt'], 't_x', '-ascii');
                fprintf('Se actualizaron los datos\n');
            else
                datos = load(['datos' medio '.txt']);
                t_x = t_x_old;
                fprintf('El nuevo ajuste es peor, así que no se actualizan los datos\n');
                r_best = r(1,2);
            end
        end
    end
end

opt = input('Elija una opción:\n1- Utilizar la período guardado.\n2- Determinar un nuevo período.\n- ');
if opt == 1
    Ts = load(['Ts' medio '.txt']);
    uT = std(Ts);
elseif opt == 2
    fprintf('Pinche en 21 picos seguidos para determinar el período\n');
    figure;
    plot(datos(:,1),datos(:,2));
    T21 = zeros(2,21);
    [T21(1,:), T21(2,:)] = ginput(21);
    Ts = zeros(1,20);
    for i = 1:20
        Ts(i) = T21(1,i+1) - T21(1,i);
    end
    opt = input('Elija una opción:\n1- Guardar si son mejores que los datos guardados.\n2- Guardar datos.\n3- Descartar.\n- ');
    if opt == 3
        Ts = load(['Ts' medio '.txt']);
        uT = std(Tsold);
    else
        uT = std(Ts);
        if opt == 2
            save(['Ts' medio '.txt'], 'Ts', '-ascii');
            fprintf('Los datos fueron guardados.\n');
        elseif opt == 1
            Tsold = load(['Ts' medio '.txt']);
            uTold = std(Tsold);
            if uT < uTold
                save(['Ts' medio '.txt'], 'Ts', '-ascii');
                fprintf('Los nuevos datos son mejores y fueron guardados.\n');
            else
                Ts = Tsold;
                uT = uTold;
                fprintf('Los nuevos datos no son mejores por lo que fueron cargados los anteriores.\n');
            end
        end
    end
end
T = mean(Ts);
showmedida('T', T, uT, 's');
fprintf('r = %f\n', r_best);
pause;

M_N = polyfit(t_x(1,:),log(t_x(2,:)),1);
uM = abs(M_N(1))*sqrt((1/r_best^2)-1)/sqrt(length(t_x(1,:))-2);
showmedida('M', M_N(1), uM, '');
uN = uM*sqrt(sum(t_x(1,:).^2)/length(t_x(1,:)));
showmedida('N', M_N(2), uN, '');

figure;
plot(datos(:,1),M_N(1).*datos(:,1) + M_N(2),'k',t_x(1,:),log(t_x(2,:)),'r.');
grid on
legend('ajuste por mínimos cuadrados','logaritmo de máximos de elongación')
title(['LOGARITMO DE LA ENVOLVENTE EN FUNCION DEL TIEMPO. (' medio ')'])
ylabel('logaritmo de la envolvente (L(m))')
xlabel('tiempo (s)')

figure;
plot(datos(:,1), datos(:,2),'r', datos(:,1), exp(M_N(1).*datos(:,1) + M_N(2)), 'k');
grid on
legend('elongación','envolvente')
title(['ELONGACIÓN EN FUNCION DEL TIEMPO. (' medio ')'])
ylabel('elongación (m)')
xlabel('tiempo (s)')

global m k b;

m = mlocal;
um = umlocal;
showmedida('m', m, um, 'kg');
b = -2*m*M_N(1);
ub = 2*sqrt((m*uM)^2+(um*M_N(1))^2);
showmedida('b', b, ub, 'kg/s');
k = 4*pi^2*m*T^(-2)+b^2/(4*m);
uk = sqrt((4*pi^2*um/T^2)^2+(8*pi^2*m*uT/T^3)^2+(b*ub/(2*m))^2+(b^2*um/(4*m^2))^2);
showmedida('k_b', k, uk, 'N/m');
kdespb = 4*pi^2*m*T^(-2);
ukdespb = 4*pi^2*sqrt((um/T^2)^2+(2*m*uT/T^3)^2);
showmedida('k', kdespb, ukdespb, 'N/m');
pause;

y0=[datos(1,2), datos(1,3)];
[t,y]=ode23('amortig2',[datos(1,1) datos(length(datos(:,1)),1)],y0);

figure;
plot(t,y(:,1),'k',datos(:,1),datos(:,2),'r');
grid on
legend('elongación teórica','elongación experimental')
title(['ELONGACION EN FUNCION DEL TIEMPO. (' medio ')'])
ylabel('elongación (m)')
xlabel('tiempo (s)')

figure;
plot(t,y(:,2),'k',datos(:,1),datos(:,3),'r')
grid on
legend('velocidad teórica','velocidad experimental')
title(['VELOCIDAD EN FUNCION DEL TIEMPO. (' medio ')'])
ylabel('velocidad (m/s)')
xlabel('tiempo (s)')

figure;
plot(y(:,1),y(:,2),'k',datos(:,2),datos(:,3),'r')
grid on
legend('velocidad teórica','velocidad experimental')
title(['VELOCIDAD EN FUNCION DE LA ELONGACIÓN. (' medio ')'])
ylabel('velocidad (m/s)')
xlabel('elongación (m)')

figure;
plot(t,k*y(:,1).^2/2+m*y(:,2).^2/2,'k',datos(:,1),k*datos(:,2).^2/2+m*datos(:,3).^2/2,'r')
grid on
legend('energía teórica','energía experimental')
title(['ENERGÍA EN FUNCION DEL TIEMPO. (' medio ')'])
ylabel('energía (J)')
xlabel('tiempo (s)')

fprintf('\n');
fprintf('\n');

end

