clear all;
close all;
clc

M12 = load('coeficiente_estatico.txt');
r = corrcoef(M12(:,1),M12(:,2));
r = r(1,2);
a_b = polyfit(M12(:,1),M12(:,2),1);
ua = abs(a_b(1))*sqrt(r^(-2) - 1)/sqrt(length(M12(:,1))-2);
ub = ua * sqrt(sum(M12(:,1))/length(M12(:,1)));

showmedida('a', a_b(1), ua, '');
showmedida('b', a_b(2), ub, 'g');

figure;
plot(M12(:,1), a_b(1)*M12(:,1) + a_b(2), 'k', M12(:,1),M12(:,2), '.r');
grid on
legend('masa teórica','masa experimental')
title('Masa 2 en función de la masa 1');
ylabel('masa 2 (g)');
xlabel('masa 1 (g)');

mus = M12(:,2)./M12(:,1);
mu = mean(mus);
umu = std(mus);
showmedida('mu', mu, umu, '');

