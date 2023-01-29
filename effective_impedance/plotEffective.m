
h = subplot(3,1,1);
set(h,'NextPlot','add');
plot(freq,ES)
xline(w_ve,'-',{'\omega_{ve}'})
xline(w_co,'--',{'\omega_{co}'})
set(gca, 'XScale', 'log')
title('Effective stiffness')
h = subplot(3,1,2); 
set(h,'NextPlot','add');
plot(freq,ED)
xline(w_ve,'-')
xline(w_co,'--')
title('Effective damping')
set(gca, 'XScale', 'log')
h = subplot(3,1,3);
plot(freq,EM)
xline(w_ve,'-')
xline(w_co,'--')
set(gca, 'XScale', 'log')
title('Effective mass')
%h = legend('Sys', '$\omega_{\mathrm{ve}}$', '$\omega_{\mathrm{co}}$'); set(h, 'Interpreter', 'latex');

