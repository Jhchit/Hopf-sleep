clc;
clear;

load('Hopf_sleep\opt-genetic\Opt-final\Simulaciones_n_normal\Test\corridatotal.mat');
FSim_N = zeros(90,90);
 mean_fval = 0;
for i = 1:100
    mean_fval = mean_fval + out(i).fval;
    FSim_N = FSim_N + out(i).FSim;
end
mean_fval =  mean_fval/100;
FSim_N = FSim_N/100;
FEmp_N = FEmp;

load('Hopf_sleep\opt-genetic\Opt-final\Simulaciones_n_unsleep\Test\corridatotal.mat');
FSim_U = zeros(90,90);
mean_fval = 0;
for i = 1:100
    mean_fval = mean_fval + out(i).fval;
    FSim_U = FSim_U + out(i).FSim;
end
mean_fval =  mean_fval/100;
FSim_U = FSim_U/100;

load('Hopf_sleep\opt-genetic\Opt-final\Simulaciones\Test\corridatotal.mat');
FSim_F = zeros(90,90);
mean_fval = 0;
n = irepe;
for i = 1:n
    mean_fval = mean_fval + out(i).fval;
    FSim_F = FSim_F + out(i).FSim;
end
mean_fval =  mean_fval/n;
FSim_F = FSim_F/n;

GoF = (ssim(FSim_N,FEmp_N)-ssim(FSim_F,FEmp_N))/(ssim(FSim_N,FEmp_N)-ssim(FSim_U,FEmp_N));