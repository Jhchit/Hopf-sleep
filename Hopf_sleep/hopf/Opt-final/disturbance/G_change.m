%改变G
clc;
clear;


load('Hopf_sleep\opt-genetic\Opt-final\Simulaciones\Test_100\corridatotal.mat','out','Cfg','SC_cell','Hopf');

RtG = out.RtG;
RtG(77) = 0;
RtG(78) = 0;

new_a = out.Rta;
new_G = RtG;
SC_xs = SC_cell{1};
TRsec = Cfg.filt.TRsec;
nNodes = Cfg.nNodes;
long_Total = Hopf.long_Total(1);
w = Hopf.w{1};
val = Hopf.val{1};
Gmethod = Hopf.Gmethod;
xs = resim_Hopf(new_a,new_G,SC_xs,TRsec,nNodes,long_Total,w,val,Gmethod);

meanfc = observable_FC(xs);