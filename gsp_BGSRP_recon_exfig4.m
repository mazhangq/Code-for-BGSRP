% This experiment is for "Fig.4" in paper:
%   References:
%     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
%     Algorithm for Recovery of Bandlimited Graph Signals Based on the 
%     Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
%     Processing, revised vesion, T-SP-25570-2019, 2020.2.    

% This experiment is conducted on sensor graphs by varying the number of vertices
% The results are depicted by boxplot form
% In this experiment, the known labeled vertices and the bandwidth
% coefficients are fixed as l=400,n=400. The eigendecomposition of graph
% Laplacian is treated as known message in Fig.4.



clear all;
close all;

TIME=[];
N=[800,1200,1600,2000,2400,2800,3200,3600,4000];
numN=length(N);
for i=1:numN
    i
 
gsp_reset_seed(0);
G = gsp_sensor(N(i));
G = gsp_compute_fourier_basis(G);
for j=1:50
% Original signal
x=2*rand(G.N,1)-1;
fx=G.U'*x;
f=G.U(:,1:fix(G.N/8))*fx(1:fix(G.N/8));

paramplot.position = [100,100,600,220];
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.1], [0.01 0.1]);

p = randperm (G.N);
labs=400;
param.bandw=400;
param.gamma=0.1;
param.sigmasquare=8;

 
    x0=p(1:labs);
    y0=f(x0);
    t1=clock;
    f1=gsp_BGSRP_recon(G,x0,y0,param);
    t2=clock;
    TIME(i,j)=etime(t2,t1);
end 

end 
 

figure('color',[1 1 1]);
set(gcf,'unit','normalized','position',[0.1,0.1,0.6,0.8]);
boxplot(TIME')
axis([-inf inf -1 3])
set(gca,'xticklabel',{'800','1200','1600','2000','2400','2800','3200','3600','4000'},'linewidth',1.5);
xlabel('graphs with different vertices','fontsize',15)
ylabel('TIMEs','fontsize',15)
