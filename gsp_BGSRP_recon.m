function y=gsp_BGSRP_recon(G,x0,y0,param)
%   Reconstruction of graph signal for BGSRP 
%   Usage:  y = gsp_BGSRP_recon(G,x0,y0,param)
%            
%   Input parameters :
%         G          : Graph structure in gspbox.
%         x0         : Locations of known labeled points
%         y0         : Labels of known labeled points
%         param      : Optional parameters
%   Output parameters:
%         y          : Reconstructed graph signal
%
%   'gsp_BGSRP_recon(G,x0,y0,param)' computes an estimation of the first 
%   k eigenvectors of the Laplacian of G using Gaussian random signal
%   filtering, following the FEARS method described in paratte2017fast.
%
%
%   Example:
%         % Generate a graph and calculate the eigendecomposition.
%         G = gsp_sensor(500);                
%         G = gsp_compute_fourier_basis(G);
%
%         % Construct an original signal with the original bandwidth N/4
%         x=2*rand(G.N,1)-1;
%         fx=G.U'*x;
%         f=G.U(:,1:fix(G.N/10))*fx(1:fix(G.N/10));
%         
%         % Locations and labels of known labeled points, parameters
%         p = randperm (G.N);
%         labs=fix(0.2*G.N);
%         x0=p(1:labs);
%         y0=f(x0);
%         param.bandw=fix(G.N/5);
%         param.gamma=0.1;
%         param.basis = 'lp-ture';

%         % Reconstruction and mean absolute error.
%         y = gsp_BGSRP_recon(G,x0,y0,param);
%         err = mean(abs(f-y));
%          
%       
%
%   Additional parameters
%   ---------------------
%    param.bandw  : The bandwidth coeffcient n in the algorithm. 
%    param.gamma  : The regularization parameter used to balance the weight 
%                   between the error and the quality index
%    param.basis  : Select the Fourier basis to be used for the computation. 
%      'lp-ture'   : Use the accurate Fourier basis of graph Laplacian
%      'lp-FEARS'  : Use the approximate Fourier basis by FEARS method
%       
% 
%   References:
%     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
%     Algorithm for Recovery of Bandlimited Graph Signals Based on the 
%     Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
%     Processing, revised vesion, T-SP-25570-2019, 2020.2.     
%
% Copyright (C) 2020 Qian Zhang and Lihua Yang. All rights reserved.
%
% This function is free software for scientific research, you can 
% redistribute it and/or modify it. But not for commercial use.
%
% This function is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.   

% This function is suggested to use in the context of the GSPBOX toolbox,
% when you use this algorithm, please kindly cite
%     N. Perraudin, J. Paratte, D. Shuman, V. Kalofolias, P. Vandergheynst,
%     and D. K. Hammond. GSPBOX: A toolbox for signal processing on graphs.
%     ArXiv e-prints, Aug. 2014.
%     and 
%     Qian Zhang, Chao Huang, Zhihua Yang, and Lihua Yang. A Fast 
%     Algorithm for Recovery of Bandlimited Graph Signals Based on the 
%     Reproducing Kernel Hilbert Space, IEEE Transactions on Signal 
%     Processing, 2020. 
%     and if needed
%     J. Paratte and L. Martin. Fast eigenspace approximation using random
%     signals. arXiv preprint arXiv:1611.00938, 2016.

% Author: Qian Zhang and Lihua Yang
% Date:   10 February 2020
% Version:2020.2.10

ell=numel(x0);
n=param.bandw;

if ~isfield(param, 'basis'), param.basis = 'lp-ture'; end

switch param.basis
    case 'lp-ture'
    % the eigendecomposition of graph Laplacian is treated as known message
        Un=G.U(:,2:n);
        mu=G.e;
        
    case 'lp-FEARS'
    % the eigendecomposition of graph Laplacian is approximeted by FEARS
        basis = gsp_eigenspace_estimation(G,n);
        mun=zeros(n,1);
        for i=n:-1:2
            mun(i)=basis(:,i)'*G.L*basis(:,i);
        end
        Un=basis(:,2:n);
        mu=mun;
        
    otherwise
        error('The Fourier basis is missing!');
end

Phin=diag(1./mu(2:n));
A=eye(ell)-ones(ell)/ell;
B=Un(x0,:);

GG=B*(param.gamma*Phin+Phin*B'*A*B*Phin)*B';
d=B*Phin*B'*A*y0;

xi=pinv(GG)*d;
g=Un*Phin*B'*xi;
z=-ones(1,ell)*(g(x0)-y0)*sqrt(G.N)/ell;
y=z*ones(G.N,1)/sqrt(G.N)+g;