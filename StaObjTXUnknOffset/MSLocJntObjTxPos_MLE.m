function sol=MSLocJntObjTxPos_MLE(rm,dm,s,Qr,Qd,initGuess)
% sol=MSLocJntObjTxPos_MLE(rm,dm,s,Qr,Qd,initGuess)
%
% This function realizes the iterative Gauss-Newton algorithm for jointly
% estimating the unknown object and transmitter positions using both 
% indirect- and direct-path range measurements.
%
% Input parameter list:
% r:         (M x 1), indirect-path range measurements.
% d:         (M x 1), direct-path range measurements.
% s:         (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% Qd:        (M x M), covariance matrix of the direct range measurements.
% initGuess: (4*Dim+2 x 1), initial solution guess.
% 
% Output parameter list:
% sol:       (2*Dim+1 x 1), final solution, 
%            [objPos;txPos;ofStTm]^T.

% The program can be used for 2D(Dim=2) or 3D(Dim=3) localization.
%
% Reference:
% Y. Zhang and K. C. Ho, "Multistatic moving object localization by a 
% moving transmitter of unknown location and offset," IEEE Trans. Signal 
% Process., vol. 68, pp. 4438-4453, 2020.
% 
% Yang Zhang and K. C. Ho   08-22-2020
%
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

[K,M]=size(s);      % M=number of receivers
                    % K=dimension
sol=initGuess;
imax=20;
Q=blkdiag(Qr,Qd);
iQ=inv(Q);

for i=1:1:imax
    u=sol(1:K);
    t=sol(K+1:2*K);
    offsetd=sol(2*K+1);
    
    do_hat=sqrt(sum((repmat(t,1,M)-s(:,1:M)).^2))'+offsetd;
    ro_hat=sqrt(sum((repmat(u,1,M)-s(:,1:M)).^2))'+norm(u-t)+offsetd;
    
    rho_ru=repmat(u,1,M)-s; rho_ru=rho_ru./(ones(K,1)*sqrt(sum(rho_ru.^2)));
    rho_rt=(u-t)/norm(u-t)*ones(1,M);
    rho_u=rho_ru+rho_rt;
    
    rho_dt=repmat(t,1,M)-s; rho_dt=rho_dt./(ones(K,1)*sqrt(sum(rho_dt.^2)));
    
    Jr=[rho_u' -rho_rt' ones(M,1)];
    Jd=[zeros(M,K) rho_dt' ones(M,1)];
    
    J=[Jr;Jd];
    err=[rm-ro_hat;dm-do_hat];
    temp=sol;
    sol=sol+inv(J'*iQ*J)*J'*iQ*err;
    
    if (norm(sol-temp)<1e-6)
        break;
    end
    
end
