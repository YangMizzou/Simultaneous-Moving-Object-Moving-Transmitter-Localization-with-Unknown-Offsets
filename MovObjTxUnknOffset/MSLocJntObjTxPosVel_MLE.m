function [sol]=MSLocJntObjTxPosVel_MLE(r,d,fr,fd,s,Qr,Qd,Qfr,Qfd,initGuess)
% [sol]=MSLocJntObjTxPosVel_MLE(r,d,fr,fd,s,Qr,Qd,Qfr,Qfd,initGuess)
%
% This function realizes the iterative Gauss-Newton algorithm for jointly
% estimating the unknown object and transmitter positions and velocities
% using both indirect- and direct-path range and range-rate measurements.
%
% Input parameter list:
% r:         (M x 1), indirect-path range measurements.
% d:         (M x 1), direct-path range measurements.
% fr:        (M x 1), indirect-path range-rate measurements.
% fd:        (M x 1), direct-path range-rate measurements.
% s:         (Dim x M), receiver position matrix, M is the number of receivers.
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% Qd:        (M x M), covariance matrix of the direct range measurements.
% Qfr:       (M x M), covariance matrix of the indirect range-rate measurements.
% Qfd:       (M x M), covariance matrix of the direct range-rate measurements.
% initGuess: (4*Dim+2 x 1), initial solution guess.
%
% Output parameter list:
% sol:       (4*Dim+2 x 1), final solution,
%            [objPos;ObjVel;txPos;txVel;ofStTm;ofStFq]^T.

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

[K,M]=size(s);      % K=dimension
                    % M=number of receivers

sol=initGuess;
imax=20;
Q=blkdiag(Qr,Qfr,Qd,Qfd);
iQ=inv(Q);

for i=1:1:imax
    u=sol(1:K);
    vu=sol(K+1:2*K);
    t=sol(2*K+1:3*K);
    vt=sol(3*K+1:4*K);
    offsetd=sol(4*K+1);
    offsetf=sol(4*K+2);
    
    do_hat=sqrt(sum((repmat(t,1,M)-s(:,1:M)).^2))'+offsetd;
    ro_hat=sqrt(sum((repmat(u,1,M)-s(:,1:M)).^2))'+norm(u-t)+offsetd;
    fdo_hat=(repmat(t,1,M)-s(:,1:M))'*vt./sqrt(sum((repmat(t,1,M)-s(:,1:M)).^2))'+offsetf;
    fro_hat=(repmat(u,1,M)-s(:,1:M))'*vu./sqrt(sum((repmat(u,1,M)-s(:,1:M)).^2))'+(u-t)'*(vu-vt)/norm(u-t)+offsetf;
    
    rho_ru=repmat(u,1,M)-s; rho_ru=rho_ru./(ones(K,1)*sqrt(sum(rho_ru.^2)));
    rho_rt=(u-t)/norm(u-t)*ones(1,M);
    rho_u=rho_ru+rho_rt;
    
    rho_dt=repmat(t,1,M)-s; rho_dt=rho_dt./(ones(K,1)*sqrt(sum(rho_dt.^2)));
    
    rho_fru=vu./sqrt(sum((repmat(u,1,M)-s).^2))-(((repmat(u,1,M)-s)'*vu)'.*(repmat(u,1,M)-s))./(sqrt(sum((repmat(u,1,M)-s).^2)).^3)+(vu-vt)/norm(u-t)-((u-t)'*(vu-vt)*(u-t))/(norm(u-t)^3);
    rho_frt=(vt-vu)/norm(u-t)-((t-u)'*(vt-vu)*(t-u))/(norm(u-t)^3)*ones(1,M);
    
    rho_fdt=vt./sqrt(sum((repmat(t,1,M)-s).^2))-(((repmat(t,1,M)-s)'*vt)'.*(repmat(t,1,M)-s))./(sqrt(sum((repmat(t,1,M)-s).^2)).^3);
    
    Jr=[rho_u' zeros(M,K) -rho_rt' zeros(M,K) ones(M,1) zeros(M,1)];
    Jd=[zeros(M,K) zeros(M,K) rho_dt' zeros(M,K) ones(M,1) zeros(M,1)];
    Jfr=[rho_fru' rho_u' rho_frt' -rho_rt' zeros(M,1) ones(M,1)];
    Jfd=[zeros(M,K) zeros(M,K) rho_fdt' rho_dt' zeros(M,1) ones(M,1)];
    J=[Jr;Jfr;Jd;Jfd];
    err=[r-ro_hat;fr-fro_hat;d-do_hat;fd-fdo_hat];
    temp=sol;
    sol=sol+inv(J'*iQ*J)*J'*iQ*err;
    
    if (norm(sol-temp)<1e-6)
        break;
    end
end
