function [sol]=MSLocObjPosInd_MLE(r,s,Qr,initGuess);
% [sol]=MSLocObjTxPosInd_MLE(r,s,Qr,initGuess);
%
% This function realizes the iterative Gauss-Newton algorithm for
% estimating the unknown object position using indirect-path range 
% measurements.
%
% Input parameter list:
% r:         (M x 1), indirect-path range measurements.
% s:         (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% initGuess: (Dim x 1), initial solution guess.
% 
% Output parameter list:
% sol:       (Dim x 1), final solution of ObjPos.
%
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

[K,M]=size(s);    % M=number of receivers
                  % K=dimension

u=initGuess;

imax=20;
H=[-ones(M-1,1),eye(M-1)]';

for j=1:1:imax
    ro_hat=sqrt(sum((repmat(u,1,M-1)-s(:,2:end)).^2))'-norm(u-s(:,1));
    rho_us=repmat(u,1,M-1)-s(:,2:end); rho_us=rho_us./(ones(K,1)*sqrt(sum(rho_us.^2)));
    rho_uso=(u-s(:,1))/norm(u-s(:,1));
    Gr=(rho_us-rho_uso*ones(1,M-1))';
    G=Gr;
    e=((r(2:end)-r(1))-ro_hat);
    u=u+0.5*inv(G'*inv(H'*Qr*H)*G)*G'*inv(H'*Qr*H)*e;
end

sol=u;

