function [sol]=MSLocObjPosVelInd_MLE(r,fr,s,Qr,Qfr,initGuess)
% [sol]=MSLocJntObjPosVelInd_MLE(r,fr,s,Qr,Qfr,initGuess)
%
% This function realizes the iterative Gauss-Newton algorithm for 
% estimating the unknown object position and velocity using indirect-
% path measurements.
%
% Input parameter list:
% r:         (M x 1), indirect-path range measurements.
% fr:        (M x 1), indirect-path range-rate measurements.  
% s:         (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% Qfr:       (M x M), covariance matrix of the indirect range-rate measurements.
% initGuess: (2*Dim x 1), initial solution guess.
% 
% Output parameter list:
% sol:       (2*Dim x 1), final solution, 
%            [objPos;ObjVel]^T.
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

theta=initGuess;
Qrfr=blkdiag(Qr,Qfr);

imax=20;
H=[-ones(M-1,1),eye(M-1)]';
Hq=[H,zeros(M,M-1);zeros(M,M-1),H];
for j=1:1:imax
    uini=theta(1:K);
    vuini=theta(K+1:2*K);
    r_hat=sqrt(sum((repmat(uini,1,M-1)-s(:,2:end)).^2))'-norm(uini-s(:,1));
    fr_hat=(repmat(uini,1,M-1)-s(:,2:end))'*vuini./sqrt(sum((repmat(uini,1,M-1)-s(:,2:end)).^2))'-((uini-s(:,1))'*vuini)/norm(uini-s(:,1));
    rho_us=repmat(uini,1,M-1)-s(:,2:end); rho_us=rho_us./(ones(K,1)*sqrt(sum(rho_us.^2)));
    rho_uso=(uini-s(:,1))/norm(uini-s(:,1));
    rho_vus=vuini./sqrt(sum((repmat(uini,1,M-1)-s(:,2:end)).^2))-(((repmat(uini,1,M-1)-s(:,2:end))'*vuini)'.*(repmat(uini,1,M-1)-s(:,2:end)))./(sqrt(sum((repmat(uini,1,M-1)-s(:,2:end)).^2)).^3);
    rho_vuso=vuini/norm(uini-s(:,1))-((uini-s(:,1))'*vuini*(uini-s(:,1)))/(norm(uini-s(:,1))^3);
    Gr=[(rho_us-rho_uso*ones(1,M-1))',zeros(M-1,K)];
    Gfr=[(rho_vus-rho_vuso*ones(1,M-1))',(rho_us-rho_uso*ones(1,M-1))'];
    G=[Gr;Gfr];
    e=[((r(2:end)-r(1))-r_hat);((fr(2:end)-fr(1))-fr_hat)];
    theta=theta+inv(G'*inv(Hq'*Qrfr*Hq)*G)*G'*inv(Hq'*Qrfr*Hq)*e;
end

sol=theta;

