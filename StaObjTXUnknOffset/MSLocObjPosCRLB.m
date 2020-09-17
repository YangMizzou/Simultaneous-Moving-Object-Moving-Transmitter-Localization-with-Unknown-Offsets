function [CRLB_Jnt,CRLB_Ind,CRLB_NoOfst]=MSLocObjPosCRLB(u,t,s,Qr,Qd)
% [CRLB_Jnt,CRLB_Ind,CRLB_NoOfst]=MSLocObjPosCRLB(u,t,s,Qr,Qd)
%
% This function obtains the CRLBs for the object position and velocity for 
% the cases of (i) using both indirect- and direct-path measurements,
% (ii) using indirect-path measurement only, and (iii) without time and 
% frequency offsets.
%
% Input parameter list:
% u:         (Dim x 1), object position.
% t:         (Dim x 1), transmitter position.
% s:         (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% Qd:        (M x M), covariance matrix of the direct range measurements.
% 
% Output parameter list:
% CRLB_Jnt:  (2 Dim x 2 Dim), CRLB of object position for joint estimation 
%            with transmitter position and offsets.
% CRLB_Ind:  (2 Dim x 2 Dim), CRLB of object position without estimation 
%            of transmitter location.
% CRLB_NoOfst: (2 Dim x 2 Dim), CRLB of object position for joint estimation 
%              with transmitter position without time and frequency
%              offsets.
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


[K,M]=size(s);
H=[-ones(M-1,1),eye(M-1)]';

rho_dt=[];
rho_rt=[];
rho_ru=[];

for i=1:M
    rho_dti=(t-s(:,i))/norm(t-s(:,i));
    rho_dt=[rho_dt;rho_dti'];
    rho_rti=(t-u)/norm(t-u);
    rho_rt=[rho_rt;rho_rti'];
    rho_rui=(u-s(:,i))/norm(u-s(:,i))+(u-t)/norm(u-t);
    rho_ru=[rho_ru;rho_rui'];
end

rho_mrphi=[rho_rt,ones(M,1)];
rho_mdphi=[rho_dt,ones(M,1)];

K_MrMd=Qr+rho_mrphi*inv(rho_mdphi'*inv(Qd)*rho_mdphi)*rho_mrphi';
K_MrMdn=Qr+rho_rt*inv(rho_dt'*inv(Qd)*rho_dt)*rho_rt';
K_Mr=H*inv(H'*Qr*H)*H';

FIM=rho_ru'*inv(K_MrMd)*rho_ru;
FIMn=rho_ru'*inv(K_MrMdn)*rho_ru;
FIMq=rho_ru'*K_Mr*rho_ru;

CRLB_Ind=inv(FIM);
CRLB_NoOfst=inv(FIMn);
CRLB_Jnt=inv(FIMq);
