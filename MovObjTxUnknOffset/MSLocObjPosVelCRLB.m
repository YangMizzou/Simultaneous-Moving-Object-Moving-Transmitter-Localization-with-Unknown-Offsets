function [CRLB_Jnt,CRLB_Ind,CRLB_NoOfst]=MSLocObjPosVelCRLB(u,t,vu,vt,s,Qr,Qd,Qfr,Qfd)
% [CRLB_Jnt,CRLB_Ind,CRLB_NoOfst]=MSLocObjPosVelCRLB(u,t,vu,vt,s,Qr,Qd,Qfr,Qfd)
%
% This function obtains the CRLBs for the object position and velocity for 
% the cases of (i) using both indirect- and direct-path measurements,
% (ii) using indirect-path measurements only, and (iii) without time and 
% frequency offsets.
%
% Input parameter list:
% u:         (Dim x 1), object position.
% t:         (Dim x 1), transmitter position.
% vu:        (Dim x 1), object velocity.
% vt:        (Dim x 1), transmitter velocity.
% s:         (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:        (M x M), covariance matrix of the indirect range measurements.
% Qd:        (M x M), covariance matrix of the direct range measurements.
% Qfr:       (M x M), covariance matrix of the indirect range-rate measurements.
% Qfd:       (M x M), covariance matrix of the direct range-rate measurements.
% 
% Output parameter list:
% CRLB_Jnt:  (2 Dim x 2 Dim), CRLB of object location for joint estimation 
%            with transmitter location and offsets.
% CRLB_Ind:  (2 Dim x 2 Dim), CRLB of object location using indirect-path
%            measurements only without estimation of transmitter location.
% CRLB_NoOfst: (2 Dim x 2 Dim), CRLB of object location for joint estimation 
%              with transmitter location without time and frequency
%              offsets.

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
Hq=[H,zeros(M,M-1);zeros(M,M-1),H];

rho_dt=[];
rho_rt=[];
rho_ru=[];
rho_fru=[];
rho_frt=[];
rho_fdt=[];

for i=1:M
    rho_dti=(t-s(:,i))/norm(t-s(:,i));
    rho_dt=[rho_dt;rho_dti'];
    rho_rti=(t-u)/norm(t-u);
    rho_rt=[rho_rt;rho_rti'];
    rho_rui=(u-s(:,i))/norm(u-s(:,i))+(u-t)/norm(u-t);
    rho_ru=[rho_ru;rho_rui'];
    rho_fdti=vt/norm(t-s(:,i))-((t-s(:,i))'*vt*(t-s(:,i)))/(norm(t-s(:,i))^3);
    rho_fdt=[rho_fdt;rho_fdti'];
    rho_frui=vu/norm(u-s(:,i))-((u-s(:,i))'*vu*(u-s(:,i)))/(norm(u-s(:,i))^3)+(vu-vt)/norm(u-t)-((u-t)'*(vu-vt)*(u-t))/(norm(u-t)^3);
    rho_fru=[rho_fru;rho_frui'];
    rho_frti=(vt-vu)/norm(u-t)-((t-u)'*(vt-vu)*(t-u))/(norm(u-t)^3);
    rho_frt=[rho_frt;rho_frti'];
end

Q_MI=blkdiag(Qr,Qfr);
Q_MD=blkdiag(Qd,Qfd);

rho_MItheta=[rho_ru,zeros(M,K);rho_fru,rho_ru];
rho_MIphi=[rho_rt,zeros(M,K),ones(M,1),zeros(M,1);rho_frt,rho_rt,zeros(M,1),ones(M,1)];
rho_MDphi=[rho_dt,zeros(M,K),ones(M,1),zeros(M,1);rho_fdt,rho_dt,zeros(M,1),ones(M,1)];


rho_MIphin=[rho_rt,zeros(M,K);rho_frt,rho_rt];
rho_MDphin=[rho_dt,zeros(M,K);rho_fdt,rho_dt];


K_MIMD=Q_MI+rho_MIphi*inv(rho_MDphi'*inv(Q_MD)*rho_MDphi)*rho_MIphi';
K_MIMDn=Q_MI+rho_MIphin*inv(rho_MDphin'*inv(Q_MD)*rho_MDphin)*rho_MIphin';
K_MI=Hq*inv(Hq'*Q_MI*Hq)*Hq';

FIM=rho_MItheta'*inv(K_MIMD)*rho_MItheta;
FIMn=rho_MItheta'*inv(K_MIMDn)*rho_MItheta;
FIMq=rho_MItheta'*K_MI*rho_MItheta;

CRLB_Jnt=inv(FIMq);
CRLB_NoOfst=inv(FIMn);
CRLB_Ind=inv(FIM);

