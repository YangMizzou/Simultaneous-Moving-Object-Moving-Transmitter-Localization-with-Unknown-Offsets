%
% This program will reproduce Fig. 6 or Fig. 8 in Y. Zhang and K. C. Ho, 
% "Multistatic moving object localization by a moving transmitter of 
% unknown location and offset," IEEE Trans. Signal Process., vol. 68, 
% pp. 4438-4453, 2020.
%
% The code generates Fig. 6 as it is.  To obtain Fig. 8, set 
%   CorrelatedNoise=1;
% on line 26.
%
% The code may not reproduce the figures exactly due to random number 
% generator seed setting or matlab version.
%
% Yang Zhang and K. C. Ho   08-22-2020
%
%       Copyright (C) 2020
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clear all;
warning('off');

CorrelatedNoise=0;

% -- configuration --
s=[0 1000;1000 0;-1000 0;0 -1000]';
to=[3000;2000];
uo=[2000;5000];
vuo=[-13 8]';
vto=[3 13]';
offsetdo=500;
offsetfo=10;
psio=[uo;to;offsetdo];

% -- simulation parameters --
K=length(uo);
M=size(s,2);    % number of receivers
L=5000;         % number of ensemble runs
nsePwrAlldB=[-10:5:30];
k=0.1;
rho=0.1;

% -- generate true measurement values --
do=sqrt(sum((repmat(to,1,M)-s(:,1:M)).^2))'+offsetdo;
ro=sqrt(sum((repmat(uo,1,M)-s(:,1:M)).^2))'+norm(uo-to)+offsetdo;
do_of=do-offsetdo;
ro_of=ro-offsetdo;
m_2=1/(2*M)*sum(do_of.^2+ro_of.^2);

% -- loop over noise power --
for nsePwrIdx=1:length(nsePwrAlldB),
    randn('state',319);
    fprintf('10 log10(nsePwr): %d\n',nsePwrAlldB(nsePwrIdx));
    nsePwr=10^(nsePwrAlldB(nsePwrIdx)/10);
    sigmar_2=nsePwr*(ro_of.^2/m_2);
    sigmad_2=nsePwr*(do_of.^2/m_2);
    Qr=diag(sigmar_2)/2;
    Qd=diag(sigmad_2)/2;
    
    %% ----- correlated noise, for Fig. 8 ----- %%
    if (CorrelatedNoise)
        for i=1:M-1
            for j=i+1:M
                Qr(i,j)=sqrt(sigmar_2(i))*sqrt(sigmar_2(j))*(norm(uo-to)^2)/((norm(uo-to)+norm(uo-s(:,i)))*(norm(uo-to)+norm(uo-s(:,j))));
                Qd(i,j)=sqrt(sigmad_2(i))*sqrt(sigmad_2(j))*rho;
            end
        end
    end;
    
    Qr=(Qr+Qr');
    Qd=(Qd+Qd');
    Qfr=k*Qr;
    Qfd=k*Qd;
    Qrfr=blkdiag(Qr,Qfr);
    
    % -- CRLBs for joint estimation using both indirect- and direct-path
    % measurements, for using indirect-path measurements only and for the
    % case without time offset --
    [CRLB_Jnt,CRLB_Ind,CRLB_NoOfst]=MSLocObjPosCRLB(uo,to,s,Qr,Qd);   
    CRLBu_Ind(nsePwrIdx)=trace(CRLB_Ind);
    CRLBu_NoOfst(nsePwrIdx)=trace(CRLB_NoOfst);
    CRLBu_Jnt(nsePwrIdx)=trace(CRLB_Jnt);
    
    SimuMseuCFS=0;
    SimuMseuIndMLE=0;
    SimuMseuMLE=0;
    
    % -- loop over ensemble runs --
    for ii=1:L,
        noiser=sqrtm(Qr)*randn(M,1);
        noised=sqrtm(Qd)*randn(M,1);
        r=ro+noiser;
        d=do+noised;
        
        % -- solution from indirect-path measurements using TDOA approach --
        initGuess=uo;     % true value initialization
        [u_hat]=MSLocObjPosInd_MLE(r,s,Qr,initGuess);
        
        % -- solution from MLE with true initialization --
        psiMLE=MSLocJntObjTxPos_MLE(r,d,s,Qr,Qd,psio);
        
        % -- proposed closed-form solution --
        [psi]=MSLocJntObjTxPos(r,d,s,Qr,Qd);
        
        % -- collecting MSE for perforamnce evaluation --
        SimuMseuIndMLE=SimuMseuIndMLE+norm((u_hat-uo),2).^2;
        SimuMseuCFS=SimuMseuCFS+norm((psi(1:K)-uo)).^2;
        SimuMseuMLE=SimuMseuMLE+norm((psiMLE(1:K)-uo),2).^2;
    end;
    
    MseuCFS(nsePwrIdx)=SimuMseuCFS/L;
    MseuIndMLE(nsePwrIdx)=SimuMseuIndMLE/L;
    MseuMLE(nsePwrIdx)=SimuMseuMLE/L;
    
end;

% -- plot result --
figure;
plot(nsePwrAlldB,10*log10(CRLBu_Jnt),'--r');
hold on;
plot(nsePwrAlldB,10*log10(CRLBu_Ind),'-r');
hold  on;
plot(nsePwrAlldB,10*log10(CRLBu_NoOfst),'-k','LineWidth',1);
hold  on;
plot(nsePwrAlldB,10*log10(MseuCFS),'dk','MarkerSize',8);
hold on;
plot(nsePwrAlldB,10*log10(MseuMLE),'*m','MarkerSize',8);
hold on;
plot(nsePwrAlldB,10*log10(MseuIndMLE),'+b','MarkerSize',8);
hold off;
grid on; xlabel('10 log10(\sigma^2(m^2))'); ylabel('10 log10(MSE(m^2))');
legend('CRLB by TDOA','CRLB for joint estimation (offsets)','CRLB for joint estimation (no offsets)','Proposed solution','IMLE for joint estimaton','IMLE by TDOA','location','northwest');
ylim([0 100]);

