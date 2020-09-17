function [solStg2,solStg1]=MSLocJntObjTxPosVel(r,d,fr,fd,s,Qr,Qd,Qfr,Qfd)
% [solStg2,solStg1]=MSLocJntObjTxPosVel(r,d,fr,fd,s,Qr,Qd,Qfr,Qfd)
%
% This function realizes the algebraic closed-form solution for jointly
% estimating the unknown object and transmitter positions and velocities
% using both indirect- and direct-path range and range-rate measurements.
%
% Input parameter list:
% r:      (M x 1), indirect-path range measurements.
% d:      (M x 1), direct-path range measurements.
% fr:     (M x 1), indirect-path range-rate measurements.
% fd:     (M x 1), direct-path range-rate measurements.
% s:      (Dim x M), receiver position matrix, M is the number of receivers.
% Qr:     (M x M), covariance matrix of the indirect range measurements.
% Qd:     (M x M), covariance matrix of the direct range measurements.
% Qfr:    (M x M), covariance matrix of the indirect range-rate measurements.
% Qfd:    (M x M), covariance matrix of the direct range-rate measurements.
%
% Output parameter list:
% solStg2: (4*Dim+2 x 1), final solution,
%          [objPos;ObjVel;txPos;txVel;ofStTm;ofStFq]^T.
% solStg1: (4*Dim+2+6 x 1), stage1 solution.

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

for i=1:M
    hr(i,:)=0.5*(r(i)^2-s(:,i)'*s(:,i));
    hfr(i,:)=r(i)*fr(i);
    hd(i,:)=0.5*(d(i)^2-s(:,i)'*s(:,i));
    hfd(i,:)=d(i)*fd(i);
    Gr(i,:)=[-s(:,i)' zeros(1,K) zeros(1,K) zeros(1,K) r(i) 0 r(i) 1 0.5 0 0 0];
    Gfr(i,:)=[zeros(1,K) -s(:,i)' zeros(1,K) zeros(1,K) fr(i) r(i) fr(i) 0 0 r(i) 1 1];
    Gd(i,:)=[zeros(1,K) zeros(1,K) -s(:,i)' zeros(1,K) d(i) 0 0 0 0.5 0 0 0];
    Gfd(i,:)=[zeros(1,K) zeros(1,K) zeros(1,K) -s(:,i)' fd(i) d(i) 0 0 0 0 0 1];
end

h1=[hr;hfr;hd;hfd];
G1=[Gr;Gfr;Gd;Gfd];
Q=blkdiag(Qr,Qfr,Qd,Qfd);

B1=eye(4*M);
W1=B1*Q*B1';

phi1=inv(G1'*inv(W1)*G1)*G1'*inv(W1)*h1;
% phi1=(G1'*inv(W1)*G1)\(G1'*inv(W1)*h1);

for k=1:3,  % repeating 3 times to improve the weighting matrix
    for i=1:M
        B1(i,i)=norm(phi1(1:K)-s(:,i));
        B1(i+M,i)=((phi1(1:K)-s(:,i))/norm(phi1(1:K)-s(:,i)))'*phi1(K+1:2*K);
        B1(i+M,i+M)=norm(phi1(1:K)-s(:,i));
        B1(i+2*M,i+2*M)=norm(phi1(2*K+1:3*K)-s(:,i));
        B1(i+3*M,i+2*M)=((phi1(2*K+1:3*K)-s(:,i))/norm(phi1(2*K+1:3*K)-s(:,i)))'*phi1(3*K+1:4*K);
        B1(i+3*M,i+3*M)=norm(phi1(2*K+1:3*K)-s(:,i));
    end
    W1=B1*Q*B1';
    phi1=inv(G1'*inv(W1)*G1)*G1'*inv(W1)*h1;
    % phi1=(G1'*inv(W1)*G1)\(G1'*inv(W1)*h1);
    
end
solStg1=phi1(1:4*K+2);

B21=eye(4*K+2);
B22=zeros(4*K+2,6);
B23=[(phi1(2*K+1:3*K)-phi1(1:K))',zeros(1,K),-(phi1(2*K+1:3*K)-phi1(1:K))',zeros(1,K),0,0;
    -phi1(2*K+1:3*K)',zeros(1,K),(2*phi1(2*K+1:3*K)-phi1(1:K))',zeros(1,K),0,0;
    zeros(1,K),zeros(1,K),-phi1(2*K+1:3*K)',zeros(1,K),phi1(4*K+1),0;
    (phi1(3*K+1:4*K)-phi1(K+1:2*K))',(phi1(2*K+1:3*K)-phi1(1:K))',-(phi1(3*K+1:4*K)-phi1(K+1:2*K))',-(phi1(2*K+1:3*K)-phi1(1:K))',0,0;
    -phi1(3*K+1:4*K)',-phi1(2*K+1:3*K)',(2*phi1(3*K+1:4*K)-phi1(K+1:2*K))',(2*phi1(2*K+1:3*K)-phi1(1:K))',0,0;
    zeros(1,K),zeros(1,K),-phi1(3*K+1:4*K)',-phi1(2*K+1:3*K)',phi1(4*K+2),phi1(4*K+1)];
B24=[2*phi1(4*K+3),0,0,0,0,0;
    2*phi1(4*K+1),2,0,0,0,0;
    0,0,1,0,0,0;
    2*phi1(4*K+6),0,0,2*phi1(4*K+3),0,0;
    2*phi1(4*K+2),0,0,2*phi1(4*K+1),2,0;
    0,0,0,0,0,2];
B2=[B21,B22;B23,B24];

h2=[phi1(1:4*K+2);phi1(4*K+3)^2;2*phi1(4*K+4);phi1(4*K+5);2*phi1(4*K+3)*phi1(4*K+6);2*phi1(4*K+7);2*phi1(4*K+8)];

G2=[eye(4*K+2);
    (phi1(1:K)-phi1(2*K+1:3*K))',zeros(1,K),-(phi1(1:K)-phi1(2*K+1:3*K))',zeros(1,K),0,0;
    phi1(2*K+1:3*K)',zeros(1,K),-(2*phi1(2*K+1:3*K)-phi1(1:K))',zeros(1,K),-2*phi1(4*K+3),0;
    zeros(1,K),zeros(1,K),phi1(2*K+1:3*K)',zeros(1,K),-phi1(4*K+1),0;
    (phi1(K+1:2*K)-phi1(3*K+1:4*K))',(phi1(1:K)-phi1(2*K+1:3*K))',-(phi1(K+1:2*K)-phi1(3*K+1:4*K))',-(phi1(1:K)-phi1(2*K+1:3*K))',0,0;
    phi1(3*K+1:4*K)',phi1(2*K+1:3*K)',(phi1(K+1:2*K)-2*phi1(3*K+1:4*K))',(phi1(1:K)-2*phi1(2*K+1:3*K))',-2*phi1(4*K+6),-2*phi1(4*K+3);
    zeros(1,K),zeros(1,K),phi1(3*K+1:4*K)',phi1(2*K+1:3*K)',-phi1(4*K+2),-phi1(4*K+1)];

W2=B2*inv(G1'*inv(W1)*G1)*B2';

phi2=inv(G2'*inv(W2)*G2)*G2'*inv(W2)*h2;

solStg2=phi2;

