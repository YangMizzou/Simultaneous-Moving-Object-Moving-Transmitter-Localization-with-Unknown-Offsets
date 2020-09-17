function [solStg2,solStg1]=MSLocJntObjTxPos(r,d,s,Qr,Qd)
% [solStg2,solStg1]=MSLocJntObjTxPos(r,d,s,Qr,Qd)
%
% This function realizes the algebraic closed-form solution for jointly
% estimating the unknown object and transmitter positions using both 
% indirect- and direct-path range measurements.
%
% Input parameter list:
% r:      (M x 1), indirect-path range measurements.
% d:      (M x 1), direct-path range measurements.
% s:      (Dim x M), receiver position matrix, M is the number of receivers.         
% Qr:     (M x M), covariance matrix of the indirect range measurements.
% Qd:     (M x M), covariance matrix of the direct range measurements.
% 
% Output parameter list:
% solStg2: (2*Dim+1 x 1), final solution, 
%          [objPos;txPos;ofStTm]^T.
% solStg1: (2*Dim+2+3 x 1), stage1 solution.

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

for i=1:M
    hr(i,:)=0.5*(r(i)^2-s(:,i)'*s(:,i));
    hd(i,:)=0.5*(d(i)^2-s(:,i)'*s(:,i));
    Gr(i,:)=[-s(:,i)' zeros(1,K) r(i) r(i) 1 0.5]; 
    Gd(i,:)=[zeros(1,K) -s(:,i)' d(i) 0 0 0.5];
end
h1=[hr;hd];
G1=[Gr;Gd];
Q=blkdiag(Qr,Qd);
B1=eye(2*M);
W1=B1*Q*B1';
phi1=inv(G1'*inv(W1)*G1)*G1'*inv(W1)*h1;
% phi1=(G1'*inv(W1)*G1)\(G1'*inv(W1)*h1);

for k=1:3,   % repeating 3 times to improve the weighting matrix
    for i=1:M
        B1(i,i)=norm(phi1(1:K)-s(:,i));
        B1(i+M,i+M)=norm(phi1(K+1:2*K)-s(:,i));
    end
    W1=B1*Q*B1';
    phi1=inv(G1'*inv(W1)*G1)*G1'*inv(W1)*h1;
    % phi1=(G1'*inv(W1)*G1)\(G1'*inv(W1)*h1);
end

solStg1=phi1;

B21=eye(2*K+1);
B22=zeros(2*K+1,3);
B23=[(phi1(K+1:2*K)-phi1(1:K))',-(phi1(K+1:2*K)-phi1(1:K))',0;
    -phi1(K+1:2*K)',(2*phi1(K+1:2*K)-phi1(1:K))',0;
    zeros(1,K),-phi1(K+1:2*K)',phi1(2*K+1)];
B24=[2*phi1(2*K+2),0,0;
    2*phi1(2*K+1),2,0;
    0,0,1];
B2=[B21,B22;B23,B24];

h2=[phi1(1:2*K+1);phi1(2*K+2)^2;2*phi1(2*K+3);phi1(2*K+4)];

G2=[eye(2*K+1);
    (phi1(1:K)-phi1(K+1:2*K))',-(phi1(1:K)-phi1(K+1:2*K))',0;
    phi1(K+1:2*K)',(phi1(1:K)-2*phi1(K+1:2*K))',-2*phi1(2*K+2);
     zeros(1,K),phi1(K+1:2*K)',-phi1(2*K+1)];

W2=B2*inv(G1'*inv(W1)*G1)*B2';

phi2=inv(G2'*inv(W2)*G2)*G2'*inv(W2)*h2;

solStg2=phi2;