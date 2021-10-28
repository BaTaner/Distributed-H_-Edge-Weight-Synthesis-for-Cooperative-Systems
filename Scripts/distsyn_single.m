function [resstp1] = distsyn_single(Hs,na, alpha_c,gamma_bound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   Hs = ss(a,b,c,d) represents the state space definition of single agent.
%   This system can be multi input multi output
%   na      defines the number of agents in the cooperative system
%   Adj_0   defines the initial value for the communication weights.
%   For this work adjacency matrix is fixed.

[a,b,c,d] = ssdata(Hs);
nws = size(b,2);
nzs = nws;
nw  = nws*na;
nd  = nzs*na;
nx  = size(a,1); 
A = zeros(size(a));
B = zeros(size(b,1), size(b,2)*(nw));
C = zeros(size(c,1)*(nw), size(c,2));
D = zeros(size(d)*(nw));
ne  =   na+1; % number of agents including Hc
sz.nw   =   nw;
sz.na   =   na;
sz.ne   =   ne;
sz.nx   =   nx;

for i =1:nw
    for j=1:nw
        if (j~=i && j~=nw)
            Adj(i,j) = (j)*alpha_c;
        elseif (j~=i && j==nw)
            Adj(i,j) = 1-(sum(Adj(i,1:nw-1)));
        end
    end
end
Adj(nw,nw-1) = 1-sum(Adj(nw,1:nw-2));

for i = 1: ne
    if i<ne
        H(i).A      = a;
        H(i).B      = B;
        H(i).B(:,i) = b;
        H(i).C      = C;
        H(i).C(i,:) = c;
        H(i).D      = D;
        H(i).D(i,i) = d;
    else
        H(i).A              = A;
        H(i).B              = B;
        H(i).C              = C;
        H(i).D              = D;
        H(i).D(1:na,1:na)   = Adj;
    end
end

% We need to extend every system H_{i} with zeros for dimensional compliance
% following matrices will give the dimensions to the extended matrices
Ae  =   zeros(nx);
Be  =   zeros(nx,nw*ne);
Ce  =   zeros(nw*ne,nx);
De  =   zeros(nw*ne);

for i = 1:ne
    Hez(i).A     = Ae;
    Hez(i).B     = Be;
    Hez(i).C     = Ce;
    Hez(i).D     = De;
    Hez(i).I     = De;
    Hez(i).A     = H(i).A;
    Hez(i).B(:,(i-1)*nw+1:i*nw)  = H(i).B;
    Hez(i).C((i-1)*nw+1:i*nw,:)  = H(i).C;
    Hez(i).D((i-1)*nw+1:i*nw,(i-1)*nw+1:i*nw)  = H(i).D;
    Hez(i).I((i-1)*nw+1:i*nw,(i-1)*nw+1:i*nw)  = eye(nw);
end


gam_lw      =   gamma_bound(1);
gam_up      =   gamma_bound(2);
gam_tol     =   0.01;
gam_err     =   inf;
gam_old     =   0;
max_iter	=   20;
gam_cnt     =   1;
time_overall    = cputime;
while true
    gam_try     =	(gam_lw + gam_up)*0.5;
    t           =   -0.00001;
    [resstp1_n]	=   solverstp_single(Hs,Hez,sz,gam_try,t);
    if resstp1_n.succeed
        gam_up	=	gam_try;
    else
        gam_lw	=	gam_try;
    end
    gam_err   =   abs(gam_old - gam_try);
    if resstp1_n.succeed && gam_err<gam_tol
        resstp1_n.gam_err      =   gam_err;
        resstp1_n.gam_iter     =   gam_cnt;
        resstp1_n.Gk           =   gam_try;
        resstp1                =   resstp1_n;
        resstp1.alpha          =   alpha_c;
        break
    elseif gam_cnt>max_iter
        break
    end
    gam_cnt           =   gam_cnt +  1;
    gam_old           =   gam_try;  
end
resstp1.time    =   cputime - time_overall;
end

