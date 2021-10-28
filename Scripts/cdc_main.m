
yalmip('clear')
close all
clc

% Inputs to the algorithm
Hs      = ss([0,1;-2,-2],[0;1],[1,0],0);
na      = 4;
alpha_c = 0.05;
gamma_bound     =   [0,100];
% Benchmark
gamma_range   =   0.01:0.01:0.21;
for i = 1:length(gamma_range)
    alpha_i     = gamma_range(i);
    [resstp,sz] = distsyn1(Hs,na, alpha_i,gamma_bound);
    bench(i).alpha      = alpha_i;
    bench(i).gamma      = resstp.Gk;
    bench(i).TimeIter   = resstp.time;
    bench(i).error      =  0;
end
% Check relative interior
[resstp0,sz0] = distsyn1(Hs,na, alpha_c,gamma_bound);
% Start main algorithm
store(1).alpha      = alpha_c;
store(1).gamma      = resstp0.Gk;
store(1).TimeIter   = resstp0.time;
store(1).error      =  0;
maxiter = 5;
tolIter = 0.1;
i = 2;
time_dist = cputime;
while i <= maxiter
    [resstp] = distsyn_single(Hs,na, store(i-1).alpha,gamma_bound);
     store(i).alpha     = resstp.alpha_syn;
     store(i).gamma     = resstp.Gk;
     store(i).TimeIter  = resstp.time;
     store(i).error     =   store(i).alpha - store(i-1).alpha;
     if abs(store(i).error) < tolIter
         TotalTime	= 0;
         for j = 1:size(store,2)
             TotalTime	= TotalTime + store(j).TimeIter;
         end
         break
     end
     i = i+1;
end
time_dist = cputime - time_dist;
% Main algorithm finishes here
%% Plotting

nw      = sz0.nw;
nStore  = size(store,2);
nBench  = size(bench,2);
alpha_0     =   store(1).alpha;
alpha_syn   =   store(nStore).alpha;
for i =1:nw
    for j=1:nw
        if (j~=i && j~=nw)
            Adj_syn(i,j) = (j)*alpha_syn;
            Adj_0(i,j) = (j)*alpha_0;
        elseif (j~=i && j==nw)
            Adj_syn(i,j) = 1-(sum(Adj_syn(i,1:nw-1)));
            Adj_0(i,j) = 1-(sum(Adj_0(i,1:nw-1)));
        end
    end
end
Adj_syn(nw,nw-1) = 1-sum(Adj_syn(nw,1:nw-2));
Adj_0(nw,nw-1) = 1-sum(Adj_0(nw,1:nw-2));

for i=1:nStore
    alp_iter(i,1) = store(i).alpha;
    gam_iter(i,1) = store(i).gamma;
    err_iter(i,1) = abs(store(i).error);
end
iter = 1:nStore;

for i=1:nBench
   alp_range(i,1) = bench(i).alpha;
   gam_range(i,1) = bench(i).gamma;
end

for i=2:nStore-1
    alp_iter(i,1) = store(i).alpha;
    gam_iter(i,1) = store(i).gamma;
end

figure
subplot(3,1,1)
plot(store(1).alpha,store(1).gamma,'Marker','o','Color','red','LineWidth',2,'LineStyle','none','MarkerSize',10);
hold on
plot(alp_range,gam_range,'Marker','d','Color','black','LineWidth',1.3);
plot(store(nStore).alpha,store(nStore).gamma,'Marker','x','Color','blue','LineWidth',2,'LineStyle','none','MarkerSize',20);
% for i=2:nStore-1
%     plot(store(i).alpha,store(i).gamma,'Marker','s','Color','cyan','LineWidth',1.5,'LineStyle','none');
% end

hold off
grid
xlabel('\alpha');
ylabel('\gamma');
title('Synthesis Results')
% lgnd = legend('\gamma_{0}','\gamma_{range}','\gamma_{syn}','\gamma_{iter}');
lgnd = legend('\gamma_{0}','\gamma_{range}','\gamma_{syn}');
set(lgnd,'location','southwest')
ylim([0,store(1).gamma+0.2]);
xlim([0,0.2]);
subplot(3,1,2)
plot(iter',gam_iter,'Marker','x','Color','black','LineWidth',1);
xlim([0.8 nStore+0.2]);
ylim([0 2]);
xlabel('Iteration number')
ylabel('\gamma')
title('Trace of \gamma over the iterations')
grid on
subplot(3,1,3)
plot(iter',alp_iter,'Marker','x','Color','black','LineWidth',1);
xlim([0.8 nStore+0.2]);
ylim([0 0.25]);
xlabel('Iteration number')
ylabel('\alpha')
title('Trace of \alpha over the iterations')
grid on
figure
plot(store(1).alpha,store(1).gamma,'Marker','o','Color','red','LineWidth',2,'LineStyle','none','MarkerSize',10);
hold on
plot(alp_range,gam_range,'Marker','d','Color','black','LineWidth',1.3);
plot(store(nStore).alpha,store(nStore).gamma,'Marker','x','Color','blue','LineWidth',2,'LineStyle','none','MarkerSize',20);
for i=2:nStore-1
    plot(store(i).alpha,store(i).gamma,'Marker','s','Color','cyan','LineWidth',1.5,'LineStyle','none');
end
%% another Figure

figure
plot(alp_range,gam_range,'color','blue','Marker','x','LineWidth',2)
hold on
plot(store(1).alpha,store(1).gamma,'color','Cyan','Marker','d','MarkerSize',14,'LineStyle','none','MarkerFaceColor','Magenta');
plot(alp_iter,gam_iter,'color','black','Marker','s','MarkerSize',6,'LineStyle','--','MarkerFaceColor','Green');

plot(store(nStore).alpha,store(nStore).gamma,'color','black','Marker','d','MarkerSize',10,'LineStyle','--','MarkerFaceColor','Red');
hold off
% title('Optimal Adjacency for given cooperative system with H_{\infty} Performance Criteria')
Leg = legend('H_{\infty} Performance along \alpha',...
        '\gamma^{0}','\gamma^{t}','\gamma^{final}');
set(Leg, 'Location', 'Best')
title('Synthesized \alpha vs \gamma')
xlabel('\alpha')
ylabel('\gamma')
% xlim([0 alfa(end)]);
% ylim([min(gamma0_bt)-0.2 max(gamma0_bt)+0.2]);
grid on