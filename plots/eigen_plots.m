%% Script to plot L2norm... and print convergence rates to screen!!
clear all

C2=load('../data_out/eigens_C_2_nonoise');
Q2=load('../data_out/eigens_Q_2_nonoise');
W2=load('../data_out/eigens_W_2_nonoise');
G2=load('../data_out/eigens_G_2_nonoise');

C4=load('../data_out/eigens_C_4_nonoise');
Q4=load('../data_out/eigens_Q_4_nonoise');
W4=load('../data_out/eigens_W_4_nonoise');
G4=load('../data_out/eigens_G_4_nonoise');

C6=load('../data_out/eigens_C_6_nonoise');
Q6=load('../data_out/eigens_Q_6_nonoise');
W6=load('../data_out/eigens_W_6_nonoise');
G6=load('../data_out/eigens_G_6_nonoise');

%% NOISE FREE...
figure(1)
subplot(1,4,1),plot(C6(:,1),C6(:,2),'k.','linewidth',2,C4(:,1),C4(:,2),'r.','linewidth',2,C2(:,1),C2(:,2),'b.','linewidth',2);axis([-2 2 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Conic')
subplot(1,4,2),plot(1e14.*Q6(:,1),Q6(:,2),'k.','linewidth',2,1e14.*Q4(:,1),Q4(:,2),'r.','linewidth',2,1e14.*Q2(:,1),Q2(:,2),'b.','linewidth',2);axis([-2 2 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Quadratic')
text(0.9,-35,'x10^{-14}','Fontsize',14,'Fontname','Times')
subplot(1,4,3),plot(1e14.*W6(:,1),W6(:,2),'k.','linewidth',2,1e14.*W4(:,1),W4(:,2),'r.','linewidth',2,1e14.*W2(:,1),W2(:,2),'b.','linewidth',2);axis([-2 2 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Wendland')
text(0.9,-35,'x10^{-14}','Fontsize',14,'Fontname','Times')
subplot(1,4,4),plot(1e9.*G6(:,1),G6(:,2),'k.','linewidth',2,1e9.*G4(:,1),G4(:,2),'r.','linewidth',2,1e9.*G2(:,1),G2(:,2),'b.','linewidth',2);axis([-2 2 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Gaussian')
text(1,-35,'x10^{-9}','Fontsize',14,'Fontname','Times')

print('-depsc','-S1000,300',"eigens_4_nonoise.eps")
print('-dpng','-S1000,300',"eigens_4_nonoise.png")

%%%%%
clear all

C2=load('../data_out/eigens_C_2_noisy');
Q2=load('../data_out/eigens_Q_2_noisy');
W2=load('../data_out/eigens_W_2_noisy');
G2=load('../data_out/eigens_G_2_noisy');

C4=load('../data_out/eigens_C_4_noisy');
Q4=load('../data_out/eigens_Q_4_noisy');
W4=load('../data_out/eigens_W_4_noisy');
G4=load('../data_out/eigens_G_4_noisy');

C6=load('../data_out/eigens_C_6_noisy');
Q6=load('../data_out/eigens_Q_6_noisy');
W6=load('../data_out/eigens_W_6_noisy');
G6=load('../data_out/eigens_G_6_noisy');


%% NOISY
figure(2)
subplot(1,4,1),plot(C6(:,1),C6(:,2),'k.','linewidth',2,C4(:,1),C4(:,2),'r.','linewidth',2,C2(:,1),C2(:,2),'b.','linewidth',2);axis([-40 40 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Conic')
subplot(1,4,2),plot(Q6(:,1),Q6(:,2),'k.','linewidth',2,Q4(:,1),Q4(:,2),'r.','linewidth',2,Q2(:,1),Q2(:,2),'b.','linewidth',2);axis([-40 40 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Quadratic')
subplot(1,4,3),plot(W6(:,1),W6(:,2),'k.','linewidth',2,W4(:,1),W4(:,2),'r.','linewidth',2,W2(:,1),W2(:,2),'b.','linewidth',2);axis([-40 40 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Wendland')
subplot(1,4,4),plot(G6(:,1),G6(:,2),'k.','linewidth',2,G4(:,1),G4(:,2),'r.','linewidth',2,G2(:,1),G2(:,2),'b.','linewidth',2);axis([-40 40 -40 40])
set(gca,'Fontsize',18)
set(gca,'Fontname','Times')
xlabel('Re \lambda');ylabel('Im \lambda');
title('Gaussian')

print('-depsc','-S1000,300',"eigens_4_noisy.eps")
print('-dpng','-S1000,300',"eigens_4_noisy.png")

