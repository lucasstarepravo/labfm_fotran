
clear all

xx=-2:0.02:2;
yy=xx;

[X,Y]=meshgrid(xx,yy);
n=length(xx);

for i=1:n
   for j=1:n
      x=xx(j);
      y=yy(i);
   
      %% Hermite in X
      x2=x*x;x3=x*x*x;x4=x2*x2;
      Hx1 = 2.0*x;
      Hx2 = 4.0*x*x - 2.0d0;
      Hx3 = 8.0*x*x*x - 12.0*x;
      Hx4 = 16.0*x*x*x*x - 48.0*x*x + 12.0;
      Hx5 = 32.0*x2*x3 - 160.0*x3 + 120.0*x;
      Hx6 = 64.0*x3*x3 - 480.0*x2*x2 + 720.0*x2 - 120.0;
      Hx7 = 128.0*x4*x3 - 1344.0*x3*x2 + 3360.0*x3 - 1680.0*x;
      Hx8 = 256.0*x4*x4 - 3584.0*x4*x2 + 13440.0*x4 - 13440.0*x2 + 1680.0;
      % Hermite in Y
      y2=y*y;y3=y*y*y;y4=y2*y2;
      Hy1 = 2.0*y;
      Hy2 = 4.0*y*y - 2.0d0;
      Hy3 = 8.0*y*y*y - 12.0*y;
      Hy4 = 16.0*y*y*y*y - 48.0*y*y + 12.0;
      Hy5 = 32.0*y2*y3 - 160.0*y3 + 120.0*y;
      Hy6 = 64.0*y3*y3 - 480.0*y2*y2 + 720.0*y2 - 120.0;
      Hy7 = 128.0*y4*y3 - 1344.0*y3*y2 + 3360.0*y3 - 1680.0*y;
      Hy8 = 256.0*y4*y4 - 3584.0*y4*y2 + 13440.0*y4 - 13440.0*y2 + 1680.0;
      
      r2=x2 + y2;
      ff1 = (9.0/pi)*exp(-r2*9.0);
      G1(i,j) = ff1*Hx1*1/sqrt(2);
      G2(i,j) = ff1*Hx1*Hy1*0.5;
      G3(i,j) = ff1*Hx2*Hy1*0.5/sqrt(2);
      G4(i,j) = ff1*Hx3*Hy1*0.25;
      G5(i,j) = ff1*Hx3*Hy2*0.25/sqrt(2);            
      G6(i,j) = ff1*Hx4*Hy2*0.125;
     
      %% Original style...
      r2=max(r2,1e-15);
      r=sqrt(r2);r3=r2*r;r5=r3*r2;r7=r5*r2;r9=r7*r2;r11=r9*r2;
      x5=x4*x;y6=y4*y2;x6=x4*x2;
      ff1 = (-6/pi)*r*exp(-r2*9.0);          

      g1(i,j) = ff1*x/r;
      g2(i,j) = -x*y*ff1/r3;
      g3(i,j) = (2.0*x2*y - y3)*ff1/r5;
      g4(i,j) = (9.0*x*y3-6.0*x3*y)*ff1/r7;
      g5(i,j) = (63.0*x3*y2 - 36.0*x*y4 - 6.0*x5)*ff1/r9;   
      g6(i,j) = (477.0*x2*y4 - 408.0*x4*y2 - 36.0*y6 + 24.0*x6)*ff1/r11;
      
      %% Symmetric, GMLS, DCPSE
      ff1 = (9.0/pi)*exp(-r2*9.0);      
      gg1(i,j) = ff1*x;
      gg2(i,j) = ff1*x*y;
      gg3(i,j) = ff1*0.5*x2*y;
      gg4(i,j) = ff1*(1/6)*x3*y;
      gg5(i,j) = ff1*(1/12)*x4*y;
      gg6(i,j) = ff1*(1/48)*x4*y2;
      

     
   endfor
endfor   


H1=log(abs(G1));
H2=log(abs(G2));
H3=log(abs(G3));
H4=log(abs(G4));
H5=log(abs(G5));
H6=log(abs(G6));

figure(1)
colormap('autumn')
subplot(1,6,1),contour(X,Y,H1,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x')
subplot(1,6,2),contour(X,Y,H2,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('xy')
subplot(1,6,3),contour(X,Y,H3,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{2}y')
subplot(1,6,4),contour(X,Y,H4,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y')
subplot(1,6,5),contour(X,Y,H5,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y^{2}')
subplot(1,6,6),contour(X,Y,H6,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{4}y^{2}')


print('-depsc','-S1600,300','abfs.eps')

h1=log(abs(g1));
h2=log(abs(g2));
h3=log(abs(g3));
h4=log(abs(g4));
h5=log(abs(g5));
h6=log(abs(g6));

figure(2)
colormap('autumn')
subplot(1,6,1),contour(X,Y,h1,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x')
subplot(1,6,2),contour(X,Y,h2,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('xy')
subplot(1,6,3),contour(X,Y,h3,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{2}y')
subplot(1,6,4),contour(X,Y,h4,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y')
subplot(1,6,5),contour(X,Y,h5,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y^{2}')
subplot(1,6,6),contour(X,Y,h6,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{4}y^{2}')


print('-depsc','-S1600,300','abfs_o.eps')


hh1=log(abs(gg1));
hh2=log(abs(gg2));
hh3=log(abs(gg3));
hh4=log(abs(gg4));
hh5=log(abs(gg5));
hh6=log(abs(gg6));

figure(3)
colormap('autumn')
subplot(1,6,1),contour(X,Y,hh1,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x')
subplot(1,6,2),contour(X,Y,hh2,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('xy')
subplot(1,6,3),contour(X,Y,hh3,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{2}y')
subplot(1,6,4),contour(X,Y,hh4,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y')
subplot(1,6,5),contour(X,Y,hh5,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{3}y^{2}')
subplot(1,6,6),contour(X,Y,hh6,20);axis('square')
set(gca,'fontsize',14);set(gca,'fontname','times');title('x^{4}y^{2}')


print('-depsc','-S1600,300','abfs_t.eps')


