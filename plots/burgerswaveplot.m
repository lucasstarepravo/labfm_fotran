
Re=100;
x=0:0.001:1;
y=-tanh((x.-0.5).*Re./2);

n=99;
A=load(strcat('../data_out/uv/uv',num2str(n)));



plot(x,y,'k-',A(:,1),A(:,3),'r.')
axis([0.4 0.6 -1.5 1.5])
set(gca,'fontname','times');set(gca,'fontsize',14)

xlabel('x');ylabel('u');


