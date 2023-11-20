%% create a run of plots so we can make an animation, just line graphs...

clear all
nn=100;

% Loop over nn output files
for n=1:nn

% Load data
fn=strcat('../data_out/uv/uv',num2str(n));
A=load(fn);
%nx = A(1,1);dx = A(1,2);xmin=A(1,3);
%A(1,:) = [];
x = A(:,1);u=A(:,3);v=A(:,4);

% Plot the data                %%% sawtooth - go -1 1
figure(1)
%%{
%hold off
%sawtooth_analytical(200,0.01*(n-1));hold on
plot(x,u,'k.','linewidth',1.5);axis([0 1 -1 1])
set(gca,'Fontsize',16);set(gca,'Fontname','Times');
axis([0 1 -1 1])
%[hleg1, hobj1] = legend("Re=10","Re=20","Re=50","Re=100","Re=500");
%set(hleg1,'position',[0.65 0.7 0.25 0.2]);
xlabel('x/\lambda');ylabel('u');
%}

pause(0.01);


% print the plots
%print('-dpng','%d.png',num2str(n))


endfor


%figure(2)
%surf(x,y,v)


