clear
clc
y1=load('8_BA50_4t005nL1B1degreecascadingfailNorefine.mat','cascLfailpercentcycleT');
    y1=(y1.cascLfailpercentcycleT);


y2=load('8_BA50_4t005nL1B1degreecascadingReTstrong.mat','cascLfailpercentcycleT');
    y2=(y2.cascLfailpercentcycleT);
y3=load('8_BA50_4t005nL1B1degreecascadingfailThalf1weak.mat','cascLfailpercentcycleT');
    y3=(y3.cascLfailpercentcycleT);

% y4=load('4_BA100_0t005nL02B1degreecascadingReTstrong_weak.mat','cascLfailpercentcycleT');
y4=load('8_BA50_4t005nL1B1degreecascadingfailThalf2strong.mat','cascLfailpercentcycleT');
    y4=(y4.cascLfailpercentcycleT);
y5=load('random8_BA50_4t005nL1B1degreecascadingReTstrong.mat','cascLfailpercentcycleT')
    y5=(y5.cascLfailpercentcycleT);
x=y1(:,1);
x0=x(2,1)-x(1,1);
y1=y1(:,2);
y2=y2(:,2);
y3=y3(:,2);
y4=y4(:,2);
y5=y5(:,2);
% y0=find(y1+y2+y3+y4==0);
y0(2)=12;
figure(3)
set(0,'defaultfigurecolor','w');
plot(x(1:y0(2)),y1(1:y0(2)),'k-*','MarkerFaceColor','k','markersize',5,'linewidth',1.5);hold on;
plot(x(1:y0(2)),y2(1:y0(2)),'m-h','MarkerFaceColor','m','markersize',5,'linewidth',1.5);hold on;
plot(x(1:y0(2)),y3(1:y0(2)),'g-s','MarkerFaceColor','g','markersize',5,'linewidth',1.5);hold on;
plot(x(1:y0(2)),y4(1:y0(2)),'b-d','MarkerFaceColor','b','markersize',5,'linewidth',1.5);hold on;
plot(x(1:y0(2)),y5(1:y0(2)),'r-p','MarkerFaceColor','r','markersize',5,'linewidth',1.5);hold on;
R=size(y1,1);
str=[num2str(x,3) num2str(y1,2)];% str=[repmat('X:',R,1) num2str(x) repmat(',Y:',R,1) num2str(y1)];
% text(x,y1,cellstr(str));
% text(x,y2,'o','color','r','fontsize',10)
xlabel('Key threshold value $ \hat{\alpha} $','Interpreter','latex','FontSize',12)
ylabel({'The complete breakdown ratio ${P_r}$ '},'Interpreter','latex','FontSize',12)
set(gca,'XTick',min(x):x0:x(y0(2)))
set(gca,'YTick',0:0.05:1)
h=legend(['Attack totally without defense strategy'],['Bi-directional rewiring strategy'],['One way rewiring  strategy 1'],['One way rewiring strategy 2'],['Bi-directional random rewiring strategy']);
set(h,'Box','off');
set(h,'Fontsize',12);
% title('8_BA50_4 nl=1 degree ∂‘±»Õº','FontSize',12);
% grid on;
% figure(2)
% plot(x,y2,'g-o','MarkerFaceColor','r','markersize',4)
% set(gca,'XTick',min(x):0.02:max(x))
% set(gca,'YTick',min(y2):0.01:max(y2))


% gtext('sin(x)')
