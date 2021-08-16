a=str2mat('Edge_betweenness_pic8_SW500_8_05t008nL1B1degreepredictionedgeEsasort.mat',...
'pic8_BA50_4t005nL1B1degreepredictionedge.mat',...
'pic20_BA100_11t005nL1B1degreepredictionedge.mat',...
'pic8_BA100_4t005nL1B1degreepredictionedge.mat',...
'pic6_BA100_3t005nL1B1degreepredictionedge.mat',...
'pic4_BA100_2t01nL1B1degreepredictionedge.mat',...
'Edge_betweenness_pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat',...
'pic4_BA100_2t005nL02B1degreepredictionedge.mat',...
'pic4_BA500_2t005nL1B1degreepredictionedge.mat',...
'pic6_BA500_3t005nL1B1degreepredictionedge.mat',...
'pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat',...
'pic8_SW500_8_05t008nL1B1degreepredictionedgeEsasort.mat',...
'pic8_BA500_4t01nL1B1degreepredictionedge.mat',...
'pic8_BA500_4t005nL1B1degreepredictionedge.mat',...
'pic4_SW100_4_05t005nL1B1degreepredictionedge.mat',...
'pic4_SW100_4_08t005nL1B1degreepredictionedge.mat',...
'pic4_SW100_4_08t005nL02B1degreepredictionedge.mat',...
'pic4_SW100_4_08t02nL02B1degreepredictionedge.mat',...
'pic6_SW500_6_08t005nL1B1degreepredictionedge.mat',...
'pic4_SW100_4_08t01nL1B1degreepredictionedge.mat',...
'pic4_SW500_4_08t01nL1B1degreepredictionedge.mat',...
'pic6_SW500_6_08t01nL1B1degreepredictionedge.mat',...
'pic8_SW500_8_05t01nL1B1degreepredictionedge.mat',...
'pic4_BA500_2t01nL1B1degreepredictionedge.mat',...
'pic4_BA500_2t01nL1B1betweenpredictionedge.mat',...
'pic8_BA500_4t009nL1B1degreepredictionedgeEsasort.mat')
 Ynumber=[];
 exasorder=1;
 if  exasorder==1
    figure(14)
else
    figure(15)
end
for ii=13%11是BA 12是SW
cascLfailtheta1=load(a(ii,:),'cascLfailtheta');
cascLfailtheta=cascLfailtheta1.cascLfailtheta;
% plot3(cascLfailtheta(:,4),cascLfailtheta(:,5),cascLfailtheta(:,3),'o')
% figure (2)10
% surf(cascLfailtheta(:,[4 5 3]))
% mesh(cascLfailtheta(:,[5 4 3]));
nn=5;
ynumber=[];
if  exasorder==1
    j=5;
else
    j=4;
end
cx=cascLfailtheta(:,5);
cy=cascLfailtheta(:,4);
cz=cascLfailtheta(:,3);
% plot(cy,cz);
sApProMax=max(cascLfailtheta(:,j));
sApProMin=min(cascLfailtheta(:,j));
if ~any(sApProMax)
    continue;
end
% x=sort(unique(cascLfailtheta(:,j)),1);
x=linspace(sApProMin,sApProMax+0.01,nn);
x=[0,x]
% x=linspace(0,3,nn+1)
cascLfailnumberunique=sort(unique(cascLfailtheta(:,3)),1);%[ 10    20    30    40    50    60    70    80    90   100]
% cascLfailnumberunique=[ 2  4   6  7  8  9  10  11  12   15 ];

for i=1:1:nn
[fx,fy]=find(cascLfailtheta(:,j)>=x(i)&cascLfailtheta(:,j)<x(i+1));
% [fx,fy]=find(cascLfailtheta(:,2)==x(i));
cascLfailnumber=cascLfailtheta(fx,3);
yy=hist(cascLfailtheta(fx,3),cascLfailnumberunique);
if ~any(cascLfailnumber)
yy=yy';
end
ynumber=[ynumber;yy];
end
Ynumber=[Ynumber; [ynumber*cascLfailnumberunique./(ynumber*ones(size(ynumber,2),1))]'];
           
% figure(ii)

        h=bar3(ynumber,'grouped ');
        view(-90,0)
%  ccolormat=[0	0	0;1	1	1;0.5 0.5 0.5;0	1	0;0	0 1;1	1	0;1	0	1;0	1	1;0.667	0.667 1;1	0	0;];
% colorh=[ccolormat ;ccolormat ;ccolormat ;ccolormat ]';
%  for n=1:numel(h)
% set(h(n),'facecolor',colorh(:,n))
% end
for i=1:length(x)-1
    xx(i)=1/2*(x(i)+x(i+1))
end
set(gca, 'YTickLabel', cascLfailnumberunique)
 set(gca,'YTickLabel', xx)
 legend(num2str(cascLfailnumberunique))% legend('1',sprintf('\n'),'suibian')% numel(h)
zlabel('Numbers of global breakdown')
% title(strcat(a(ii,(4:20)),'The correlation between variable factor and Cascading failure of Periodic numbers'))
% ylabel({'Triggered the first stage neighbouring edge' ,'maximized α'})
%  ylabel('Edge weight ratio w_{ij}k_{j}/s_{j}')

%         % Periodic numbers
        if  exasorder==1
        ylabel('Edge weight ratio w_{ij}k_{j}/s_{j}')
        else
        ylabel({'Triggered the first stage neighbouring edge' ,'maximized α'})
        end
end
% ccolormat=[0	0	0;1	1	1;0.5 0.5 0.5;0	1	0;0	0 1;1	1	0;1	0	1;0	1	1;0.667	0.667 1;1	0	0;];
% hot1=hot;
% % hot1=[ccolormat;ccolormat;ccolormat];
% figure(size(a,1)+1)
% Ynumber1=Ynumber;
%  Ynumber1=Ynumber-[Ynumber(:,1),Ynumber(:,1),Ynumber(:,1)]
% for ij=1:1:21
% plot(Ynumber1(ij,:),'ro-','LineWidth',1,'MarkerSize',1.5,'MarkerFaceColor','r')
% hold on
% end
% legend(num2str([1:size(Ynumber1,1)]'))
% ss=tabulate(cascLfailtheta(:,1));
% ss=sortrows(ss,-3);
% ss=ss(:,1);
% bb=textread('F:\py\Mtb\6_BA500_3\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
% bb=[[1:size(bb,1)]',bb];
% bb=sortrows(bb,-2);
% bc=textread('F:\py\Mtb\6_BA500_3\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
% bc=[[1:size(bc,1)]',bc];
% bc=sortrows(bc,-2);
% bd=textread('F:\py\Mtb\6_BA500_3\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
% bd=[[1:size(bd,1)]',bd];
% bd=sortrows(bd,-2);
% [resumbb,rerowsbb]=ismember(ss,bb(:,1),'rows');
% [resumbc,rerowsbc]=ismember(ss,bc(:,1),'rows');
% [resumbd,rerowsbd]=ismember(ss,bd(:,1),'rows');
% rss=[1:size(ss,1)]';
% bbcd=[rss,rerowsbb,rerowsbc,rerowsbd];
% plot(rss(1:10,1),rss(1:10,1),'m',rss(1:10,1),rerowsbb(1:10,1),' r',rss(1:10,1),rerowsbc(1:10,1),' g',rss(1:10,1),rerowsbd(1:10,1),' b')
% legend('正比例函数','betweenness_centrality','closeness_centrality','degree_centrality');
% subplot(2,1,1)
% figure(22)
% [X,Y]=meshgrid(1:3,1:21);
% meshz(X,Y,Ynumber1);
% 
% % plot3(X',Y',Ynumber1','.r--')
% set(gca,'xtick',[0:20:360])
% axis ([1 3 1 22 -2 5]); 
% colormap(mymap)
% subplot(2,1,2)


% figure(23)
% % % 画循环周期的图
% [X,Y]=meshgrid(1:3,1:21);
% mesh(X',1*ones(size(Y')),Ynumber1')
% axis ([1 3 1 2 -2 5]); 
% view(0,0)
% view(-90,0)





% h=bar3(ynumber);
% co=[0 1 1;0 1 0]
% for n=1:numel(h)
% set(h(n),'facecolor',co(n,:))
% end
