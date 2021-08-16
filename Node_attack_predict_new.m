% a=str2mat('pic8_BA50_4t005nL1B1degreepredictionedge.mat',...%一般
% 'pic20_BA100_11t005nL1B1degreepredictionedge.mat',...%效果可以 很好 bu
% 'pic8_BA100_4t005nL1B1betweenpredictionedge.mat',...%效果可以 很好 不
% 'pic6_BA100_3t005nL1B1degreepredictionedge.mat',...%效果可以  很好 bu
%'pic6_BA100_3t005nL02B1betweenpredictionedge.mat'   %很好
% 'pic4_BA100_2t005nL1B1degreepredictionedge.mat',...%效果可以
% 'pic4_BA100_2t005nL02B1degreepredictionedge.mat',...%很好
% 'pic4_BA500_2t005nL1B1degreepredictionedge.mat',...一般
% 'pic6_BA500_3t005nL1B1degreepredictionedge.mat',...一般
% 'pic8_BA500_4t005nL1B1degreepredictionedge.mat',...一般 中
% 'pic4_SW100_4_05t005nL1B1degreepredictionedge.mat',...
% 'pic4_SW100_4_08t005nL1B1degreepredictionedge.mat',...
%'pic4_SW100_4_08t01nL1B1degreepredictionedge.mat' 一般
% 'pic4_SW100_4_08t005nL02B1degreepredictionedge.mat',...
% 'pic4_SW100_4_08t02nL02B1degreepredictionedge.mat',...
% 'pic6_SW500_6_08t005nL1B1degreepredictionedge.mat',...
% 'pic4_SW100_4_08t01nL1B1degreepredictionedge.mat',...
% 'pic4_SW500_4_08t01nL1B1degreepredictionedge.mat',...
% 'pic6_SW500_6_08t01nL1B1degreepredictionedge.mat',...
% 'pic8_SW500_8_05t01nL1B1degreepredictionedge.mat',...
% 'pic6_BA500_3t01nL1B1degreepredictionedge.mat',...
% 'pic4_BA100_2t01nL1B1degreepredictionedge.mat') %一般
%'pic4_BA100_2t005nL02B1degreepredictionedge.mat' %很好
% aa=[ 'pic8_BA500_4t008nL1B1degreepredictionedge.mat'];
% aa=['pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat']

% aa=['pic8_BA50_4t008nL1B1betweenpredictionedge.mat']10 09
% aa=['pic8_SW500_8_05t008nL1B1degreepredictionedgeEsasort.mat'] 

% aa=['pic8_SW500_8_05t010nL1B1degreepredictionedgeEsasort.mat'] 
% aa=['pic8_BA500_4t010nL1B1degreepredictionedgeEsasort.mat']
% aa=['pic8_BA500_4t005nL1B1degreepredictionedgeEsasort.mat']
% aa=['pic8_BA500_4t012nL1B1degreepredictionedgeEsasort.mat']
% aa=['ADD_Edge_betweenness_pic8_BA500_4t008nL1B1degreepredictionedgeEsasort.mat']
aa=['ADD_Edge_betweenness_pic6_BA500_3t008nL1B1degreepredictionedgeEsasort.mat']
for ii=1:1:1
cascLfailtheta1=load(aa,'cascLfailtheta');
cascLfailtheta=cascLfailtheta1.cascLfailtheta;
ss=tabulate([cascLfailtheta(:,1);cascLfailtheta(:,2)]);
ss1=sortrows(ss,-3);
ss1(find(ss1(:,3)==0),:)=[];
ss=ss1(:,1);

% ssA=tabulate([cascLfailtheta11(1:40,1);cascLfailtheta11(1:40,2)]);
% ss1A=sortrows(ssA,-3);
% ss1A(find(ss1A(:,3)==0),:)=[];
% ssA=ss1A(:,1);




bb=textread('F:\py\Mtb\6_BA500_3\nx.closeness_centrality.txt','' , 'headerlines', 0);% 中介度
bb=[[1:size(bb,1)]',bb];
bb_c=bb;
bbp=bb(:,2);
bbp(find(bbp>=0.002))=100*bbp(find(bbp>=0.002));
bbp(find(bbp<0.01))=1*bbp(find(bbp<0.01));
bb=sortrows(bb,-2);

bc=textread('F:\py\Mtb\6_BA500_3\nx.betweenness_centrality.txt','' , 'headerlines', 0);% 中介度
bc=[[1:size(bc,1)]',bc];
bc_c=bc;
bcp=bc(:,2);
bc=sortrows(bc,-2);
bd=textread('F:\py\Mtb\6_BA500_3\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
bd=[[1:size(bd,1)]',bd];
bd_c=bd;
bdp=bd(:,2);
bd=sortrows(bd,-2);
avg_centrailty=[bd_c(:,1),(1/5.*bc_c(:,2)+1/5.*bd_c(:,2)+3/5.*bb_c(:,2))];
% avg_centrailty(:,1)=floor(avg_centrailty(:,1))
avg_centrailty=sortrows(avg_centrailty,-2);

% be=textread('F:\py\Mtb\8_BA50_4\nx.eigenvector_centrality.txt','' , 'headerlines', 0);% 中介度
% be=[[1:size(be,1)]',be];
% bep=be(:,2);
% be=sortrows(be,-2);
% figure(ii+2)
% plot([1:size(be,1)]',bb(:,1),[1:size(be,1)]',bc(:,1),[1:size(be,1)]',bd(:,1),[1:size(be,1)]',be(:,1))
% legend({'betweenness centrality','closeness centrality','degree centrality','eigenvector centrality'},'FontSize',12);

[resumbb,rerowsbb]=ismember(ss,bb(:,1),'rows');
[resumbc,rerowsbc]=ismember(ss,bc(:,1),'rows');
[resumbd,rerowsbd]=ismember(ss,bd(:,1),'rows');
[resumavg,rerowsavg]=ismember(ss,avg_centrailty(:,1),'rows');

% [resumbe,rerowsbe]=ismember(ss,be(:,1),'rows');
rss=[1:size(ss,1)]';
bbcd=[rss,rerowsbc,rerowsbd,rerowsavg];

% avg1=[1/3*rerowsbd+1/3*rerowsbc+1/3*rerowsbb];
figure(ii+2)
% plot(rss,rss ,'m',rss,rerowsbc,' bo--',rss ,rerowsbd,' g*--',rss,rerowsavg,'r*--')%(1:10,1)
plot(rss,rss ,'m-','linewidth',1.5)
hold on
plot(rss ,rerowsbb,' go--',rss,rerowsbc,' b.:',rss ,rerowsbd,' k.--',rss,rerowsavg,'r*--')
axis([0 90 0 205])   
set(gca,'xtick',0:5:90)
set(gca,'ytick',0:20:205)
legend({'Positive proportion function','Closeness centrality','Betweenness centrality','Degree centrality','Comprehensive Factor'},'FontSize',12);

% hold on 
% plot(rss([1:2,4:6]),rss([1:2,4:6]),'ko','Markersize',12)
% hold on
% plot(rss([1:2,4:6]),rss([1:2,4:6]),'ko','Markersize',13)
xlabel('Sorted position 1 of nodes in special links set of network');
ylabel('Node Sorted position');
figure(ii+1)
% plot(rss ,rss ,'m',rss ,rerowsbb,' r*--',rss ,rerowsbc,' go--',rss ,rerowsbd,' b.--',rss,rerowsbe,'cd--')%(1:10,1)
% legend({'linear function','betweenness centrality','closeness centrality','degree centrality','eigenvector centrality'},'FontSize',12);
% plot(rss ,rss ,'m',rss ,rerowsbb,' r*--',rss ,rerowsbc,' go--',rss ,rerowsbd,' b.--',rss ,avg1,'ko-')%(1:10,1)
% plot(rss ,rss ,'m',rss ,rerowsbc,' go--',rss ,rerowsbd,' k*--',rss,rerowsavg,'bo:','linewidth',1)%(1:10,1)
% rss,rss ,'m',rss ,rerowsbb,' go:',rss,rerowsbc,' b.--',rss ,rerowsbd,' r.--',rss,rerowsavg,'ko--','linewidth',1)%(1:10,1)

% plot(rss ,rerowsbb,' go--',rss,rerowsbc,' b.:',rss ,rerowsbd,' k.--',rss,rerowsavg,'r*--')

% plot(rss ,rerowsbb,' k.--',rss,rerowsbc,' bo:',rss ,rerowsbd,' g*--',rss,rerowsavg,'r*--')
plot(rerowsavg-rss,'k.--')
hold on
dd=rerowsavg-rss
plot(rss([1:17,19:26]), dd([1:17,19:26]),'rh','Markersize',13,'linewidth',1.5)
hold on
axis([0 30 -15 35])   
set(gca,'xtick',0:2:30)
set(gca,'ytick',-15:5:35)

plot(0.*rss ,'m--','linewidth',1.5)
legend({'Comprehensive Factor'},'FontSize',12);
% legend({'Positive proportion function','Closeness centrality','Betweenness centrality','Degree centrality','Comprehensive Factor'},'FontSize',12);
% title('8_BA100_4 BA network important nodes prediction','FontSize',12)
xlabel('Sorted position 1 of nodes in special links set of network');
ylabel({'Deviation of sorted position','of nodes in special links set of network'});
figure(ii+3)
plot(rerowsbb-rss,'gs:')
hold on
plot(rerowsbc-rss,'cp-.')
hold on
plot(rerowsbd-rss,'bo--')
hold on
plot(rerowsavg-rss,'k.--')
hold on
plot(0.*rerowsavg,'m--')
legend({'Betweenness centrality','Closeness centrality','Degree centrality','Comprehensive Factor'},'FontSize',12);
xlabel('Sorted position 1 of nodes in special links set of network');
ylabel({'Deviation of sorted position','of nodes in special links set of network'});
axis([0 55 -35 205])   
set(gca,'xtick',0:5:55)
set(gca,'ytick',-35:20:205)


figure(ii+8)
bar(rerowsavg-rss)
axis([0 50 -35 205])   
set(gca,'xtick',0:5:50)
set(gca,'ytick',-35:10:205)
xlabel('Sorted position 1 of nodes in special links set of network');
ylabel({'Deviation of sorted position','of nodes in special links set of network'});
figure(ii+7)
bar(rerowsbd-rss)
axis([0 50 -35 205])   
set(gca,'xtick',0:5:50)
set(gca,'ytick',-35:10:205)
figure(ii+9)
bar(rerowsbb-rss)
axis([0 50 -35 205])   
set(gca,'xtick',0:5:50)
set(gca,'ytick',-35:10:205)
end




%strcat(aa7(ii,(4:40)))