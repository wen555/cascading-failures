B=textread('6_BA500_3\nx.degree_centrality.txt','' , 'headerlines', 0);% 中介度
D=textread('F:\py\Mtb\6_BA500_3\nx.betweenness_centrality.txt','', 'headerlines', 0);% DU
L=textread('F:\py\Mtb\6_BA500_3\nx.closeness_centrality.txt','', 'headerlines', 0);%%连接关系
plot(B,'r')
hold on
plot(D,'g')
hold on
plot(L,'k')
figure(2)
plot(B,D)