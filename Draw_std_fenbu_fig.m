clear
figure(1)
stddata=[1900.931958
1832.01663
1592.896365
1222.066137
1179.85312
1106.175373
741.0711135
401.334367
380.7048016
318.4205987
297.3748936
284.2424283
260.1665215
172.9562109
156.1085905
153.3376179
90.12045434
76.08711493
60.25589374
42.92
41.32
38.46
37.14
36.34378244
35.1702971
10.02
8.211
7.5762
0;]'
stddata(2,:)=zeros(1,length(stddata))
legenddata={'1.4_BA50_2', '2.6_BA50_3',  '3.8_BA50_4','4.6_BA100_3','5.4_BA100_2','6.8_BA100_4',...
        '7.20_BA100_11','8.4_SW100_4_05' ,'9.6_SW100_6_08','10.8_SW100_8_08','11.4_BA500_2','12.8_BA500_4',...
'13.6_BA500_3','14.4_BA1000_2','15.8_BA1000_4','16.6_BA1000_3','17.4_SW500_4_08','18.6_SW500_6_08','19.8_SW500_8_05','20.4_SW1000_4_05',...
'21.4_BA5000_2','22.8_BA5000_4','23.6_SW1000_6_05','24.6_BA5000_3','25.8_SW1000_8_08','26.4_SW5000_4_05','27.6_SW5000_6_05','28.8_SW5000_8_08','29.6_SW500_6_00'}
hh=bar(stddata)
color1=jet
color11=color1(1:length(stddata),:)
for i=1:length(stddata)
 set(hh(i),'FaceColor',color11(i,:))  
end


xlimdown=0.6;
xlimup=1.4;
xlim([xlimdown xlimup])
% xlable(legenddata)1/2*(xlimup-xlimdown)/length(stddata)
set(gca,'XTickLabel',legenddata)
xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄 
xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。                    
ytextp=-0.1*yt(3)*ones(1,length(xtb));
% xx=[0.6:2/55:1.4]
legenddatacell={}
for ii=1:length(xtb)
   legenddatacell{ii}=num2str(ii);
end
text([xlimdown+1/2*(xlimup-xlimdown)/length(stddata):(xlimup-xlimdown)/length(stddata):xlimup-1/2*(xlimup-xlimdown)/length(stddata)],ytextp,legenddatacell,'HorizontalAlignment','right','rotation',46)
set(gca,'XTickLabel',[]); %将原坐标(1,2,3,..)去掉
xlabel('Network type')
ylabel('Standard deviation of edge initial load')
legend(legenddata)
% clear
% y3=[4,2,5,5,5,5,8,8,5,8];
% x3str={'<=10','10-25','25-50','50-90','90-300','300-350','350-450','450-550','550-700','>700'}; %新坐标的值
% bar(y3)  %先 bar 后 set
% xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
% xt = get(gca,'XTick');% 获取横坐标轴刻度句柄
% yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄 
% xtextp=xt;%每个标签放置位置的横坐标，这个自然应该和原来的一样了。                    
% ytextp=-0.1*yt(3)*ones(1,length(xt));
% text(xtextp,ytextp,x3str,'HorizontalAlignment','right','rotation',46)
% set(gca,'XTickLabel',[]); %将原坐标(1,2,3,..)去掉
% % --添加坐标提示
% xlabel('Rating Counts');
% ylabel('Number of Ratings'); 
% legend('FilmTrust');
% % ----设置 xlabel在右边而非中间
% h=xlabel('Rating Counts');
% xlim = get(gca,'XLim');
% ylim = get(gca,'YLim');
% set(h,'Position',[xlim(2)+(xlim(2)-xlim(1))*0.05,ylim(1)])
% x=[1.2 3.1 2.2; 0 0 0];
% 
% subplot(2,2,1);
% bar(x);
% xlim([0 2]);% another way: axis([0 2 0 4])
% 
% subplot(2,2,2);
% bar(x);
% xlim([0 2]);
% set(gca,'xticklabel',{'A   B   C',''});
% 
% subplot(2,2,3);
% bar(x);
% xlim([0 2])
% legend('A','B','C','0')
% 
% subplot(2,2,4)
% x0=[1.2 3.1 2.2]
% bar(x0)











