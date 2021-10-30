%%%一个旅行商人要拜访全国31个省会城市，需要选择最短的路径%%%%

%%%蚁群算法解决TSP问题%%%%%%%
 
clear all; %清除所有变量
close all; %清图
clc;       %清屏

m = 50;       %% m 蚂蚁个数
Alpha = 1;    %% Alpha 表征信息素重要程度的参数
Beta = 5;     %% Beta 表征启发式因子重要程度的参数
Rho = 0.1;    %% Rho 信息素蒸发系数
NC_max = 100; %%最大迭代次数
Q = 100;      %%信息素增加强度系数,也即蚂蚁循环一次所释放的信息素总量

length_each = [];

C=[
1304 2312;
3639 1315;
4177 2244;
3712 1399;
3488 1535;
3326 1556;
3238 1229;
4196 1004;
4312 790;
4386 570;
3007 1970;
2562 1756;
2788 1491;
2381 1676;
1332 695;
3715 1678;
3918 2179;
4061 2370;
3780 2212;
3676 2578;
4029 2838;
4263 2931;
3429 1908;
3507 2367;
3394 2643;
3439 3201;
2935 3240;
3140 3550;
2545 2357;
2778 2826;
2370 2975
];%%31个省会坐标
%%-------------------------------------------------------------------------
%% 主要符号说明
%% C n个城市的坐标，n×2的矩阵
%% NC_max 最大迭代次数
%% m 蚂蚁个数
%% Alpha 表征信息素重要程度的参数
%% Beta 表征启发式因子重要程度的参数
%% Rho 信息素蒸发系数
%% Q 信息素增加强度系数
%% Route_best 各代最佳路线
%% L_best 各代最佳路线的长度
%%=========================================================================
%% 循环控制
%control = 1;
%while control <= 30
    %% 第一步：变量初始化
    n = size(C,1);%n表示问题的规模（城市个数）
    Distance = zeros(n,n);%Distance表示完全图的赋权邻接矩阵,初始化为全0

    for i = 1 : n
        for j = 1 : n
            if i ~= j
                Distance(i, j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
            else
                Distance(i, j) = eps; %i=j时不计算，应该为0，但后面的启发因子要取倒数，用eps（浮点相对精度）表示
            end
            Distance(j, i) = Distance(i, j);	%矩阵对称    
        end
    end

    Eta = 1./Distance;          %Eta为启发因子，这里设为距离的倒数
    Tau = ones(n,n);     %Tau为信息素矩阵,全部初始化为1
    Tabu = zeros(m,n);   %存储并记录路径的生成
    NC = 1;              %迭代计数器，记录迭代次数
    Route_best = zeros(NC_max,n);       %各代最佳路线
    L_best = inf.*ones(NC_max,1);   %各代最佳路线的长度
    L_ave = zeros(NC_max,1);        %各代路线的平均长度

    while NC <= NC_max  %停止条件之一：达到最大迭代次数，停止
        %%第二步：随机将m只蚂蚁放到n个城市上
        Randpos = [];   %随即存取
        for i = 1 : (ceil(m/n))
            Randpos = [Randpos, randperm(n)]; %执行两次操作之后Randpos的维度为1*62，分开做的原因在于取值范围需要在1-31
        end
        Tabu(:,1) = (Randpos(1,1:m))';%按行存储每只蚂蚁的访问点信息

        %%第三步：m只蚂蚁按概率函数选择下一座城市，完成各自的周游
        for j = 2 : n     %所在城市不计算，共计循环（城市数-1）次
            for i = 1 : m
                visited = Tabu(i, 1 : (j - 1));  %记录已访问的城市，避免重复访问，第i只蚂蚁的访问记录被存储在第i行，1～j-1列,每次需从Tabu表中取出
                J = zeros(1, (n - j + 1));       %存储待访问的城市，初始化为1行n-j+1列的矩阵，n-j+1表示剩下还未被访问的城市个数
                P = J;                           %待访问城市的选择概率分布
                Jc = 1;                          %一个列索引值
                for k = 1 : n                            %找到未访问的城市，并存在数组J中
                    if isempty(find(visited == k, 1))  %开始时置0 find函数返回在visited数组中k所在的位置 没有则返回0 1表示只找1次
                        J(Jc) = k;
                        Jc = Jc + 1;                     %访问的城市个数自加1
                    end
                end
                %下面计算待选城市的概率分布
                for k = 1 : length(J)
                    P(k) = (Tau(visited(end),J(k))^Alpha) * (Eta(visited(end),J(k))^Beta);  %计算每一项的分子，分母即为所有的分子的和
                end
                P = P / (sum(P));           %更新待访问城市概率数组中元素的值
                %按轮盘赌法选取下一个城市
                Pcum = cumsum(P);           %cumsum，元素的逐次累加和,返回值为和P维度相同的行矩阵
                Select = find(Pcum >= rand);%选择概率相对较大的那个节点
                to_visit = J(Select(1));    %即将被访问的点    
                Tabu(i,j) = to_visit;       %访问该点，存进记录中
            end
        end

        if NC >= 2
            Tabu(1,:) = Route_best(NC-1,:);   %保留一下上次最优路线至tabu第一行，保障本次迭代情况不至于太差
        end

        %%第四步：记录本次迭代每只蚂蚁所走距离L，记录每次迭代最佳路线距离L_best和最佳路线信息R_best
        L = zeros(m,1);                       %开始距离为0，m*1的列向量
        for i = 1 : m
            Route = Tabu(i, : );
            for j = 1 : (n - 1)
                L(i) = L(i) + Distance(Route(j), Route(j + 1));   %原距离加上第j个城市到第j+1个城市的距离
            end
            L(i) = L(i) + Distance(Route(1), Route(n));         %一轮下来后走过的距离
        end

        L_best(NC) = min(L);                  %最佳距离取最小
        L_ave(NC) = mean(L);                  %此轮迭代后的平均距离

        pos = find(L == L_best(NC));
        Route_best(NC, : ) = Tabu(pos(1), : );        %此轮迭代后的最佳路线

        NC = NC + 1;                          %迭代继续



        %%第五步：更新信息素
        Delta_Tau = zeros(n,n);               %开始时信息素为n*n的0矩阵
        
        %---------------------------------------------------------------------------------
        %改进前
%         for i=1:m
%             for j=1:(n-1)
%                 Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
%                 %此次循环在路径（i，j）上的信息素增量
%             end
%             Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
%             %.此次循环在整个路径上的信息素增量
%         end
        %---------------------------------------------------------------------------------
        %%优化改进，取路径长度在前百分之二十的路线来更新信息素，其余路线不产生信息素
        Update_Info = [];
        quantity = ceil(m * 0.2);
        [L_seq, L_pos] = sort(L);
        L_pos_0X = L_pos(1 : quantity, : );
        for i = 1 : length(L_pos_0X)
            Update_Info = [Update_Info ;Tabu(i, : )];
        end

        for i = 1 : quantity
            if i <= quantity * 0.5
            part = 1;
            elseif i > quantity * 0.5 && i <= quantity *0.7
                part = 0.5;
            elseif i > quantity * 0.7
            part = 0.2;
            end
            for j = 1 : (n - 1)
                Delta_Tau(Update_Info(i,j),Update_Info(i,j+1))=Delta_Tau(Update_Info(i,j),Update_Info(i,j+1)) + part * Q/L(i);
                %此次循环在路径（i，j）上的信息素增量
            end
            Delta_Tau(Update_Info(i,n),Update_Info(i,1))=Delta_Tau(Update_Info(i,n),Update_Info(i,1)) + part * Q/L(i);
            %.此次循环在整个路径上的信息素增量
        end


        Tau=(1-Rho).*Tau+Delta_Tau; %考虑信息素挥发，更新后的信息素

        %%第六步：禁忌表清零
        Tabu=zeros(m,n);            %直到最大迭代次数
        
        %%图像实时更新
        figure(1)
        clf
        N=size(Route_best,2);
        scatter(C(:,1),C(:,2));
        for i = 1:size(C,1)
            text(C(i,1),C(i,2),['   ' num2str(i)]);
        end
        hold on
        plot([C(Route_best(NC-1,1),1),C(Route_best(NC-1,N),1)],[C(Route_best(NC-1,1),2),C(Route_best(NC-1,N),2)],'g')
        hold on
        for ii=2:N
            plot([C(Route_best(NC-1,ii-1),1),C(Route_best(NC-1,ii),1)],[C(Route_best(NC-1,ii-1),2),C(Route_best(NC-1,ii),2)],'g')
            hold on
        end
        title(['Iteration:',num2str(NC-1),'Length:',num2str(L_best(NC-1))])
        
    end


    
    
    %%第七步：输出结果
    Pos=find(L_best==min(L_best)); %找到最佳路径（非0为真）
    Shortest_Route=Route_best(Pos(1),:);%最大迭代次数后最佳路径
    Shortest_Length=L_best(Pos(1)) %最大迭代次数后最短距离
    %% 循环控制
    %control = control + 1;
    %length_each = [length_each ; Shortest_Length];%运行30次的时候存储最短路径长度的矩阵
%end
%% 画出路线图，和L_best，L_ave迭代曲线
figure(1)
subplot(1,2,1)                  %绘制第一个子图形
N=length(Shortest_Route);
scatter(C(:,1),C(:,2));
for i = 1:size(C,1)
	text(C(i,1),C(i,2),['   ' num2str(i)]);
end
hold on
plot([C(Shortest_Route(1),1),C(Shortest_Route(N),1)],[C(Shortest_Route(1),2),C(Shortest_Route(N),2)],'g')
hold on
for ii=2:N
    plot([C(Shortest_Route(ii-1),1),C(Shortest_Route(ii),1)],[C(Shortest_Route(ii-1),2),C(Shortest_Route(ii),2)],'g')
    hold on
end
title('旅行商问题优化结果 ')

subplot(1,2,2)                  %绘制第二个子图形
plot(L_best)
hold on                         %保持图形
plot(L_ave,'r')
title('平均距离和最短距离')      %标题