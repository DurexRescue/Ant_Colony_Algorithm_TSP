%%%һ����������Ҫ�ݷ�ȫ��31��ʡ����У���Ҫѡ����̵�·��%%%%

%%%��Ⱥ�㷨���TSP����%%%%%%%
 
clear all; %������б���
close all; %��ͼ
clc;       %����

m = 50;       %% m ���ϸ���
Alpha = 1;    %% Alpha ������Ϣ����Ҫ�̶ȵĲ���
Beta = 5;     %% Beta ��������ʽ������Ҫ�̶ȵĲ���
Rho = 0.1;    %% Rho ��Ϣ������ϵ��
NC_max = 100; %%����������
Q = 100;      %%��Ϣ������ǿ��ϵ��,Ҳ������ѭ��һ�����ͷŵ���Ϣ������

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
];%%31��ʡ������
%%-------------------------------------------------------------------------
%% ��Ҫ����˵��
%% C n�����е����꣬n��2�ľ���
%% NC_max ����������
%% m ���ϸ���
%% Alpha ������Ϣ����Ҫ�̶ȵĲ���
%% Beta ��������ʽ������Ҫ�̶ȵĲ���
%% Rho ��Ϣ������ϵ��
%% Q ��Ϣ������ǿ��ϵ��
%% Route_best �������·��
%% L_best �������·�ߵĳ���
%%=========================================================================
%% ѭ������
%control = 1;
%while control <= 30
    %% ��һ����������ʼ��
    n = size(C,1);%n��ʾ����Ĺ�ģ�����и�����
    Distance = zeros(n,n);%Distance��ʾ��ȫͼ�ĸ�Ȩ�ڽӾ���,��ʼ��Ϊȫ0

    for i = 1 : n
        for j = 1 : n
            if i ~= j
                Distance(i, j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
            else
                Distance(i, j) = eps; %i=jʱ�����㣬Ӧ��Ϊ0�����������������Ҫȡ��������eps��������Ծ��ȣ���ʾ
            end
            Distance(j, i) = Distance(i, j);	%����Գ�    
        end
    end

    Eta = 1./Distance;          %EtaΪ�������ӣ�������Ϊ����ĵ���
    Tau = ones(n,n);     %TauΪ��Ϣ�ؾ���,ȫ����ʼ��Ϊ1
    Tabu = zeros(m,n);   %�洢����¼·��������
    NC = 1;              %��������������¼��������
    Route_best = zeros(NC_max,n);       %�������·��
    L_best = inf.*ones(NC_max,1);   %�������·�ߵĳ���
    L_ave = zeros(NC_max,1);        %����·�ߵ�ƽ������

    while NC <= NC_max  %ֹͣ����֮һ���ﵽ������������ֹͣ
        %%�ڶ����������mֻ���Ϸŵ�n��������
        Randpos = [];   %�漴��ȡ
        for i = 1 : (ceil(m/n))
            Randpos = [Randpos, randperm(n)]; %ִ�����β���֮��Randpos��ά��Ϊ1*62���ֿ�����ԭ������ȡֵ��Χ��Ҫ��1-31
        end
        Tabu(:,1) = (Randpos(1,1:m))';%���д洢ÿֻ���ϵķ��ʵ���Ϣ

        %%��������mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
        for j = 2 : n     %���ڳ��в����㣬����ѭ����������-1����
            for i = 1 : m
                visited = Tabu(i, 1 : (j - 1));  %��¼�ѷ��ʵĳ��У������ظ����ʣ���iֻ���ϵķ��ʼ�¼���洢�ڵ�i�У�1��j-1��,ÿ�����Tabu����ȡ��
                J = zeros(1, (n - j + 1));       %�洢�����ʵĳ��У���ʼ��Ϊ1��n-j+1�еľ���n-j+1��ʾʣ�»�δ�����ʵĳ��и���
                P = J;                           %�����ʳ��е�ѡ����ʷֲ�
                Jc = 1;                          %һ��������ֵ
                for k = 1 : n                            %�ҵ�δ���ʵĳ��У�����������J��
                    if isempty(find(visited == k, 1))  %��ʼʱ��0 find����������visited������k���ڵ�λ�� û���򷵻�0 1��ʾֻ��1��
                        J(Jc) = k;
                        Jc = Jc + 1;                     %���ʵĳ��и����Լ�1
                    end
                end
                %��������ѡ���еĸ��ʷֲ�
                for k = 1 : length(J)
                    P(k) = (Tau(visited(end),J(k))^Alpha) * (Eta(visited(end),J(k))^Beta);  %����ÿһ��ķ��ӣ���ĸ��Ϊ���еķ��ӵĺ�
                end
                P = P / (sum(P));           %���´����ʳ��и���������Ԫ�ص�ֵ
                %�����̶ķ�ѡȡ��һ������
                Pcum = cumsum(P);           %cumsum��Ԫ�ص�����ۼӺ�,����ֵΪ��Pά����ͬ���о���
                Select = find(Pcum >= rand);%ѡ�������Խϴ���Ǹ��ڵ�
                to_visit = J(Select(1));    %���������ʵĵ�    
                Tabu(i,j) = to_visit;       %���ʸõ㣬�����¼��
            end
        end

        if NC >= 2
            Tabu(1,:) = Route_best(NC-1,:);   %����һ���ϴ�����·����tabu��һ�У����ϱ��ε������������̫��
        end

        %%���Ĳ�����¼���ε���ÿֻ�������߾���L����¼ÿ�ε������·�߾���L_best�����·����ϢR_best
        L = zeros(m,1);                       %��ʼ����Ϊ0��m*1��������
        for i = 1 : m
            Route = Tabu(i, : );
            for j = 1 : (n - 1)
                L(i) = L(i) + Distance(Route(j), Route(j + 1));   %ԭ������ϵ�j�����е���j+1�����еľ���
            end
            L(i) = L(i) + Distance(Route(1), Route(n));         %һ���������߹��ľ���
        end

        L_best(NC) = min(L);                  %��Ѿ���ȡ��С
        L_ave(NC) = mean(L);                  %���ֵ������ƽ������

        pos = find(L == L_best(NC));
        Route_best(NC, : ) = Tabu(pos(1), : );        %���ֵ���������·��

        NC = NC + 1;                          %��������



        %%���岽��������Ϣ��
        Delta_Tau = zeros(n,n);               %��ʼʱ��Ϣ��Ϊn*n��0����
        
        %---------------------------------------------------------------------------------
        %�Ľ�ǰ
%         for i=1:m
%             for j=1:(n-1)
%                 Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
%                 %�˴�ѭ����·����i��j���ϵ���Ϣ������
%             end
%             Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
%             %.�˴�ѭ��������·���ϵ���Ϣ������
%         end
        %---------------------------------------------------------------------------------
        %%�Ż��Ľ���ȡ·��������ǰ�ٷ�֮��ʮ��·����������Ϣ�أ�����·�߲�������Ϣ��
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
                %�˴�ѭ����·����i��j���ϵ���Ϣ������
            end
            Delta_Tau(Update_Info(i,n),Update_Info(i,1))=Delta_Tau(Update_Info(i,n),Update_Info(i,1)) + part * Q/L(i);
            %.�˴�ѭ��������·���ϵ���Ϣ������
        end


        Tau=(1-Rho).*Tau+Delta_Tau; %������Ϣ�ػӷ������º����Ϣ��

        %%�����������ɱ�����
        Tabu=zeros(m,n);            %ֱ������������
        
        %%ͼ��ʵʱ����
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


    
    
    %%���߲���������
    Pos=find(L_best==min(L_best)); %�ҵ����·������0Ϊ�棩
    Shortest_Route=Route_best(Pos(1),:);%���������������·��
    Shortest_Length=L_best(Pos(1)) %��������������̾���
    %% ѭ������
    %control = control + 1;
    %length_each = [length_each ; Shortest_Length];%����30�ε�ʱ��洢���·�����ȵľ���
%end
%% ����·��ͼ����L_best��L_ave��������
figure(1)
subplot(1,2,1)                  %���Ƶ�һ����ͼ��
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
title('�����������Ż���� ')

subplot(1,2,2)                  %���Ƶڶ�����ͼ��
plot(L_best)
hold on                         %����ͼ��
plot(L_ave,'r')
title('ƽ���������̾���')      %����