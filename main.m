
CostFunction=@(x) F17(x);       

nVar=2;            

VarSize=[1 nVar];   %

VarMin=-100;        
VarMax= 100;         

%% ABC Settings

MaxIt=50;             

nPop=500;              

nOnlooker=nPop;        

L=round(0.6*nVar*nPop); 

a=1;                   

%%
%改进ABC
% 参数设置
D = 2;          % 优化问题的维度
N = 500;          % 种群大小
MaxDT = 50;    % 最大迭代次数

%% Initialization

% Empty Bee Structure
empty_bee.Position=[];
empty_bee.Cost=[];

% Initialize Population Array
pop=repmat(empty_bee,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Population
for i=1:nPop
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Cost=CostFunction(pop(i).Position);
    if pop(i).Cost<=BestSol.Cost
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% ABC Main Loop

for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    MeanCost = mean([pop.Cost]);
    for i=1:nPop
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        k=K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        newbee.Cost=CostFunction(newbee.Position);
        
        % Comparision
        if newbee.Cost<=pop(i).Cost
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
            pop(i).Cost=CostFunction(pop(i).Position);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if pop(i).Cost<=BestSol.Cost
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Display Iteration Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end
    
%% Results


%%%%%%%%%%%%%



U = 100 * ones(1, D);
L = -100 * ones(1, D);
limit = 0.8 * N * D;

% 初始化数组和变量
trial = zeros(1, N);
MR_max = 0.9;
fitness_values = zeros(1, MaxDT);

tic;

% 初始化种群
X = (U - L) .* rand(N, D) + L;
val = zeros(1, N);
for i = 1:N
    val(i) = F17(X(i, :)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% 计算初始全局最优
[gscore, Index] = min(val);
gsolution = X(Index, :);

position = zeros(1, N);
Improveposition = zeros(1, N);
Probability = zeros(1, N);

% 主循环
for iter = 1:MaxDT
    lamta = (MaxDT - iter + 1) / MaxDT;
    MR = exp(-iter / MaxDT) * MR_max;

    % 迭代更新种群
    for i = 1:N
        sol = X(i, :);
        neighbor = fix(rand * N) + 1;
        while neighbor == i
            neighbor = fix(rand * N) + 1;
        end

        position(neighbor) = position(neighbor) + 1;

        for k = 1:D
            if rand < MR
                sol(k) = X(i, k) + (rand - 0.5) * 2 * (X(i, k) - X(neighbor, k));
            end
        end

        Funval = F17(sol);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if Funval < val(i)
            X(i, :) = sol;
            val(i) = Funval;
            Improveposition(neighbor) = Improveposition(neighbor) + 1;
        else
            trial(i) = trial(i) + 1;
        end
    end

    ind = find(val == min(val));
    ind = ind(end);
    gsolution = X(ind, :);
    gscore = val(ind);
    fprintf('第 %d 轮的适应度值：%f\n', iter, gscore);

    % 计算每个个体的选择概率
    for i = 1:N
        Probability(i) = Improveposition(i) / (position(i) + 0.01);
    end

    if sum(Probability == 0) == N
        for i = 1:N
            Probability(i) = 1 / N;
        end
    end

    % 轮盘赌选择
    t = 0;
    i = 1;

    while (t < N)
        if (rand < Probability(i))
            t = t + 1;
            sol = X(i, :);
            neighbor = fix(rand * N) + 1;
            while neighbor == i
                neighbor = fix(rand * N) + 1;
            end

            if val(neighbor) < val(i)
                for k = 1:D
                    if rand < MR
                        if X(neighbor, k) > X(i, k)
                            sol(k) = (1 - lamta) * X(i, k) + lamta * rand * (X(neighbor, k) - X(i, k));
                        else
                            sol(k) = (1 - lamta) * X(i, k) - lamta * rand * (X(i, k) - X(neighbor, k));
                        end
                    end
                end
            else
                for k = 1:D
                    if rand < MR
                        sol(k) = (1 - lamta) * gsolution(k) + lamta * (rand - 0.5) * 2 * (gsolution(k) - X(i, k));
                    end
                end
            end

            Funval = F17(sol); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if Funval < val(i)
                X(i, :) = sol;
                val(i) = Funval;
                Improveposition(neighbor) = Improveposition(neighbor) + 1;
            else
                trial(i) = trial(i) + 1;
            end
        end

        i = i + 1;
        if (i == N + 1)
            i = 1;
        end
    end

    % 处理被丢弃的个体
    ind = find(trial == max(trial));
    ind = ind(end);

    if (trial(ind) > limit)
        trial(ind) = 0;
        sol = (1 - lamta) * gsolution + lamta * (rand - 0.5) * 2 * (gsolution - X(ind, :));
        X(ind, :) = sol;
        val(ind) = F17(sol); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    ind = find(val == min(val));
    ind = ind(end);

    if (val(ind) < gscore)
        gscore = val(ind);
        gsolution = X(ind, :);
    end
    
    fitness_values(iter) = gscore;
    

end

toc;
disp(['原始ABC最小值为：' num2str(min(BestCost))]);
disp(['改进ABC目标函数的最小值为' num2str(gscore)]);
gsolution;


% 绘制迭代次数与适应度图

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);

grid on;
hold on

plot(1:MaxDT, fitness_values, 'r', 'LineWidth', 2);
title('迭代次数与适应度');
xlabel('迭代次数');
ylabel('适应度');
legend('原始ABC','改进ABC')
grid on;

