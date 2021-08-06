function varargout = mtsp_ga(xy,dmat,salesmen,min_tour,pop_size,num_iter,show_prog,show_res)
% 输入：
%     XY：各个城市坐标的N*2矩阵，N为城市的个数
%     DMAT：各个城市之间的距离矩阵
%     SALESMEN ：旅行商人数
%     MIN_TOUR：每个人所经过的最少城市点
%     POP_SIZE：种群大小
%     NUM_ITER：迭代次数
% 输出：
%     最优路线
%     总距离
 % 初始化
nargs = 8;
for k = nargin:nargs-1
    switch k
        case 0
            xy = 10*rand(40,2);
        case 1
            N = size(xy,1);
            a = meshgrid(1:N);
            dmat = reshape(sqrt(sum((xy(a,:)-xy(a',:)).^2,2)),N,N);
        case 2
            salesmen = 5; %     SALESMEN：人数
        case 3
            min_tour = 3; %     MIN_TOUR :每个人所经过的最少城市点
        case 4
            pop_size = 80; % POP_SIZE：种群大小
        case 5
            num_iter = 5e3; % NUM_ITER：迭代次数
        case 6
            show_prog = 1;
        case 7
            show_res = 1;
        otherwise
    end
end
%调整输入数据 
N = size(xy,1);
[nr,nc] = size(dmat);
if N ~= nr || N ~= nc
    error('Invalid XY or DMAT inputs!')
end
n = N;
 %约束条件
salesmen = max(1,min(n,round(real(salesmen(1)))));
min_tour = max(1,min(floor(n/salesmen),round(real(min_tour(1)))));
pop_size = max(8,8*ceil(pop_size(1)/8));
num_iter = max(1,round(real(num_iter(1))));
show_prog = logical(show_prog(1));
show_res = logical(show_res(1));
 % 路线的初始化
num_brks = salesmen-1;
dof = n - min_tour*salesmen;       
addto = ones(1,dof+1);
for k = 2:num_brks
    addto = cumsum(addto);
end
cum_prob = cumsum(addto)/sum(addto);
 % 种群初始化
pop_rte = zeros(pop_size,n);          
pop_brk = zeros(pop_size,num_brks);   
for k = 1:pop_size
    pop_rte(k,:) = randperm(n);
    pop_brk(k,:) = randbreaks();
end
tmp_pop_rte = zeros(8,n);
tmp_pop_brk = zeros(8,num_brks);
new_pop_rte = zeros(pop_size,n);
new_pop_brk = zeros(pop_size,num_brks);
 %图形中各个路线的颜色
clr = [1 0 0; 0 0 1; 0.67 0 1; 0 1 0; 1 0.5 0];
if salesmen > 5
    clr = hsv(salesmen);
end
 % 遗传算法实现
global_min = Inf;
total_dist = zeros(1,pop_size);
dist_history = zeros(1,num_iter);
if show_prog
    pfig = figure('Name','MTSP_GA | Current Best Solution','Numbertitle','off');
end
for iter = 1:num_iter
    for p = 1:pop_size
        d = 0;
        p_rte = pop_rte(p,:);
        p_brk = pop_brk(p,:);
        rng = [[1 p_brk+1];[p_brk n]]';
        for s = 1:salesmen
            d = d + dmat(p_rte(rng(s,2)),p_rte(rng(s,1)));
            for k = rng(s,1):rng(s,2)-1
                d = d + dmat(p_rte(k),p_rte(k+1));
            end
        end
        total_dist(p) = d;
    end
    % 种群中的最优路线
    [min_dist,index] = min(total_dist);
    dist_history(iter) = min_dist;
    if min_dist < global_min
        global_min = min_dist;
        opt_rte = pop_rte(index,:);
        opt_brk = pop_brk(index,:);
        rng = [[1 opt_brk+1];[opt_brk n]]';
        if show_prog
            % 画最优路线图
            figure(pfig);
            for s = 1:salesmen
                rte = opt_rte([rng(s,1):rng(s,2) rng(s,1)]);
                plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
                title(sprintf('Total Distance = %1.4f, Iterations = %d',min_dist,iter));
                hold on
            end
            hold off
        end
    end
   %遗传算法的进程操作 
    rand_grouping = randperm(pop_size);
    for p = 8:8:pop_size
        rtes = pop_rte(rand_grouping(p-7:p),:);
        brks = pop_brk(rand_grouping(p-7:p),:);
        dists = total_dist(rand_grouping(p-7:p));
        [ignore,idx] = min(dists);
        best_of_8_rte = rtes(idx,:);
        best_of_8_brk = brks(idx,:);
        rte_ins_pts = sort(ceil(n*rand(1,2)));
        I = rte_ins_pts(1);
        J = rte_ins_pts(2);
        for k = 1:8 % 产生新的解
            tmp_pop_rte(k,:) = best_of_8_rte;
            tmp_pop_brk(k,:) = best_of_8_brk;
            switch k
                case 2 % 选择交叉算子
                    tmp_pop_rte(k,I:J) = fliplr(tmp_pop_rte(k,I:J));
                case 3 %变异算子
                    tmp_pop_rte(k,[I J]) = tmp_pop_rte(k,[J I]);
                case 4 %重排序
                    tmp_pop_rte(k,I:J) = tmp_pop_rte(k,[I+1:J I]);
                case 5 %修改城市点
                    tmp_pop_brk(k,:) = randbreaks();
                case 6 %选择交叉
                    tmp_pop_rte(k,I:J) = fliplr(tmp_pop_rte(k,I:J));
                    tmp_pop_brk(k,:) = randbreaks();
                case 7 % 变异，修改城市点
                    tmp_pop_rte(k,[I J]) = tmp_pop_rte(k,[J I]);
                    tmp_pop_brk(k,:) = randbreaks();
                case 8 %重排序，修改城市点
                    tmp_pop_rte(k,I:J) = tmp_pop_rte(k,[I+1:J I]);
                    tmp_pop_brk(k,:) = randbreaks();
                otherwise 
            end
        end
        new_pop_rte(p-7:p,:) = tmp_pop_rte;
        new_pop_brk(p-7:p,:) = tmp_pop_brk;
    end
    pop_rte = new_pop_rte;
    pop_brk = new_pop_brk;
end
 if show_res
% 画图
    figure('Name','MTSP_GA | Results','Numbertitle','off');
    subplot(2,2,1);
    plot(xy(:,1),xy(:,2),'k.');
    title('City Locations');
    subplot(2,2,2);
    imagesc(dmat(opt_rte,opt_rte));
    title('Distance Matrix');
    subplot(2,2,3);
    rng = [[1 opt_brk+1];[opt_brk n]]';
    for s = 1:salesmen
        rte = opt_rte([rng(s,1):rng(s,2) rng(s,1)]);
        plot(xy(rte,1),xy(rte,2),'.-','Color',clr(s,:));
        title(sprintf('Total Distance = %1.4f',min_dist));
        hold on;
    end
    subplot(2,2,4);
    plot(dist_history,'b','LineWidth',2);
    title('Best Solution History');
    set(gca,'XLim',[0 num_iter+1],'YLim',[0 1.1*max([1 dist_history])]);
end
 % 结果展示
if nargout
    varargout{1} = opt_rte;
    varargout{2} = opt_brk;
    varargout{3} = min_dist;
end
     % 随机城市点的生成
    function breaks = randbreaks()
        if min_tour == 1 
            tmp_brks = randperm(n-1);
            breaks = sort(tmp_brks(1:num_brks));
        else 
            num_adjust = find(rand < cum_prob,1)-1;
            spaces = ceil(num_brks*rand(1,num_adjust));
            adjust = zeros(1,num_brks);
            for kk = 1:num_brks
                adjust(kk) = sum(spaces == kk);
            end
            breaks = min_tour*(1:num_brks) + cumsum(adjust);
        end
    end
end


