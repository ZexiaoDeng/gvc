function [ x, gamma ] = gvcsolve( oracle, Aineq, bineq, opts )
%   GVCSOLVE finds a constrained minimum of a function of several variables.
%   GVCSOLVE attempts to solve concave programming problems of the form:
%
%       min    f(x) = x'*Q*x + 2*p'*x 
%       s.t.   Aineq*x <= bineq        ( linear constraints )
%    
%   GVCSOLVE implements four different algorithms: 
%   (1) gamma-valid cut method ( gvc ) ( concave cut )
%
%   Input:
%       Q    : negative semidefinite matrix of the concave function f(x)
%       p    : column vector of the concave function f(x)
%       Aineq: system matrix of the linear inequalities
%       bineq: column vector of the linear inequalities
%   Output:
%      gamma: optimal value
%      x    : an optimal solution (column vector)
%
%    see also 全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社

maxiter = opts.maxiter ;            % 最大迭代步数
epsilon = opts.epsilon ;            % 精度控制
% ==============================
% 初始化
% ==============================
row   = size( Aineq, 2 ) ;
n     = row ;
rep.M = eye( n ) ;
rep.B = Aineq ;
rep.b = bineq ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P  = eval( polyh( rep ) ) ;      % 求出多胞体 P 的 H-rep, V-rep, P-rep
CH    = vrep( P ) ;              % 获取多胞体 P 的 V-rep
Adj     = adj( P ) ;             % 获取顶点对应的链表( 邻接表表示形式 )

flag    = 1 ;                       % 循环结束标志
k       = 1 ;                       % 迭代步数计数器

while k < maxiter
    
    % 找出当前最好的可行点
    fval = zeros( n, 1 ) ;
    for i = 1: length( Adj )
         fval( i, 1 ) = feval( oracle, CH.V( :, i ) ) ;
    end
    [ gamma, idxopt ] = min( fval ) ;
    x     = CH.V( :, idxopt ) ;
    idx   = randi( n ) ;                % 随机选取其中一个顶点
    while idx == idxopt
        idx = randi( n ) ;
    end
    v     = CH.V( :, idx ) ;
    
    column = length( Adj{idx} ) ;
    Y   = zeros( row, column ) ;
    for i = 1: column
        d = CH.V( :, Adj{idx}(i) ) - v ;     % 给出射线方向
        
        % 计算 gamma-扩张
        [ y, theta ] = gamma_extension( gamma, v, d, oracle ) ;
        Y( :, i ) = y - v ;
%         Y( :, i ) = theta*d ;
    end
    % 缺点: 还得主元归一化
    e     = ones( row, 1 ) ;
    H_vec = e'*Y^( -1 ) ;           % H_plus 正半空间
    H_sca = 1 + H_vec*v ;
    
    % 停止准则 P <= H_minus, 剩余多胞体为 H_minus 负半空间的子集
    rep1.M = eye( row ) ;
    rep1.B = H_vec ;
    rep1.b = H_sca ;
    H_minus = polyh( rep1 ) ;
    
    flag = ~le( P, H_minus, epsilon ) ;      % 判别 P <= H_minus
    if flag == 0
        break ;
    end
    
    % 更新多胞体信息
    rep.M = eye( row ) ;
    rep.B = [ rep.B ; -H_vec ; ] ;
    rep.b = [ rep.b ; -H_sca ; ] ;
    rep.l = zeros( 2, 1 ) ;
    rep.u = [ inf ; inf ; ] ;

    P = eval( polyh( rep ) ) ;
    CH   = vrep( P ) ;              % 获取多胞体 P 的 V-rep
    Adj    = adj( P ) ;               % 获取顶点对应的链表( 邻接表表示形式 )
    n    = size( CH.V, 2 ) ;
    
    k = k + 1 ;
end

k
    
    
    
end











