function [ lhs, rhs ] = gamma_valid_cuts( gamma, v, Vi, Adj, oracle )
%   	GAMMA_VALID_CUTS  计算凹性割
%
%       H_+ = { x: e'*Y^(-1 )*x >= 1 + e'*Y^(-1 )*v }
%
%   输入:
%       gamma  : 当前最好目标值
%       v      : 单纯形 P 非退化顶点
%       Vi     : v 的邻接顶点集合
%       Adj    : v 邻接顶点的指标集合 
%       oracle : 目标函数
%
%   输出:
%       lhs    : 凹性割左边
%       rhs    : 凹性割右边
%
%               lhs*x >= rhs
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P151
%

rows    = length( v   ) ;      % 决策变量个数 
columns = length( Adj ) ;      % 邻接顶点个数

Y     = zeros( rows, columns ) ;
for i = 1: columns
    d = Vi( :, Adj(i) ) - v ;     % 给出射线方向

    % 计算 gamma-扩张
    [ y, theta ] = gamma_extension( gamma, v, d, oracle ) ;
    Y( :, i ) = y - v ;
%   Y( :, i ) = theta*d ;
end

% 缺点: 还得主元归一化
e   = ones( rows, 1 ) ;
lhs = e'*Y^( -1 ) ;           % H_plus 正半空间
rhs = 1 + lhs*v ;









end



