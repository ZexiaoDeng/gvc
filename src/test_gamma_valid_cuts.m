% ===========================================
% 全局优化引论, 
% R. Horst, P.M. Pardalos, N.V. Thoai 著, 
% 清华大学出版社, 
% P152
% 测试 gamma_valid_cuts 子函数
% ===========================================

clc ;
clear ;
close all ;

path = './bt-1.3' ;
addpath( path ) ;

Aineq = [ -.5,  1 ; ...
            1,  1 ; ...
            1, -1 ; ...
           -1,  0 ; ...     % lb
            0, -1 ; ] ;
bineq = [ 1 ; ...
          4 ; ...
          2 ; ...
          0 ; ...
          0 ; ] ;
      
oracle = @(x) -( x(1) - 1 )^2 - ( x(2) - 1 )^2 + 2 ;

rep.B = Aineq ;
rep.b = bineq ;

P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep
CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
Adj = adj( P ) ;                  % 获取顶点对应的链表( 邻接表表示形式 )

opt.color = [ 0, 1, 1 ] ;
plot( P, opt ) ;
axis equal ;
grid on ;
hold on ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', ...
        CH.V( 1, i ), ...
        CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', Adj{i}(1), Adj{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

% % 找出当前最好的可行点
% fval  = zeros( length( Adj ), 1 ) ;
% for i = 1: length( Adj )
%      fval( i, 1 ) = feval( oracle, CH.V( :, i ) ) ;
% end
% [ gamma, idxopt ] = min( fval )
% v                 = CH.V( :, idxopt ) ;
% [ lhs, rhs ] = gamma_valid_cuts( gamma, v, CH.V, Adj{idxopt}, oracle )

% 找出当前最好的可行点
gamma = feval( oracle, CH.V( :, 1 ) ) ;
v     = CH.V( :, 1 ) ;
[ lhs, rhs ] = gamma_valid_cuts( gamma, v, CH.V, Adj{1}, oracle )

rep.B = [ rep.B ; -lhs ; ] ;
rep.b = [ rep.b ; -rhs ; ] ;

P   = eval( polyh( rep, 'h' ) ) ; % 求出多胞体 P 的 H-rep, V-rep, P-rep
CH  = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
Adj = adj( P ) ;                  % 获取顶点对应的链表( 邻接表表示形式 )

opt.color = [ 1, 0, 0 ] ;
plot( P, opt ) ;
axis equal

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', ...
        CH.V( 1, i ), ...
        CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', Adj{i}(1), Adj{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

