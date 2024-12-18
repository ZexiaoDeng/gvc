% ===========================================
% 全局优化引论, 
% R. Horst, P.M. Pardalos, N.V. Thoai 著, 
% 清华大学出版社, 
% P161
% 测试 gvc_solver 子函数
% ===========================================

clc ;
clear ;
close all ;

path = './bt-1.3' ;
addpath( path ) ;

format long ;

Aineq = [ -4,  2 ; 
           0,  1 ; 
           1,  1 ;
           1,  0 ;
           1, -4 ; ] ;
bineq = [ 1 ; 
          2 ; 
          4 ; 
          3 ; 
          1 ;] ;
opts.maxiter = 15 ;            % 最大迭代步数
opts.epsilon = 1e-9 ;          % 精度控制

[ x, gamma ] = gvc_solver( @oracle, Aineq, bineq, opts )

function f = oracle( x )
    Q = [ -1,  2 ;
           2, -4 ; ] ;
    p = [ 1 ; 2 ; ] ;
    f = x'*Q*x + 2*p'*x ;
end

  





