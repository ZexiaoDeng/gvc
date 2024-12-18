% ===========================================
% 全局优化引论, 
% R. Horst, P.M. Pardalos, N.V. Thoai 著, 
% 清华大学出版社, 
% P148
% 测试 simplicial_subdivision 子函数
% =========================================


clc ;
clear ;
close all ;

format long

path = './bt-1.3' ;
addpath( path ) ;

w = [ 1.2 ; ...
      0   ; ...
      0.8 ; ] ;     % 单纯形 P 的一个内点

lambda = [ 0.5 ; ...
           0.3 ; ...
           0.0 ; ...
           0.2 ; ] ;       % 有限指标集合
  
% w = [ 1.2 ; ...
%       1.0 ; ...
%       0.8 ; ] ;     % 单纯形 P 的一个内点
%   
% I = [ 1 ; ...
%       1 ; ...
%       1 ; ...
%       1 ; ] ;       % 有限指标集合

% 包含可行集合的单纯形 P
V1 = [ 0 ; 0 ; 0 ; ] ;
V2 = [ 4 ; 0 ; 0 ; ] ;
V3 = [ 0 ; 4 ; 0 ; ] ;
V4 = [ 0 ; 0 ; 4 ; ] ;

P = [ V1, V2, V3, V4 ; ] ;              % 顶点组成的单纯形

Pi = simplicial_subdivision( w, lambda, P ) ;    % 单纯形辐射状剖分

rep.V = P ;
P    = polyh( rep, 'v' ) ;

opt.color = [0.5 0.6 0.1] ;
plot( P, opt ) ;
axis equal
grid on
hold on
plot3( w(1), w(2), w(3), '-rs', 'LineWidth', 2 ) ;



for idx = 1: length( Pi )
    rep.V = Pi{ idx } ;
    P     = polyh( rep, 'v' ) ;
    opt.color = [ 0.5 0.3 0.1*idx ] ;
    plot( P, opt ) ;
end






