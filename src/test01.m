

clc ;
clear ;
close all ;

path = './bt-1.3' ;
addpath( path ) ;

global v1 d1 v2 d2 v3 d3 v4 d4
global gamma1 gamma2 gamma3 gamma4 Q p

% ========================
% min    f(x) = x'*Q*x + 2*p'*x 
% s.t.   A*x <= b
% 
Q = [ -1,  2 ;
       2, -4 ; ] ;
p = [ 1 ; 2 ; ] ;

A = [ -4,  2 ; 
       0,  1 ; 
       1,  1 ;
       1,  0 ;
       1, -4 ; ] ;
b = [ 1 ; 
      2 ; 
      4 ; 
      3 ; 
      1 ;] ;

rep.M = eye( 2 ) ;
rep.B = A ;
rep.b = b ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P = polyh( rep ) ;

plot( P ) ;
grid on ;
hold on ;


P=eval(P);
CH = vrep(P) ;
A = adj( P ) ;
for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', CH.V( 1, i ), CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', A{i}(1), A{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

% 初始化
flag = 1 ;
k = 1 ;

v1 = CH.V( :, 5 ) ;
f1  = v1'*Q*v1 + 2*p'*v1
gamma1 = f1 ;


d1 = CH.V( :, A{5}(1) ) - v1
d2 = CH.V( :, A{5}(2) ) - v1


fun = @(t) -t ;
t0 = 0 ;
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[ t1, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon1, options )
[ t2, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon2, options )

theta1 = t1 ;
theta2 = t2 ;

y1 = v1 + theta1*d1
y2 = v1 + theta2*d2

% 凹性割
e = [ 1 ; 1 ; ] ;
Y = [ y1 - v1, y2 - v1 ] ;
H_vec = e'*Y^( -1 )
H_sca = 1 + e'*Y^( -1 )*v1

ezplot( '0.5*x1 + 1*x2 - 1 = 0', [-1 4 -1 3] ) ;

rep.M = eye( 2 ) ;
rep.B = [ rep.B ; -H_vec ; ] ;
rep.b = [ rep.b ; -H_sca ; ] ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P = polyh( rep ) ;

opt.color = [0.5 0.5 0.5] ;
plot( P, opt ) ;

P=eval(P);
CH = vrep(P) ;
A = adj( P ) ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', CH.V( 1, i ), CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', A{i}(1), A{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

% ===============================
% 第二次迭代
% ==============================
flag = 1 ;
k = 2 ;
v2 = CH.V( :, 5 ) ;

d1 = CH.V( :, A{5}(1) ) - v2
d2 = CH.V( :, A{5}(2) ) - v2

% 找出最好的可行点
for i = 1: length( A )
     fval( i, 1 ) = CH.V( :, i )'*Q*CH.V( :, i ) + 2*p'*CH.V( :, i ) ;
end
[ f2, idx ] = min( fval )
x2 = CH.V( :, idx )
gamma2 = f2 ;

fun = @(t) -t ;
t0 = 0 ;
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[ t1, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon3, options )
[ t2, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon4, options )

theta1 = t1
theta2 = t2

y1 = v2 + theta1*d1
y2 = v2 + theta2*d2

% 凹性割
e = [ 1 ; 1 ; ] ;
Y = [ y1 - v2, y2 - v2 ] ;
H_vec = e'*Y^( -1 )
H_sca = 1 + e'*Y^( -1 )*v2

ezplot( '-0.800000016000427*x1 - 1.387355764807656*x2 + 3.374711561616165 = 0', [-1 4 -1 3] ) ;

rep.M = eye( 2 ) ;
rep.B = [ rep.B ; -H_vec ; ] ;
rep.b = [ rep.b ; -H_sca ; ] ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P = polyh( rep ) ;

opt.color = [0.2 0.2 1] ;
plot( P, opt ) ;

P=eval(P);
CH = vrep(P) ;
A = adj( P ) ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', CH.V( 1, i ), CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', A{i}(1), A{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end


% ===============================
% 第 3 次迭代
% ==============================
k = 3 ;
v3 = CH.V( :, 3 ) ;

d1 = CH.V( :, A{3}(1) ) - v3
d2 = CH.V( :, A{3}(2) ) - v3

% 找出最好的可行点
for i = 1: length( A )
     fval( i, 1 ) = CH.V( :, i )'*Q*CH.V( :, i ) + 2*p'*CH.V( :, i ) ;
end
[ f3, idx ] = min( fval )
x3 = CH.V( :, idx )
gamma3 = f3 ;

fun = @(t) -t ;
t0 = 0 ;
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[ t1, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon5, options )
[ t2, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon6, options )

theta1 = t1
theta2 = t2

y1 = v3 + theta1*d1
y2 = v3 + theta2*d2

% 凹性割
e = [ 1 ; 1 ; ] ;
Y = [ y1 - v3, y2 - v3 ] ;
H_vec = e'*Y^( -1 )
H_sca = 1 + e'*Y^( -1 )*v3

ezplot( '0.779220816572490*x1 + 0.519480575900015*x2 - 1.623376674351069 = 0', [-1 4 -1 3] ) ;

rep.M = eye( 2 ) ;
rep.B = [ rep.B ; -H_vec ; ] ;
rep.b = [ rep.b ; -H_sca ; ] ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P = polyh( rep ) ;

opt.color = [0 0.2 1] ;
plot( P, opt ) ;

P=eval(P);
CH = vrep(P) ;
A = adj( P ) ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f', CH.V( 1, i ), CH.V( 2, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', A{i}(1), A{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

% ===============================
% 第 4 次迭代
% ==============================
flag = 1 ;
k    = 4 ;

% 找出最好的可行点
for i = 1: length( A )
     fval( i, 1 ) = CH.V( :, i )'*Q*CH.V( :, i ) + 2*p'*CH.V( :, i ) ;
end
[ gamma4, idx ] = min( fval )

v4 = CH.V( :, 2 ) ;
d1 = CH.V( :, A{2}(1) ) - v4
d2 = CH.V( :, A{2}(2) ) - v4

fun = @(t) -t ;
t0 = 0 ;
options = optimoptions('fmincon','Display','off','Algorithm','interior-point');
[ theta1, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon7, options )
[ theta2, fval ] = fmincon( fun, t0, [], [], [], [], [], [], @nonlcon8, options )

y1 = v4 + theta1*d1
y2 = v4 + theta2*d2

% 凹性割
e = [ 1 ; 1 ; ] ;
Y = [ y1 - v4, y2 - v4 ] ;
H_vec = e'*Y^( -1 )
H_sca = 1 + e'*Y^( -1 )*v4

ezplot( '-0.014013761141702*x1 + 0.556314097605284*x2 - 1.102117800043021 = 0', [-1 4 -1 3] ) ;

rep.M = eye( 2 ) ;
rep.B = [ rep.B ; -H_vec ; ] ;
rep.b = [ rep.b ; -H_sca ; ] ;
rep.l = zeros( 2, 1 ) ;
rep.u = [ inf ; inf ; ] ;

P = polyh( rep ) ;


function [ c, ceq ] = nonlcon1( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v1 + t*d1 )'*Q*( v1 + t*d1 ) + 2*p'*( v1 + t*d1 ) - gamma1 ) ;
    ceq = [] ;
end
function [ c, ceq ] = nonlcon2( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v1 + t*d2 )'*Q*( v1 + t*d2 ) + 2*p'*( v1 + t*d2 ) - gamma1 ) ;
    ceq = [] ;
end

function [ c, ceq ] = nonlcon3( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v2 + t*d1 )'*Q*( v2 + t*d1 ) + 2*p'*( v2 + t*d1 ) - gamma2 ) ;
    ceq = [] ;
end
function [ c, ceq ] = nonlcon4( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v2 + t*d2 )'*Q*( v2 + t*d2 ) + 2*p'*( v2 + t*d2 ) - gamma2 ) ;
    ceq = [] ;
end

function [ c, ceq ] = nonlcon5( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v3 + t*d1 )'*Q*( v3 + t*d1 ) + 2*p'*( v3 + t*d1 ) - gamma3 ) ;
    ceq = [] ;
end
function [ c, ceq ] = nonlcon6( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v3 + t*d2 )'*Q*( v3 + t*d2 ) + 2*p'*( v3 + t*d2 ) - gamma3 ) ;
    ceq = [] ;
end

function [ c, ceq ] = nonlcon7( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v4 + t*d1 )'*Q*( v4 + t*d1 ) + 2*p'*( v4 + t*d1 ) - gamma4 ) ;
    ceq = [] ;
end
function [ c, ceq ] = nonlcon8( t )
    global v1 d1 v2 d2 v3 d3 v4 d4
    global gamma1 gamma2 gamma3 gamma4 Q p
    c = -( ( v4 + t*d2 )'*Q*( v4 + t*d2 ) + 2*p'*( v4 + t*d2 ) - gamma4 ) ;
    ceq = [] ;
end








