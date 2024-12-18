function Pi = conical_subdivision( w, lambda, P )
%   CONICAL_SUBDIVISION  计算锥形剖分
%
%       P = Union Pi(w) i in I
%       w in S
%       w = sum( lambda_i*di ), lambda_i >= 0, i = 1, ..., n
%       sum( lambda_i ) = 1, i = 1, ..., n
%
%   输入:
%       w      : 单纯形 S 的内方向
%       lambda : 一个有限指标集, int( Pi )非空
%       P      : 原多面锥, 
%
%   输出:
%       Pi     : 子多面锥 cell
%                Pi = [ d1, d2, ..., di-1, w, di+1, ..., dn ]
%
%    see also 
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P148
%

Pi    = {} ;
icout = 1  ;

for idx = 1: size( P, 2 )
    if lambda( idx ) > 0
        Pi{ icout } = P ;
        Pi{ icout }( : , idx ) = w ;
        icout = icout + 1 ;
    else
        continue ;
    end
end

return ;










end



