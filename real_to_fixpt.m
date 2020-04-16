function y = real_to_fixpt(x,a,b)
%a=1;
% if nargin < 1
%  a = 0;
% end
% if isempty(a)
%  a = 0;
% end
%b=15;
% if nargin < 1
%  b = 16;
% end
% if isempty(b)
%  b = 16;
% end
% total bits
N = a + b;
y = [];
%S = 0;
for i = 1:length(x)
    if x(i)>=0
     d = dec2bin((x(i)*2^b),N);
    else
     d = dec2bin(2^N-ceil(-x(i)*2^b),N);
    end
y(i,:) = d;
end
y= char(y);
