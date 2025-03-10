function dp=poly_dla(t,c)
%% values of derivative of lagrange polynomials through t at c
% function dp=poly_dla(t,c);
% INPUT:
%       t lagrange points in R^m+1
%       c evaluation point(s) in R 
% OUTPUT:
%       dp values of derivative of lagrange polynomials through t at c
%

% (c) DDE-BIFTOOL v. 3.1.1(19), 11/04/2014
%
% 
%
%%
m=length(t)-1;
nc=length(c);
% compute dp:
dp=zeros(m+1,nc);
for j=1:m+1
    for k=1:m+1
        if k~=j
            f=ones(1,nc);
            for l=1:m+1
                if l~=k && l~=j
                    f=f.*(c-t(l))/(t(j)-t(l));
                end;
            end;
            dp(j,:)=dp(j,:)+f/(t(j)-t(k));
        end
    end
end
% reshape to keep compatibility with original version
if nc==1
    dp=dp';
end
end
