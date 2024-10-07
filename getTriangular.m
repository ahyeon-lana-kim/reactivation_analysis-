function [vec]=getTriangilar(mat,Lower)
if Lower==1
    mat=tril(mat);
else
    mat=triu(mat);
end
mat(mat==0)=nan;
vec=mat(~isnan(mat));
end
