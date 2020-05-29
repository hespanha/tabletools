function W=myPCA(X,d)

% centering

if 0
    mu=ones(size(X,1),1)*(sum(X,1)/size(X,1));
    X=X-mu;
    CV=X*X';
else
    mu=sum(X,1)/size(X,1);
    CV=X*X'-mu*mu';
end

if issparse(CV)
    [U,S,V]=svds(CV,d);
    W=U*S;
else
    [U,S,V]=svd(CV);
    W=U*S(1:d,1:d);
end

