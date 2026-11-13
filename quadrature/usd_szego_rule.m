function [nodes,weights,info]=usd_szego_rule(U,v,m)
%USD_SZEGO_RULE m-node unit-circle quadrature from an isometric Arnoldi run.
%   Implements fully reorthogonalized Arnoldi on v,Uv,... and completes the
%   first m-1 unitary Hessenberg columns with a unimodular last column.
%   The resulting vector spectral measure matches v'*U^k*v for
%   |k|<=m-1, subject to rounding. No eigenvalues of U are computed.
sz=U([],'size');n=sz(1);
assert(numel(sz)==2 && sz(2)==n,'U must be square.');
validateattributes(m,{'numeric'},{'scalar','integer','positive','<=',n});
v=v(:);
assert(numel(v)==n && norm(v)>0,'Starting vector must be nonzero.');
v=v/norm(v);
Q=zeros(n,m);Q(:,1)=v;
H=zeros(m,m);
calls=0;used=m;closed=false;
for j=1:m-1
    w=U(Q(:,j),'notransp');calls=calls+1;
    wnorm=norm(w);
    for pass=1:2
        for i=1:j
            a=Q(:,i)'*w;
            H(i,j)=H(i,j)+a;
            w=w-Q(:,i)*a;
        end
    end
    b=norm(w);
    if b<=1e-12*max(1,wnorm)
        % Krylov subspace is invariant: projected matrix is already unitary.
        used=j;H=H(1:j,1:j);closed=true;break
    end
    H(j+1,j)=b;
    Q(:,j+1)=w/b;
end
if ~closed
    h=Q(:,1:m)'*U(Q(:,m),'notransp');calls=calls+1;
    if m==1
        last=1;
    else
        last=null(H(:,1:m-1)');
        assert(size(last,2)==1,'Arnoldi columns lost independence.');
    end
    phase=last'*h;
    if abs(phase)>1e-14
        last=last*(phase/abs(phase));
    end
    H(:,m)=last;
end
[V,D]=eig(H);
nodes=angle(diag(D));
weights=real(abs(V(1,:)).^2).';
weights=weights/sum(weights);
info=struct('nodes',used,'matvecs',calls, ...
    'unitarity_residual',norm(H'*H-eye(used),'fro'));
end
