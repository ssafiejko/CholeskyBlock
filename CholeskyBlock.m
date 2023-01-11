function [L,t] = CholeskyBlock(A)
tic

if size(A,1)~=size(A,2) | rem(size(A,1),2)~=0 | ~isequal(A', A) | ~isequal(A(1:size(A)/2, (floor(size(A,1)/2)+1):size(A)), eye(size(A,1)/2))
    ME = MException("Cholesky:wrongInput", "Macierz A nie jest w postaci z polecenia");
        throw(ME)
end

if ~all(eig(A)>0)
    ME = MException("Cholesky:wrongInput", "Macierz A nie jest dodatnio okreslona");
        throw(ME)
end

L=zeros(size(A,1));
n=size(A,1)/2;

    A11= A(1:size(A,1)/2,1:size(A,1)/2);
    A22= A((floor(size(A,1)/2)+1):size(A,1), (floor(size(A,1)/2)+1):size(A,1));

    I=A(1:size(A)/2, (floor(size(A,1)/2)+1):size(A));

    A11= (basiccholesky(A11));
    A21= I*inv(A11');
    A22= basiccholesky(A22- A21*(A21'));
    L=[A11,zeros(n);A21,A22];
   
 t=toc; 
end