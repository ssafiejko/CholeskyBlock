function [L] = basiccholesky(A)

%funkcja zwraca macierz L, która jest dekompozycją macierzy A
%spełniająca warunek postaci A = LL^T 

%A - macierz symetryczna, dodatnio określona
%L - macierz dolnotrójkątna

tic
if ~isequal(A', A)
    ME = MException("Cholesky:wrongInput", "Macierz A nie jest symetryczna");
        throw(ME)
end

if size(A,1) ~= size(A,2)
    ME = MException("Cholesky:wrongInput", "Złe wymiary macierzy");
        throw(ME)
end

n = size(A,1);
L = zeros(n,n);

for k = 1:n
    L(k,k) = sqrt(A(k,k) - sum(L(k,1:k-1).^2));

    L(k+1:n,k) = (A(k+1:n,k) - (sum(L(k+1:n,1:k-1).*L(k,1:k-1),2)))/L(k,k);
    
end
end