function [X,t] = choleskysolve(A,B)
tic
% funkcja rozwiązuje równanie
% AX = B korzystając z metody blokowej Cholesky'ego
% przyjmuje dwa argumenty:
% macierz symetryczną A dodatnio określoną o n wierszach,
% macierz B o m kolumnach i n wierszach
% funkcja zwraca macierz X o m kolumnach i n wierszach będącą rozwiązaniem
% równania AX = B


L = CholeskyBlock(A);

% LY = B
Y = solve_upperorlower(L, B, "lower");

% L'X = Y
X = solve_upperorlower(L', Y, "upper");

% X spełnia równanie L(L'X) = LY = B
t=toc;
end