function [blad_wzgledny] = check_result(A,B)
%funkcja sprawdza poprawność przykładów

X = choleskysolve(A, B);
z = inv(A)*B;
L=CholeskyBlock(A);
wskaznik_uwar=cond(A);
blad_rozkladu= norm(A-L*L')/norm(A);
blad_wzgledny= norm(z-X)/norm(z);
wsp_stabilnosci= norm(X-z)/(norm(z)*cond(A));
wsp_poprawnosci= norm(B-A*X)/(norm(A)*norm(X));


tabelka=table(wskaznik_uwar, blad_rozkladu, blad_wzgledny, wsp_stabilnosci, wsp_poprawnosci, 'VariableNames', ["Wskaźnik uwarunkowania", "Błąd rozkładu", "Błąd względny", "Współczynnik stabilności", "Współczynnik poprawności"]);
end