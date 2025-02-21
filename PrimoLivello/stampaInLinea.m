function [stringa] = stampaInLinea(x)
% STAMPAINLINEA Ritorna una rappresentazione come stringa dell'array o
% matrice.
%   La rappresentazione di x è stata ideata in modo che sia simile alla
%   rappresentazione di un array in matlab, ovvero con le parentesi
%   quadre attorno, gli spazi per separare gli elementi di una riga e i
%   punto e virgola per separare le righe. Innanzitutto x viene
%   trasformato in un array/matrice di stringhe usando string(). I
%   valori <missing> vengono trasformati nella stringa "NaN". Dopodiché
%   l'array di stringhe viene stampato in accordo alle sue dimensioni.
%   Si noti che sia i vettori riga che i vettori colonna verranno
%   stampati come dei vettori riga, ovvero con gli elementi separati da
%   spazi, però il vettore colonna avrà ' alla fine per indicare che è
%   un vettore colonna.
%   PARAMETRI
%   x: array da stampare in una riga
%   OUTPUT
%   stringa: rappresentazione stringa corrispondente a x
%   ESEMPIO
%   c = [4; 5; 2];
%   A = [0 0.6 0.8;
%       -1 2 0;
%       1 0 -1;
%       -1 0 0;
%       0 -1 0;
%       0 0 -1];
%   b = [500; 0; 0; 0; 0; 0];
%   B = [1 4 5];
%   res = simplessoPrimale(c,A,b,B);
%   stampaStruct(res.steps);
    if isempty(x)
        stringa = "[]";
        return;
    end
    stringa = string(x);
    stringa(ismissing(stringa)) = "NaN";
    if isrow(x)
        stringa = sprintf("[%s]",join(stringa,' '));
    elseif iscolumn(x)
        stringa = sprintf("[%s].'",join(stringa,' '));
    else
        stringa = sprintf("[%s]",join(join(stringa,' '),'; '));
    end
end

