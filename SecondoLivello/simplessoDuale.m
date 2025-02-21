function [res] = simplessoDuale(c,A,b,B,maxIter)
% SIMPLESSODUALE Risolve un problema di programmazione lineare in forma
% primale standard.
%   Segue l'algoritmo del simplesso duale per risolvere il problema
%   di massimizzazione di c'*x con vincoli A*x<=b. Usa B come base di
%   partenza. B deve essere ammissibile per il duale. Se B non viene
%   passata, genera una base ammissibile con ausiliarioDuale, se
%   esiste. Ritorna res.val=NaN se non esiste soluzione, res.val=-inf
%   se la soluzione è infinita, e ritorna un valore finito se la
%   soluzione esiste. Le altre info sono settate in accordo. Ritorna
%   inoltre steps, ovvero una struct contenente informazioni per ogni
%   passo dell'algoritmo che è stato svolto. Steps può essere stampato
%   in modo elegante sulla console utilizzando stampaStruct.
%   PARAMETRI
%   c: vettore colonna della funzione obiettivo
%   A: matrice dei vincoli
%   b: vettore colonna dei termini noti
%   B: base ammissibile da cui partire
%   maxIter: numero massimo di iterazioni eseguibili
%   OUTPUT
%   res: struct:
%     - x: soluzione primale di base ottima
%     - y: soluzione duale di base ottima
%     - base: base ottima
%     - val: valore della soluzione di base ottima
%     - steps: struct:
%         - base: base dello step i-esimo
%         - invAB: l'inversa di A(base,:)
%         - x: soluzione primale della base
%         - y: soluzione duale della base
%         - k: indice entrante
%         - Ak: colonna di indice h della matrice W=-A(base,:)^-1
%         - Ak_times_W: prodotto tra Ak e W
%         - rapporti: i rapporti usati per calcolare h
%         - h: indice entrante
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
    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var") 
        error("Errore. c, A e b devono essere tutti e 3 specificati.");
    end
    c = sym(c);
    A = sym(A);
    b = sym(b);
    controlloPrimale(c, A, b);
    [nVincoli,nVariabili] = size(A);
    if ~exist("B", "var")
        B = ausiliarioDuale(A,c);
        if isnan(B) % se ausiliarioDuale torna NaN il problema è vuoto
            res.x = NaN;
            res.y = NaN;
            res.base = NaN;
            res.val = NaN;
            res.steps = NaN;
            return;
        end
    end
    controlloBase(B, A);
    if ~exist("maxIter","var")
        maxIter = -1;
    end
    if maxIter ~= -1 && (maxIter<0 || round(maxIter) ~= maxIter)
        error("Errore. MaxIter deve essere un numero naturale oppure -1.");
    end
    iter = 0;
    yB = (c'*A(B,:)^-1)';
    if any(yB<0)
        error("Errore. Base [%s] non ammissibile per il duale.", join(string(B),' '));
    end
    steps = [];
    step = struct;
    step.base = NaN;
    step.invAB = NaN;
    step.x = NaN;
    step.y = NaN;
    step.k = NaN;
    step.Ak = NaN;
    step.Ak_times_W = NaN;
    step.rapporti = NaN;
    step.h = NaN;
    while maxIter == -1 || iter ~= maxIter
        iter = iter + 1;
        step.base = B;
        N = setdiff(1:nVincoli, B);
        step.invAB = A(B,:)^-1;
        W = -A(B,:)^-1;
        x = - W * b(B,:); % calcola x
        step.x = x;
        yB = - c' * W; % calcola yB
        y = sym(zeros(nVincoli, 1));
        y(B,:) = yB;
        step.y = y;
        if all(A(N,:)*x-b(N) <= 0) % se i vincoli non di base sono tutti rispettati, siamo all'ottimo
            steps = [steps; step];
            res.x = x;
            res.y = y;
            res.base = B;
            res.val = c' * x;
            res.steps = steps;
            return;
        end
        k = find(A*x-b>0,1); % altrimenti calcola k
        step.k = k;
        Ak = A(k,:);
        step.Ak = Ak;
        step.Ak_times_W = Ak*W;
        if all(Ak*W >= 0) % se tutti gli Ak*Wi con i in B sono positivi, la soluzione è -infinito
            steps = [steps; step];
            eyeN = eye(nVincoli);
            eyeN = eyeN(N,:);
            M = [A';eyeN]^-1;
            Mk = M(:,nVariabili+find(N == k));
            Mk(Mk~=0) = Mk(Mk~=0)*inf;
            res.x = NaN;
            res.y = y + Mk;
            res.base = NaN;
            res.val = -inf;
            res.steps = steps;
            return;
        end
        BnegIndex = find(Ak*W < 0);
        rapporti = (yB(BnegIndex)./(-Ak*W(:,BnegIndex))); % altrimenti calcola i rapporti
        step.rapporti = rapporti;
        rapportoMinimo = min(rapporti);
        h = B(BnegIndex(find(rapporti == rapportoMinimo,1))); % calcola h
        step.h = h;
        steps = [steps; step];
        step = struct;
        step.base = NaN;
        step.invAB = NaN;
        step.x = NaN;
        step.y = NaN;
        step.k = NaN;
        step.Ak = NaN;
        step.Ak_times_W = NaN;
        step.rapporti = NaN;
        step.h = NaN;
        B(B==h) = k; % metti k al posto di h
        B = sort(B);
    end
    res.x = x;
    res.y = y;
    res.base = B;
    res.val = c' * x;
    res.steps = steps;
end

