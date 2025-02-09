function [baseAmm,info] = ausiliarioDuale(A,c)
    %AUSILIARIODUALE Trova una base ammissibile nel duale per un problema
    %di programmazione lineare in forma primale standard.
    %   Trova una base ammissibile costruendo il problema duale ausiliario,
    %   e procedendo a cercare la soluzione. Se la soluzione ottimale ha 
    %   valore diverso da 0, ritorna NaN al posto di baseAmm, altrimenti
    %   ricava la base ammissibile dalla soluzione ottimale del problema
    %   ausiliario. Ritorna inoltre info varie sullo svolgimento
    %   dell'algoritmo.
    %   PARAMETRI
    %   A: matrice dei vincoli
    %   c: funzione obiettivo
    %   OUTPUT
    %   baseAmm: una base ammissibile per il problema duale, oppure NaN
    %   info: struct:
    %     - cAux: funzione obiettivo del problema ausiliario
    %     - AAux: matrice dei vincoli del problema ausiliario
    %     - bAux: matrice dei termini noti del problema ausiliario
    %     - baseAux: base ammissibile da cui partire per risolvere il
    %       problema ausiliario con il simplesso duale
    %     - resSimplesso: output di simplessoDuale applicato al problema
    %       duale ausiliario
    %   ESEMPIO
    %   c = [4; 5; 2];
    %   A = [0 0.6 0.8;
    %       -1 2 0;
    %       1 0 -1;
    %       -1 0 0;
    %       0 -1 0;
    %       0 0 -1];
    %   b = [500; 0; 0; 0; 0; 0];
    %   B = [1 2 6];
    %   [baseAmm,info] = ausiliarioDuale(A,b);
    %   stampaStruct(info.resSimplesso.steps);
    %   baseAmm
    
    if ~exist("c", "var") || ~exist("A", "var")
        error("Errore. c e A devono essere tutti e due specificati.");
    end
    controlloAc(A,c);
    A = sym(A);
    c = sym(c);
    [nVincoli, nVariabili] = size(A);
    % costruisco il problema ausiliario duale min bAux'*y con vincoli y'*AAux=cAux' e y>=0
    cAux = c.*sign(c);
    AAux = [A.*sign(c)'; eye(nVariabili)];
    bAux = sym([zeros(nVincoli,1); ones(nVariabili,1)]);
    baseAux = nVincoli+1:nVincoli+nVariabili;
    res = simplessoDuale(cAux, AAux, bAux, baseAux); % risolvo il problema ausiliario duale
    if res.val ~= 0 % se il valore ottimo non è 0, non esiste base ammissibile per il duale
        baseAmm = NaN;
    else % se il valore ottimo è 0, si può trovare una base ammissibile a partire dalla base ottima del duale ausiliario
        baseAmm = res.base;
        % algoritmo per trovare una base ammissibile a partire da una base ottima dell'ausiliario (ipotizziamo che la dimostrazione nel libro del prof sia giusta)
        while any(baseAmm>nVincoli)
            yB = baseAmm(baseAmm<=nVincoli);
            yN = setdiff(1:nVincoli,yB);
            epsilonB = baseAmm(baseAmm>nVincoli);
            M = AAux(baseAmm,:)^-1;
            H = A(yN,:)*M(:,baseAmm>nVincoli);
            [row,col] = find(H,1);
            k = yN(row);
            h = epsilonB(col);
            baseAmm(baseAmm==h) = k;
        end
    end
    % salvo informazioni importanti
    info.cAux = cAux;
    info.AAux = AAux;
    info.bAux = bAux;
    info.baseAux = baseAux;
    info.resSimplesso = res;
end

