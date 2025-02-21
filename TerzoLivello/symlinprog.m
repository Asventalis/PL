function [x,fval] = symlinprog(c,Aleq,bleq,Aeq,beq,l,u)
    %SYMLINPROG Risolve un problema di programmazione lineare generale
    %usando il metodo del simplesso primale
    %   Risolve il problema di minimizzazione
    %   min c'*x
    %   s.t. Aleq*x <= bleq
    %        Aeq*x = beq
    %        l <= x <= u
    %   con il metodo del simplesso
    %   INPUT
    %   c: vettore dei coefficienti della funzione obiettivo
    %   Aleq: matrice dei coefficienti delle disequazioni
    %   bleq: vettore dei termini noti delle disequazioni
    %   Aeq: matrice dei coefficienti delle equazioni
    %   beq: vettore dei termini noti delle equazioni
    %   l: vettore dei limiti inferiori
    %   u: vettore dei limiti superiori
    %   OUTPUT
    %   res: struct ottenuta con la funzione linprog
    
    % Controllo parametri
    if nargin < 3
        error('Inserire almeno una funzione obiettivo e dei vincoli');
    elseif nargin == 3
        Aeq = [];
        beq = [];
        l = [];
        u = [];
    elseif nargin == 4
        error('Non si è inserito il vettore dei termini noti per le equazioni');
    elseif nargin == 5
        l = [];
        u = [];
    elseif nargin == 6
        u = [];
    elseif nargin > 7
        error('Troppi parametri inseriti');
    end
    
    % Solo uno tra Aleq e Aeq può essere vuoto
    if isempty(Aleq) && isempty(Aeq)
        error('Inserire almeno un vincolo');
    end
    
    % Convertire i parametri vuoti in parametri di default
    if isempty(Aleq)
        Aleq = zeros(0,length(c));
        bleq = zeros(0,1);
    end
    
    if isempty(Aeq)
        Aeq = zeros(0,length(c));
        beq = zeros(0,1);
    end
    
    if isempty(l)
        l = -inf*ones(length(c),0);
    end
    
    if isempty(u)
        u = Inf*ones(length(c),0);
    end
    
    % Ruota vettori riga in vettori colonna
    if size(c,2) > 1
        c = c';
    end

    if size(bleq,2) > 1
        bleq = bleq';
    end

    if size(beq,2) > 1
        beq = beq';
    end

    if size(l,2) > 1
        l = l';
    end

    if size(u,2) > 1
        u = u';
    end

    % Controllo dimensioni
    if length(c) ~= max(size(l))
        error('La lunghezza di l deve essere uguale alla lunghezza di c');
    end
    
    if length(c) ~= max(size(u))
        error('La lunghezza di u deve essere uguale alla lunghezza di c');
    end
    
    if length(c) ~= size(Aleq,2)
        error('Il numero di colonne di Aleq deve essere uguale alla lunghezza di c');
    end
    
    if length(c) ~= size(Aeq,2)
        error('Il numero di colonne di Aeq deve essere uguale alla lunghezza di c');
    end
    
    if length(bleq) ~= size(Aleq,1)
        error('Il numero di righe di Aleq deve essere uguale alla lunghezza di bleq');
    end
    
    if length(beq) ~= size(Aeq,1)
        error('Il numero di righe di Aeq deve essere uguale alla lunghezza di beq');
    end
    
    % Controllo del tipo di variabili
    controlloNumSym(c);
    controlloNumSym(Aleq);
    controlloNumSym(bleq);
    controlloNumSym(Aeq);
    controlloNumSym(beq);
    controlloNumSym(l);
    controlloNumSym(u);
    
    % Converte i limiti in vincoli
    A_l = -eye(length(l));
    A_u = eye(length(u));
    b_l = -l;
    b_u = u;
    
    % Converte il problema in forma Primale Standard
    c_primale = -c;
    A_primale = [Aleq; Aeq; -Aeq; A_l; A_u];
    b_primale = [bleq; beq; -beq; b_l; b_u];
    
    controlloPrimale(c_primale,A_primale,b_primale);
    
    % Risolve il problema in forma Primale Standard
    res = simplessoDuale(c_primale,A_primale,b_primale);
    x = res.x;
    fval = -res.val;
    end
    
    