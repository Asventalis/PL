function [x,fval] = symintlinprog(c,Aleq,bleq,Aeq,beq,l,u)
    %SYMINTLINPROG Risolve un problema di programmazione lineare intera in forma generale
    %usando il metodo dei tagli di gomory
    %   Risolve il problema di minimizzazione
    %   min c'*x
    %   s.t. Aleq*x <= bleq
    %        Aeq*x = beq
    %        l <= x <= u
    %        x interi
    %   con il metodo dei tagli di gomory
    %   INPUT
    %   c: vettore dei coefficienti della funzione obiettivo
    %   Aleq: matrice dei coefficienti delle disequazioni
    %   bleq: vettore dei termini noti delle disequazioni
    %   Aeq: matrice dei coefficienti delle equazioni
    %   beq: vettore dei termini noti delle equazioni
    %   l: vettore dei limiti inferiori
    %   u: vettore dei limiti superiori
    %   OUTPUT
    %   res: struct ottenuta con la funzione simplessoDuale
    
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
    
    % Controllo dimensioni
    if length(c) ~= max(size(l))
        error('Il numero di elementi di c deve essere uguale alla lunghezza di l');
    end
    
    if length(c) ~= max(size(u))
        error('Il numero di elementi di c deve essere uguale alla lunghezza di u');
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
    
    % Converte i limiti delle variabili in vincoli
    A_l = -eye(length(l));
    A_u = eye(length(u));
    b_l = -l;
    b_u = u;
    % Costruisce il problema in forma Primale standard
    c_primale = -c;
    A_primale = [Aleq; Aeq; -Aeq; A_l; A_u];
    b_primale = [bleq; beq; -beq; b_l; b_u];
    
    % Controllo se il problema finale rappresenta effettivamente un problema in forma Primale standard risolvibile
    controlloPrimale(c_primale,A_primale,b_primale);
    
    % Calcola i limiti per il politopo
    [lb,ub] = limitiPrimale(A_primale,b_primale);
    if any(lb == -Inf) || any(ub == Inf)
        error('Il politopo deve essere limitato');
    end
    % Traslo con valori interi il politopo in modo da rendere tutte le variabili definite positive e mantenere la soluzione ottima intera
    traslazione = ceil(-lb);
    c_p = sym(c_primale);
    A_p = sym(A_primale);
    b_p = sym(b_primale+A_primale*traslazione);
    
    % Rendo i vincoli interi in modo da avere anche le variabili di scarto intere
    [A_p,b_p] = semplificaVincoli(A_p,b_p);
    
    [n_vincoli_p, n_variabili_p] = size(A_p);
    
    % visualizzaPolitopo2D(double(c_p),double(A_p),double(b_p));
    % visualizzaPolitopo2D(double(c_primale),double(A_primale),double(b_primale));
    
    count = 1;
    B_p = ausiliarioDuale(A_p,c_p);
    while 1
        fprintf('Iterazione %d\n', count);
        res = simplessoDuale(c_p,A_p,b_p,B_p);
        z = res.x;
        B_p = res.base;
    
        x = z - traslazione;
        fval = c(1:n_variabili_p).'*x;
        % fprintf('z = %s\n', stampaInLinea(z));
        % fprintf('z =~ %s\n', stampaInLinea(vpa(z)));
        fprintf('x = %s\n', stampaInLinea(x));
        fprintf('x =~ %s\n', stampaInLinea(round(x,3)));
        % fprintf('round(x) = %s\n', stampaInLinea(round(x)));
        % fprintf('fract(x) = %s\n', stampaInLinea(parteFrazionaria(x)));
        % dist = norm(x - round(x));
        % fprintf('Euclidean distance =~ %s\n\n', char(vpa(dist)));
        if all(x == round(x))
            break;
        end
    
        % [~,A_d,b_d,B_d] = dualeConMosse(c_p,A_p,b_p,B_p);
        % [~, n_variabili_d] = size(A_d);
        % N_d = setdiff(1:n_variabili_d,B_d);
        n_variabili_d = n_variabili_p+n_vincoli_p;
        A_d = [A_p eye(n_vincoli_p)];
        b_d = b_p;
        N_d = B_p + n_variabili_p;
        B_d = setdiff(1:n_variabili_d,N_d);
    
    
        b_tilde = A_d(:,B_d)\b_d;
        A_tilde = A_d(:,B_d)\A_d(:,N_d);
    
        mask_fract = logical(b_tilde ~= round(b_tilde)); % maschera per prendere le variabili di base non intere
    
        A_tagli_d = sym(zeros(sum(mask_fract), n_variabili_d));
        A_tagli_d(:,N_d) = -parteFrazionaria(A_tilde(mask_fract,:));
    
        b_tagli_d = -parteFrazionaria(b_tilde(mask_fract));
    
        A_tagli_p = A_tagli_d(:,1:n_variabili_p) - A_tagli_d(:,n_variabili_p+1:end)*A_d(:,1:n_variabili_p);
        b_tagli_p = b_tagli_d - A_tagli_d(:,n_variabili_p+1:end)*b_d;
    
        [A_tagli_p,b_tagli_p] = semplificaVincoli(A_tagli_p,b_tagli_p);
    
        % fprintf('B = %s\n', stampaInLinea(B_p));
        % fprintf('xd = %s\n', stampaInLinea(xd));
        % fprintf('Bd = %s\n', stampaInLinea(B_d));
        % fprintf('Nd = %s\n', stampaInLinea(N_d));
        % fprintf('R = %s\n', stampaInLinea(R));
        % fprintf('RIndice = %s\n', stampaInLinea(RIndice));
        % fprintf('Ad =\n');
        % disp(A_d);
        % fprintf('Ad(:,Bd) =\n');
        % disp(A_d(:,B_d));
        % fprintf('Ad(:,Bd)^-1 =\n');
        % disp(A_d(:,B_d)^-1);
        % fprintf('Ad(:,Nd) =\n');
        % disp(A_d(:,N_d));
        % fprintf('ATilde =\n');
        % disp(ATilde);
        % fprintf('ATagliDuale =\n');
        % disp(A_tagli_d);
        % fprintf('bTagliDuale = %s\n', stampaInLinea(b_tagli_d));
        % fprintf('ATagli =\n');
        % disp(A_tagli_p);
        % fprintf('bTagli = %s\n', stampaInLinea(bTagli));
        count = count + 1;
    
        % Tra i piani di taglio, trova quello più distante da z, e lo aggiunge ai vincoli
        distanze = abs(A_tagli_p*z - b_tagli_p).^2./sum(A_tagli_p.^2,2);
        [~, indice] = max(distanze);
    
        
        A_p = [A_p; A_tagli_p(indice,:)];
        b_p = [b_p; b_tagli_p(indice,:)];
        n_vincoli_p = n_vincoli_p + 1;
    
        % visualizzaPianiDiGomory(double(A_tagli_p(indice,:)),double(b_tagli_p(indice,:)));
        % visualizzaPianiDiGomory(double(A_tagli_p(indice,:)),double(b_tagli_p(indice,:) - A_tagli_p(indice,:)*traslazione));
        drawnow;
    end
    
    
    end