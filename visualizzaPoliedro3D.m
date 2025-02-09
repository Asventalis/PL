function visualizzaPoliedro3D(c, A, b)
    % Funzione per visualizzare un poliedro 3D definito da Ax <= b e la direzione del vettore c
    %   PARAMETRI
    %   c: vettore colonna della funzione obiettivo
    %   A: matrice dei vincoli
    %   b: vettore colonna dei termini noti
    if ~exist("c", "var") || ~exist("A", "var") || ~exist("b", "var") 
        error("Errore. c, A e b devono essere tutti e 3 specificati.");
    end
    controlloPrimale(c, A, b);
    [nVincoli,nVariabili] = size(A);

    if nVariabili ~= 3
        error("Errore. Questa funzione funziona solo con problemi a 3 variabili.");
    end
    
    warning('off', 'MATLAB:rankDeficientMatrix');

    c = sym(c);
    A = sym(A);
    b = sym(b);
    vertices = [];
    infinity_starts = [];
    infinity_directions = [];
    vertex_bases = {}; % Per memorizzare le basi dei vertici
    combinations = nchoosek(1:nVincoli, 3);
    for i = 1:size(combinations, 1)
        B = combinations(i, :);
        N = setdiff(1:nVincoli, B);
        if rank(A(B, :)) == 3
            W = - A(B, :)^-1;
            vertex = - W * b(B,:);
            % se il vertice rispetta tutti i vincoli non di base, lo aggiungiamo
            if all(A(N, :) * vertex <= b(N))
                vertices = [vertices; vertex'];
                vertex_bases{end+1} = B; % Memorizza la base
                % Se A(:, N) * W_i ha solo valori negativi, allora aggiungiamo W_i a infinity_directions
                for j = 1:3
                    if all(A(N, :) * W(:, j) <= 0)
                        infinity_starts = [infinity_starts; vertex'];
                        infinity_directions = [infinity_directions; W(:, j)'];
                    end
                end
            end
        else
            % Se il sistema ha infinite soluzioni, prendine una a caso e aggiungila ai vertici
            if rank(A(B, :)) == rank([A(B, :) b(B)])
                % trova soluzione specifica di A(B, :) * x = b(B)
                vertex = A(B, :) \ b(B);
                if all(A(N, :) * vertex <= b(N))
                    vertices = [vertices; vertex'];
                    vertex_bases{end+1} = B; % Memorizza la base
                    W = null(A(B, :));
                    for j = 1:size(W, 2)
                        if all(A(N, :) * W(:, j) >= 0)
                            infinity_starts = [infinity_starts; vertex'];
                            infinity_directions = [infinity_directions; W(:, j)'];
                        end
                    end
                end
            end
        end
    end

    
    % Rimuovi vertici duplicati
    [vertices, unique_idx] = unique(vertices, 'rows');
    vertex_bases = vertex_bases(unique_idx); % Aggiorna le basi per i vertici unici

    scaling_factor = 1;
    % Calcola la lunghezza del vettore c in base alla dimensione del poliedro
    if ~isempty(vertices)
        max_vertex_distance = max(sqrt(sum(vertices.^2, 2)));
        scaling_factor = max_vertex_distance / 5;
    end

    % add a vertex for each infinity direction
    for i = 1:size(infinity_starts, 1)
        start = infinity_starts(i, :);
        direction = infinity_directions(i, :);
        direction = direction / double(norm(direction)) * scaling_factor;
        infinity_starts(i, :) = start + direction;
        vertices = [vertices; start + direction];
        vertex_bases{end+1} = [];
    end

    % Disegna il poliedro
    figure;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Poliedro 3D con vettore c');
    grid on;
    
    % imposta la vista in prospettiva e fai in modo che gli assi non cambino scala
    view(3);
    axis equal;
    daspect([1 1 1]);
    camproj('perspective');
    camva(10); % Set a fixed camera view angle
    hold on;
    % Aggiungi luci da tutte le direzioni
    light_distance = double(5*max_vertex_distance);
    east = [1 0 0]*light_distance;
    front = [0 1 0]*light_distance;
    north = [0 0 1]*light_distance;
    west = -east;
    back = -front;
    south = -north;
    light('Position', east, 'Style', 'local');
    light('Position', north, 'Style', 'local');
    light('Position', front, 'Style', 'local');
    light('Position', west, 'Style', 'local');
    light('Position', south, 'Style', 'local');
    light('Position', back, 'Style', 'local');
    if size(vertices, 1) > 3
        % Trova i convessi del poliedro
        num_vertices = double(vertices);
        K = convhull(num_vertices);
        trisurf(K, num_vertices(:,1), num_vertices(:,2), num_vertices(:,3), 'FaceAlpha', 1, 'EdgeColor', 'k', 'FaceColor', '#005500');
    elseif size(vertices, 1) == 3
        % Disegna il singolo faccia
        fill3(vertices(:,1), vertices(:,2), vertices(:,3), 'k', 'FaceAlpha', 1, 'FaceColor', '#005500');
    elseif size(vertices, 1) == 2
        % Disegna i due punti e la linea tra di loro
        plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'ko', 'MarkerFaceColor', 'k');
        line(vertices(:,1), vertices(:,2), vertices(:,3), 'Color', 'k', 'LineWidth', 2);
    elseif size(vertices, 1) == 1
        % Disegna il singolo vertice
        scatter3(vertices(:,1), vertices(:,2), vertices(:,3), 'filled');
    end
    
    c = c / norm(c) * scaling_factor;

    % Disegna le direzioni all'infinito
    for i = 1:size(infinity_starts, 1)
        start = infinity_starts(i, :);
        direction = infinity_directions(i, :);
        direction = direction / norm(direction) * scaling_factor;
        quiver3(start(1), start(2), start(3), direction(1), direction(2), direction(3), 'b', 'LineWidth', 2, 'MaxHeadSize', 2, 'ShowArrowHead', 'on');
    end
    
    % Aggiungi annotazioni ai vertici
    for i = 1:size(vertices, 1)
        if ~isempty(vertex_bases{i})
            base_str = sprintf('Base: [%s]', num2str(vertex_bases{i}));
            text(vertices(i, 1), vertices(i, 2), vertices(i, 3), base_str, 'FontSize', 8, 'Color', 'r');
        end
    end
    
    % Disegna il vettore c
    quiver3(0, 0, 0, c(1), c(2), c(3), 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    
    hold off;

    warning('on', 'MATLAB:rankDeficientMatrix');
end