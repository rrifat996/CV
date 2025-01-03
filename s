function visibility_graph(agent, goal, polygons)
    % Create a new figure window
    figure;
    hold on;
    
    % Plot obstacles (polygons)
    for i = 1:length(polygons)
        fill(polygons{i}(:, 1), polygons{i}(:, 2), 'k', 'FaceAlpha', 0.5);
    end
    
    % Plot agent (start) and goal nodes
    plot(agent(1), agent(2), 'go', 'MarkerFaceColor', 'g');
    text(agent(1), agent(2), ' Agent', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    
    plot(goal(1), goal(2), 'ro', 'MarkerFaceColor', 'r');
    text(goal(1), goal(2), ' Goal', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    
    % Collect all the vertices of the polygons
    all_vertices = [];
    for i = 1:length(polygons)
        all_vertices = [all_vertices; polygons{i}];
    end
    
    % Check visibility between agent, goal, and vertices
    nodes = [agent; goal; all_vertices];
    
    % Loop through all pairs of nodes and check for line-of-sight
    edges = {};  % Store valid edges to keep track of the lines to plot
    for i = 1:size(nodes, 1)
        for j = i+1:size(nodes, 1)
            if line_of_sight(nodes(i, :), nodes(j, :), polygons)
                % If there's a line-of-sight, store the edge
                edges{end+1} = [nodes(i, :), nodes(j, :)];  % Save the edge as a 4-element vector [x1, y1, x2, y2]
            end
        end
    end
    
    % Plot edges
    for i = 1:length(edges)
        plot([edges{i}(1), edges{i}(3)], [edges{i}(2), edges{i}(4)], 'k--');
    end
    
    % Remove diagonals (edges that do not belong to the boundary of the polygon)
    for i = 1:length(polygons)
        poly = polygons{i};
        n = size(poly, 1);  % Number of vertices in the polygon
        for j = 1:n
            % Get the boundary edges of the polygon
            boundary_start = poly(j, :);
            boundary_end = poly(mod(j, n) + 1, :);
            
            % Identify and remove diagonals (edges not on the boundary)
            for k = length(edges):-1:1
                % Check if the edge corresponds to a diagonal (not on the boundary)
                if is_diagonal(edges{k}, poly)
                    edges(k) = [];  % Remove the diagonal from edges
                end
            end
        end
    end
    
    % Plot updated edges
    for i = 1:length(edges)
        plot([edges{i}(1), edges{i}(3)], [edges{i}(2), edges{i}(4)], 'k-');
    end
    
    % Set the axis limits for better visualization
    axis([0 12 0 12]);
    grid on;
    xlabel('X-axis');
    ylabel('Y-axis');
    title('Visibility Graph with Obstacles');
    
    hold off;
end

% Function to check if there is a line-of-sight between two points
function visible = line_of_sight(p1, p2, polygons)
    visible = true;
    
    % Loop through all polygons to check if the line segment intersects any edge
    for i = 1:length(polygons)
        poly = polygons{i};
        n = size(poly, 1);  % Number of vertices in the polygon
        
        for j = 1:n
            % Define the edge of the polygon
            edge_start = poly(j, :);
            edge_end = poly(mod(j, n) + 1, :);
            
            % Check if the line segment (p1 to p2) intersects the polygon edge
            if line_intersects(p1, p2, edge_start, edge_end)
                visible = false;
                return; % No visibility
            end
        end
    end
end

% Function to check if two line segments intersect
function intersect = line_intersects(p1, p2, p3, p4)
    % Calculate direction of the lines
    denom = (p4(2) - p3(2)) * (p2(1) - p1(1)) - (p4(1) - p3(1)) * (p2(2) - p1(2));
    if denom == 0
        intersect = false; % Lines are parallel or collinear
        return;
    end
    
    % Calculate the intersection point
    ua = ((p4(1) - p3(1)) * (p1(2) - p3(2)) - (p4(2) - p3(2)) * (p1(1) - p3(1))) / denom;
    ub = ((p2(1) - p1(1)) * (p1(2) - p3(2)) - (p2(2) - p1(2)) * (p1(1) - p3(1))) / denom;

    % Check if the intersection point is within the bounds of both segments
    intersect = (ua >= 0 && ua <= 1 && ub >= 0 && ub <= 1);
end

% Function to check if an edge is a diagonal (not part of the polygon's boundary)
function is_diag = is_diagonal(edge, polygon)
    is_diag = true;
    
    % Extract start and end points of the edge
    edge_start = edge(1:2);
    edge_end = edge(3:4);
    
    % Loop through the edges of the polygon to see if the edge is part of the boundary
    n = size(polygon, 1);  % Number of vertices in the polygon
    for i = 1:n
        boundary_start = polygon(i, :);
        boundary_end = polygon(mod(i, n) + 1, :);
        
        % If the edge matches the boundary edge, it's not a diagonal
        if (isequal(edge_start, boundary_start) && isequal(edge_end, boundary_end)) || ...
           (isequal(edge_start, boundary_end) && isequal(edge_end, boundary_start))
            is_diag = false;  % It's a boundary edge, not a diagonal
            break;
        end
    end
end