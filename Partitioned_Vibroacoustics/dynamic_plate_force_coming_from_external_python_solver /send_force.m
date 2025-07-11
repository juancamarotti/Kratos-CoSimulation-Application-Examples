function force = send_force(num_nodes)
    % Example: apply a constant force in Z direction
    force = repmat([0.0, 0.0, -0.01], num_nodes, 1);
    force = reshape(force', 1, []);  % Return as flat row vector
end