function [nodes,triangles] = make_triangle_mesh(depth)

nodes     = eye(3);
n_nodes   = 3;
triangles = [1 2 3]';

for it_depth = 1:depth,
  my_triangles = [];
  n_triangles = size(triangles,2);
  for it = 1:n_triangles,
    old_indices = triangles(:,it);
    new_indices = n_nodes + [1 2 3];
    old_nodes   = nodes(:,old_indices);
    new_nodes   = 0.5 * [ old_nodes(:,2)+old_nodes(:,3), old_nodes(:,3)+old_nodes(:,1), old_nodes(:,1)+old_nodes(:,2)];
    nodes   = [nodes, new_nodes];
    n_nodes = n_nodes + 3;
    new_triangles =  [[old_indices(1); new_indices(3); new_indices(2)],...
                      [old_indices(2); new_indices(1); new_indices(3)],...
                      [old_indices(3); new_indices(2); new_indices(1)],...
                      [new_indices(1); new_indices(2); new_indices(3)]];
    my_triangles = [my_triangles, new_triangles];
  end
  triangles = my_triangles;
end
