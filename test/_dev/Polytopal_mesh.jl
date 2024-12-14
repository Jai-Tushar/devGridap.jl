using Gridap
using Gridap.Geometry
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using Gridap.ReferenceFEs

using FillArrays


function orient_nodes(v)
  _angle(c,v) = atan(v[1]-c[1],v[2]-c[2])
  c = mean(v)
  sortperm(v, by = x -> _angle(c,x), rev = true)
end

D = 2
partition = (0,1,0,1)
cells = (2,2)

# primal (Simplicial mesh)

model = simplexify(CartesianDiscreteModel(partition, cells))

pgrid = get_grid(model)
ptopo = get_grid_topology(model)
pface = get_face_labeling(model)


pn2c = get_faces(ptopo,0,D)
pn2e = get_faces(ptopo,0,D-1)

pc2n = get_faces(ptopo,D,0)
pe2n = get_faces(ptopo,D-1,0)
pnodes = get_vertex_coordinates(ptopo)

is_pbnode = Geometry.get_isboundary_face(ptopo,0)
is_pbface = Geometry.get_isboundary_face(ptopo,D-1)

pbnodes = findall(is_pbnode)
pbfaces = findall(is_pbface)
pnode_reindex = find_inverse_index_map(pbnodes,num_faces(ptopo,0))    # Not counting the interior nodes for the global dual mesh
pface_reindex = find_inverse_index_map(pbfaces,num_faces(ptopo,D-1))  # Not counting the interior faces for the global dual mesh


pbulk_cen = map(mean,lazy_map(Broadcasting(Reindex(pnodes)), pc2n))
pedge_cen = map(mean,lazy_map(Broadcasting(Reindex(pnodes)), pe2n))

# Dual (Polytopal mesh)

Dnodes = vcat(pbulk_cen, pedge_cen[is_pbface], pnodes[is_pbnode])
n_Dnodes = length(Dnodes)
n_Dcells = num_faces(ptopo,0)

ptrs = zeros(Int,n_Dcells+1)
data = zeros(Int,n_Dcells * 3*2^(D-1))    # rough upper-bound

f_offset = num_faces(ptopo,D)
n_offset = f_offset + sum(is_pbface)

ptrs[1] = 1
for n in 1:n_Dcells
  if !is_pbnode[n]
      # Interior node of the dual mesh are given by the neighboring cells of the primal mesh
      nbr_cells = view(pn2c,n)
      ptrs[n+1] = ptrs[n] + length(nbr_cells)
      Dnode_connectivity = nbr_cells
      n_Dnodes = length(nbr_cells)
  else
      # Boundary node of the dual mesh are given by the boundary faces and the boundary nodes of the primal mesh
      nbr_cells = view(pn2c,n)
      nbr_faces = pface_reindex[filter(f -> is_pbface[f], view(pn2e,n))]
      new_node  = pnode_reindex[n]
      Dnode_connectivity = vcat(nbr_cells, nbr_faces .+ f_offset, new_node + n_offset)
      n_Dnodes = length(nbr_faces) + length(nbr_cells) + 1
  end
  Dcoords = view(Dnodes,Dnode_connectivity)
  cc_perm = orient_nodes(Dcoords)
  ptrs[n+1] = ptrs[n] + n_Dnodes
  data[ptrs[n]:ptrs[n+1]-1] .= Dnode_connectivity[cc_perm]
end

Dc2n = Gridap.Arrays.Table(data,ptrs)

Dcell_types = collect(1:n_Dcells)
poly = map(Polygon,lazy_map(Broadcasting(Reindex(Dnodes)),Dc2n))

Dtopo = Geometry.UnstructuredGridTopology(Dnodes,Dc2n,Dcell_types,poly)

Dcell_map = Fill(identity,n_Dcells)



# In the physical space

# Facet measure
face_coords = get_face_coordinates(Dtopo,2)

get_face_coordinates(Poly)