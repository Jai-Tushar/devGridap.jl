using Gridap
using Gridap.Geometry
using Gridap.Fields
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Helpers
using Gridap.ReferenceFEs

using FillArrays

# @enter CartesianDiscreteModel((0,1,0,1), (2,2))

grid = CartesianGrid((0,1,0,1), (2,2))

@enter _grid =  UnstructuredGrid(grid)

# collect1d(get_node_coordinates(grid))

# struct VoronoiDiscreteModel{D,T,F} <: DiscreteModel{D,D}
#   grid::PhysicalGrid{D,T,F}
#   grid_topology::UnstructuredGridTopology{D,D,T,Oriented}
#   face_labeling::FaceLabeling
# end




function orient_nodes(v)
  _angle(c,v) = atan(v[1]-c[1],v[2]-c[2])
  c = mean(v)
  sortperm(v, by = x -> _angle(c,x), rev = true)
end

# D = 2
# partition = (0,1,0,1)
# cells = (2,2)

# # primal (Simplicial mesh)

# model = simplexify(CartesianDiscreteModel(partition, cells))

# pgrid = get_grid(model)
# ptopo = get_grid_topology(model)
# pface = get_face_labeling(model)


# pn2c = get_faces(ptopo,0,D)
# pn2e = get_faces(ptopo,0,D-1)

# pc2n = get_faces(ptopo,D,0)
# pe2n = get_faces(ptopo,D-1,0)
# pnodes = get_vertex_coordinates(ptopo)

# is_pbnode = Geometry.get_isboundary_face(ptopo,0)
# is_pbface = Geometry.get_isboundary_face(ptopo,D-1)

# pbnodes = findall(is_pbnode)
# pbfaces = findall(is_pbface)
# pnode_reindex = find_inverse_index_map(pbnodes,num_faces(ptopo,0))    # Not counting the interior nodes for the global dual mesh
# pface_reindex = find_inverse_index_map(pbfaces,num_faces(ptopo,D-1))  # Not counting the interior faces for the global dual mesh


# pbulk_cen = map(mean,lazy_map(Broadcasting(Reindex(pnodes)), pc2n))
# pedge_cen = map(mean,lazy_map(Broadcasting(Reindex(pnodes)), pe2n))

# # Dual (Polytopal mesh)

# Dnodes = vcat(pbulk_cen, pedge_cen[is_pbface], pnodes[is_pbnode])
# n_Dnodes = length(Dnodes)
# n_Dcells = num_faces(ptopo,0)

# ptrs = zeros(Int,n_Dcells+1)
# data = zeros(Int,n_Dcells * 3*2^(D-1))    # rough upper-bound

# f_offset = num_faces(ptopo,D)
# n_offset = f_offset + sum(is_pbface)

# ptrs[1] = 1
# for n in 1:n_Dcells
#   if !is_pbnode[n]
#       # Interior node of the dual mesh are given by the neighboring cells of the primal mesh
#       nbr_cells = view(pn2c,n)
#       ptrs[n+1] = ptrs[n] + length(nbr_cells)
#       Dnode_connectivity = nbr_cells
#       n_Dnodes = length(nbr_cells)
#   else
#       # Boundary node of the dual mesh are given by the boundary faces and the boundary nodes of the primal mesh
#       nbr_cells = view(pn2c,n)
#       nbr_faces = pface_reindex[filter(f -> is_pbface[f], view(pn2e,n))]
#       new_node  = pnode_reindex[n]
#       Dnode_connectivity = vcat(nbr_cells, nbr_faces .+ f_offset, new_node + n_offset)
#       n_Dnodes = length(nbr_faces) + length(nbr_cells) + 1
#   end
#   Dcoords = view(Dnodes,Dnode_connectivity)
#   cc_perm = orient_nodes(Dcoords)
#   ptrs[n+1] = ptrs[n] + n_Dnodes
#   data[ptrs[n]:ptrs[n+1]-1] .= Dnode_connectivity[cc_perm]
# end

# Dc2n = Gridap.Arrays.Table(data,ptrs)

# Dcell_types = collect(1:n_Dcells)
# poly = map(Polygon,lazy_map(Broadcasting(Reindex(Dnodes)),Dc2n))

# Dtopo = Geometry.UnstructuredGridTopology(Dnodes,Dc2n,Dcell_types,poly)

# Dcell_map = Fill(identity,n_Dcells)

# facet_normals = map(get_facet_normal, poly)



# Physical grid



# struct PhysicalGrid{D,T,F} <: Grid{D,D}
#   node_coords::PhysicalCoordinates{D,T,F}
#   cell_node_ids::PhysicalCellNodes{D}
#   cell_type::Fill{identity,Tuple{Base.OneTo{Int}}}
#   @doc """
#       PhysicalGrid(desc::CartesianDescriptor)
#   """
#   function PhysicalGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
#     node_coords = CartesianCoordinates(desc)
#     cell_node_ids = CartesianCellNodes(desc.partition)
#     cell_type = Fill(identity,length(cell_node_ids))
#     new{D,T,F}(node_coords,cell_node_ids,cell_type)
#   end
# end


# function VoronoiDiscreteModel(desc::VoronoiDescriptor{D,T,F})
#   grid  = PhysicalGrid(desc)

# Orientation ?
struct PhysicalGrid{Dc,Dp,Tp,O,Tn} <: Grid{Dc,Dp} 
  node_coordinates::Vector{Point{Dp,Tp}}
  cell_node_ids::Table{Int64, Vector{Int64}, Vector{Int64}}
  cell_types::Vector{Int8}
  orientation_style::O
  facet_normal
  cell_map

  function PhysicalGrid(
    node_coordinates::Vector{Point{Dp,Tp}},
    cell_node_ids::Table{Ti},
    cell_types::Vector{Int8},
    orientation_style,
    facet_normal,
    cell_map) where {Dp,Tp,Ti}
    B = typeof(orientation_style)
    Tn = typeof(facet_normal)
    new{Dc,Dp,Tp,B,Tn}(
    node_coordinates,
    cell_node_ids,
    cell_types,
    orientation_style,
    facet_normal,
    cell_map)
  end
end

# @doc """
# function PhysicalGrid(
#   node_coordinates::Vector{Point{Dp,Tp}},
#   cell_node_ids::Table{Ti},
#   cell_types::Vector,
#   orientation_style::OrientationStyle=NonOriented(),
#   facet_normal::Tn) where {Dc,Dp,Tp,Ti}
# """


# """
#     PhysicalGrid(args...;kwargs...)

#   PhysicalGrid born out of  `CartesianDescriptor`
# """

function PhysicalGrid(args...;kwargs...)
  desc = CartesianDescriptor(args...;kwargs...)
  PhysicalGrid(desc)
end

# @doc """
#       PhysicalGrid(desc::VoronoiDescriptor)
# """

function PhysicalGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}

  # Primal simplicial model
    pmodel = simplexify(CartesianDiscreteModel(desc))
    ptopo = get_grid_topology(pmodel)

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

    # Dual (Polygopal mesh)
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
    Dcell_reffes = Fill(identity,n_Dcells)  
    polys = map(Polygon,lazy_map(Broadcasting(Reindex(Dnodes)),Dc2n))
    facet_normals = map(get_facet_normal, polys)
    Dcell_map = Fill(identity,n_Dcells)
    orientation = NonOriented()

    PhysicalGrid(Dnodes, Dc2n, Dcell_types, orientation, facet_normals, Dcell_map) # Dual cell to node connectivity
end

# @enter _grid = PhysicalGrid((0,1,0,1), (2,2))

# Coordinates

# facet_normal


# # Facet measure
# face_coords = get_face_coordinates(Dtopo,2)

# get_face_coordinates(Poly)



####################################

grid = simplexify(CartesianGrid((0,1,0,1), (2,2)))


# function PhysicalGrid(grid::Grid)
#   @assert is_regular(grid) "This parent grid has is not regular"
