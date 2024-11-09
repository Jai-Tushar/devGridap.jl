
# MomentBasedDofBasis

struct Moment <: Dof end

struct MomentBasedDofBasis{P,V} <: AbstractVector{Moment}
  nodes::Vector{P}
  face_moments::Vector{Array{V}}
  face_nodes::Vector{UnitRange{Int}}

  function MomentBasedDofBasis(nodes,f_moments,f_nodes)
    P = eltype(nodes)
    V = eltype(eltype(f_moments))
    new{P,V}(nodes,f_moments,f_nodes)
  end

  function MomentBasedDofBasis(f_nodes,f_moments)
    P = eltype(eltype(f_nodes))
    V = eltype(eltype(f_moments))
    nodes = P[]
    face_nodes = UnitRange{Int}[]
    nfaces = length(f_nodes)
    n = 1
    for fi in 1:nfaces
      nodes_fi = f_nodes[fi]
      nini = n
      for node_fi in nodes_fi
        push!(nodes,node_fi)
        n += 1
      end
      nend = n-1
      push!(face_nodes,nini:nend)
    end
    new{P,V}(nodes,f_moments,face_nodes)
  end
end

Base.size(a::MomentBasedDofBasis) = (length(a.nodes),)
Base.axes(a::MomentBasedDofBasis) = (axes(a.nodes,1),)
Base.getindex(a::MomentBasedDofBasis,i::Integer) = Moment()
Base.IndexStyle(::MomentBasedDofBasis) = IndexLinear()

get_nodes(b::MomentBasedDofBasis) = b.nodes
get_face_moments(b::MomentBasedDofBasis) = b.face_moments
get_face_nodes_dofs(b::MomentBasedDofBasis) = b.face_nodes

function num_dofs(b::MomentBasedDofBasis)
  n = 0
  for m in b.face_moments
    n += size(m,2)
  end
  n
end

function return_cache(b::MomentBasedDofBasis{P,V},field) where {P,V}
  alloc_cache(vals::AbstractVector,T,ndofs) = zeros(T,ndofs)
  alloc_cache(vals::AbstractMatrix,T,ndofs) = zeros(T,ndofs,size(vals,2))
  cf = return_cache(field,b.nodes)
  vals = evaluate!(cf,field,b.nodes)
  T = typeof(dot(zero(V),zero(eltype(vals))))
  r = alloc_cache(vals,T,num_dofs(b))
  c = CachedArray(r)
  (c, cf)
end

function evaluate!(cache,b::MomentBasedDofBasis,field)
  c, cf = cache
  vals = evaluate!(cf,field,b.nodes)
  dofs = c.array
  _eval_moment_dof_basis!(dofs,vals,b)
  dofs
end

function _eval_moment_dof_basis!(dofs,vals::AbstractVector,b)
  o = 1
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in eachindex(face_moments)
    moments = face_moments[face]
    if !iszero(length(moments))
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        dofs[o] = z
        for i in 1:ni
          dofs[o] += moments[i,j]⋅vals[nodes[i]]
        end
        o += 1
      end
    end
  end
end

function _eval_moment_dof_basis!(dofs,vals::AbstractMatrix,b)
  o = 1
  na = size(vals,2)
  z = zero(eltype(dofs))
  face_nodes = b.face_nodes
  face_moments = b.face_moments
  for face in eachindex(face_moments)
    moments = face_moments[face]
    if !iszero(length(moments))
      nodes = face_nodes[face]
      ni,nj = size(moments)
      for j in 1:nj
        for a in 1:na
          dofs[o,a] = z
          for i in 1:ni
            dofs[o,a] += moments[i,j]⋅vals[nodes[i],a]
          end
        end
        o += 1
      end
    end
  end
end

# MomentBasedReferenceFE

mutable struct FaceMeasure{Df,Dc}
  face ::Int
  cpoly::Polytope{Dc}
  fpoly::Polytope{Df}
  quad ::Quadrature
  fmaps::Vector{<:Field}
  function FaceMeasure(
    cpoly::Polytope{Dc},fpoly::Polytope{Df},order::Int
  ) where {Df,Dc}
    # Quadrature on the face
    quad = Quadrature(fpoly,order)
    # Face to cell coordinate map
    if Df == Dc
      fmaps = [GenericField(identity)]
    else # TODO: Could this be an AffineMap?
      fcoords = get_face_coordinates(cpoly,Df)
      basis = get_shapefuns(LagrangianRefFE(Float64,fpoly,1))
      fmaps = map(c -> linear_combination(c,basis),fcoords)
    end
    new{Df,Dc}(1,cpoly,fpoly,quad,fmaps)
  end
end

function set_face!(m::FaceMeasure{Df},face::Int) where {Df}
  @assert 0 < face <= num_faces(m.cpoly, Df)
  m.face = face
  return m
end

# TODO: Normals are accesed, but tangent are computed on demand. This means 
# that we will be repeating work unless we cache them.
function get_facet_normal(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  n = get_facet_normal(m.cpoly)
  return ConstantField(n[m.face])
end

function get_edge_tangent(m::FaceMeasure{1,Dc}) where {Dc}
  t = get_edge_tangent(m.cpoly)
  return ConstantField(t[m.face])
end

function get_extension(m::FaceMeasure{Df,Dc}) where {Df,Dc}
  @assert Df == Dc - 1
  vs = ReferenceFEs._nfaces_vertices(Float64,m.cpoly,Df)[m.face]
  return ConstantField(TensorValue(hcat([vs[2]-vs[1]...],[vs[3]-vs[1]...])))
end

function Arrays.return_cache(
  σ::Function,                # σ(φ,μ,ds) -> Field/Array{Field}
  φ::AbstractArray{<:Field},  # φ: prebasis (defined on the cell)
  μ::AbstractArray{<:Field},  # μ: polynomial basis (defined on the face)
  ds::FaceMeasure             # ds: face measure
)
  fmap = ds.fmaps[ds.face]
  φf = transpose(Broadcasting(Operation(∘))(φ,fmap))
  f = σ(φf,μ,ds)

  xf = get_coordinates(ds.quad)
  w = get_weights(ds.quad)
  fmap_cache = return_cache(fmap,xf)
  
  f_cache = return_cache(f,xf)
  return fmap_cache, f_cache, xf, w
end

function Arrays.evaluate!(
  cache,
  σ::Function,                # σ(φ,μ,ds) -> Field/Array{Field}
  φ::AbstractArray{<:Field},  # φ: prebasis (defined on the cell)
  μ::AbstractArray{<:Field},  # μ: polynomial basis (defined on the face)
  ds::FaceMeasure             # ds: face measure
)
  fmap_cache, f_cache, xf, w = cache

  fmap = ds.fmaps[ds.face]
  φf = transpose(Broadcasting(Operation(∘))(φ,fmap))
  f = σ(φf,μ,ds)

  xc = evaluate!(fmap_cache,fmap,xf) # quad pts on the cell
  fx = evaluate!(f_cache,f,xf) # f evaluated on the quad pts
  fx .= w .* fx
  return fx, xc
end

component_basis(T::Type{<:Real}) = [one(T)]
function component_basis(V::Type{<:MultiValue})
  T = eltype(V)
  n = num_components(V)
  z, o = zero(T), one(T)
  return [V(ntuple(i -> ifelse(i == j, o, z),Val(n))) for j in 1:n]
end

"""
A moment is given by a triplet (f,σ,μ) where 
  - f is vector of ids if faces Fk
  - σ is a function σ(φ,μ,ds) that returns a Field-like object to be integrated over each Fk
  - μ is a polynomials basis on Fk

We are assuming that all the faces in a moment are of the same type.
"""
function MomentBasedReferenceFE(
  name::ReferenceFEName,
  p::Polytope{D},
  prebasis::AbstractVector{<:Field},
  moments::AbstractVector{<:Tuple{Vector{Int},<:Function,<:AbstractArray{<:Field}}},
  conformity::Conformity;
) where D

  n_faces = num_faces(p)
  n_moments = length(moments)
  face_dims = get_facedims(p)
  face_offsets = get_offsets(p)
  reffaces, face_types = _compute_reffaces_and_face_types(p)

  T = return_type(prebasis)
  order = get_order(prebasis)
  φ_vec = component_basis(T)
  φ = map(constant_field,φ_vec)

  # Create face measures for each moment
  measures = Vector{FaceMeasure}(undef,n_moments)
  for (k,(faces,σ,μ)) in enumerate(moments)
    ftype = face_types[first(faces)]
    @assert all(face_types[faces] .== ftype)
    qdegree = order + get_order(μ) + 1
    fp = reffaces[ftype]
    measures[k] = FaceMeasure(p,fp,qdegree)
  end

  # Count number of dofs and quad pts per face
  face_n_dofs = zeros(Int,n_faces)
  face_n_nodes = zeros(Int,n_faces)
  for ((faces,σ,μ),ds) in zip(moments,measures)
    face_n_dofs[faces] += size(μ,2)
    face_n_nodes[faces] += num_points(ds.quad)
  end
  
  # Compute face moment and node indices
  n_dofs = 1
  n_nodes = 1
  face_own_dofs = Vector{Vector{Int}}(undef,n_faces)
  face_nodes = Vector{UnitRange{Int}}(undef, n_faces)
  face_moments = Vector{Array{T}}(undef, n_faces)
  for face in 1:n_faces
    n_dofs_i = face_n_dofs[face]
    n_nodes_i = face_n_nodes[face]
    face_own_dofs[face] = collect(n_dofs:(n_dofs+n_dofs_i-1))
    face_nodes[face] = n_nodes:(n_nodes+n_nodes_i-1)
    face_moments[face] = zeros(T,n_nodes_i,n_dofs_i)
    n_dofs += n_dofs_i
    n_nodes += n_nodes_i
  end
  nodes = Vector{Point{D}}(undef,n_nodes)

  # Compute face moments and nodes
  fill!(face_n_dofs,0)
  fill!(face_n_nodes,0)
  for ((faces,σ,μ),ds) in zip(moments,measures)
    cache = return_cache(σ,φ,μ,ds)
    for face in faces
      d = face_dims[face]
      lface = face - face_offsets[d]
      set_face!(ds,lface)

      # vals : (nN, nμ, nφ), coords : (nN)
      vals, coords = evaluate!(cache,σ,φ,μ,ds)

      dof_offset = face_n_dofs[face]
      node_offset = first(face_nodes[face]) + face_n_nodes[face]
      for i in axes(vals,1)
        for j in axes(vals,2)
          face_moments[face][i,j+dof_offset] = dot(vals[i,j,:],φ_vec)
        end
        nodes[i+node_offset] = coords[i]
      end

      face_n_nodes[face] += size(vals,1)
      face_n_dofs[face] += size(vals,2)
    end
  end

  dof_basis = MomentBasedDofBasis(nodes, face_moments, face_nodes)
  metadata = nothing
  return GenericRefFE{typeof(name)}(
    n_dofs, p, prebasis, dof_basis, conformity, metadata, face_own_dofs
  )
end
