
function similar_range(r::Base.OneTo,n::Integer)
  Base.OneTo(Int(n))
end

function similar_range(r::BlockedUnitRange,n::Integer)
  blockedrange([n])
end

function similar_range(r::MultiLevelBlockedUnitRange,n::Integer)
  r = similar_range(first(r.local_ranges),n)
  append_ranges([r])
end

# Restricted for a single non zero block
struct BlockFieldArrayCoo{T,N,A,X} <: AbstractBlockArray{T,N}
  axes::X
  blockids::Vector{NTuple{N,Int}}
  block::A
  function BlockFieldArrayCoo(
    _axes::NTuple{N},
    blockids::Vector{NTuple{N,Int}},
    block::A) where {T,N,A<:AbstractArray{T,N}}

    @assert length(blockids) == 1
    #I = first(blockids)
    #@check blocks_equal(axes(block),map(local_range,_axes,I)) "The given block and axes are incompatible."

    X = typeof(_axes)
    new{T,N,A,X}(_axes,blockids,block)
  end
end

struct BlockFieldArrayCooMap{N} <: Map
  blocksize::NTuple{N,Int}
  blockids::Vector{NTuple{N,Int}}
  ptrs::Array{Int,N}
  function BlockFieldArrayCooMap(blocksize::NTuple{N,Int}, blockids::Vector{NTuple{N,Int}}) where N
    ptrs = fill(-1,blocksize)
    for (p,I) in enumerate(blockids)
      ptrs[I...] = p
      for i in 1:N
        @check blocksize[i] >= I[i]
      end
    end
    new{N}(blocksize,blockids,ptrs)
  end
end

@inline function lazy_map(k::BlockFieldArrayCooMap,T::Type,f::AbstractArray...)
  s = Arrays._common_size(f...)
  N = length(s)
  LazyArray(T,Val(N),Fill(k,s),f...)
end

BlockArrays.blocksize(a::BlockFieldArrayCooMap) = a.blocksize
is_zero_block(a::BlockFieldArrayCooMap,i::Integer) = a.ptrs[i] < 0
is_zero_block(a::BlockFieldArrayCooMap{N},i::Vararg{Integer,N}) where N = a.ptrs[i...] < 0
is_zero_block(a::BlockFieldArrayCooMap{N},i::Vararg{Block,N}) where N = is_zero_block(a,map(Int,i)...)
is_zero_block(a::BlockFieldArrayCooMap,i::Block) = is_zero_block(a,convert(Tuple,i)...)
is_zero_block(a::BlockFieldArrayCooMap,i::CartesianIndex) = is_zero_block(a,Tuple(i)...)

function evaluate!(cache,k::BlockFieldArrayCooMap,axes,block)
  @check map(i->first(blocksize(i)),axes) == k.blocksize "The given axes are not compatible with the given BlockFieldArrayCooMap"
  BlockFieldArrayCoo(axes,k.blockids,block)
end

testitem(a::BlockFieldArrayCoo) = testitem(a.block)

# Specific API

function is_zero_block(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  i != first(a.blockids)
end

# AbstractBlockArray

@inline function BlockArrays.getblock(a::BlockFieldArrayCoo{T,N}, i::Vararg{Integer, N}) where {T,N}
  if i == first(a.blockids)
    a.block
  else
    @notimplemented "Cannot get a zero block from a BlockFieldArrayCoo"
  end
end

# AbstractArray

Base.size(a::BlockFieldArrayCoo) = map(length,Base.axes(a))
Base.axes(a::BlockFieldArrayCoo) = a.axes
Base.IndexStyle(::Type{<:BlockFieldArrayCoo}) = IndexCartesian()

function Base.getindex(a::BlockFieldArrayCoo{T,N},i::Vararg{Integer,N}) where {T,N}
  s = map(findblockindex,a.axes,i)
  ai = a[s...]
  ai
end

# Evaluation

function return_cache(a::BlockFieldArrayCoo,x::Point)
  fc = return_cache(a.block,x)
  fx = return_value(a.block,x)
  bsize = map(i->first(blocksize(i)),a.axes)
  k = BlockArrayCooMap(bsize,a.blockids)
  cr = return_cache(k,a.axes,fx)
  (k,fc,cr)
end

@inline function evaluate!(cache,a::BlockFieldArrayCoo,x::Point)
  k,fc,cr = cache
  fx = evaluate!(fc,a.block,x)
  evaluate!(cr,k,a.axes,fx)
end

function return_cache(a::BlockFieldArrayCoo,x::AbstractVector{<:Point})
  fc = return_cache(a.block,x)
  fx = return_value(a.block,x)
  pr = similar_range(first(a.axes),length(x))
  axs = (pr,a.axes...)
  blockids = map(i->(1,i...),a.blockids)
  bsize = map(i->first(blocksize(i)),axs)
  k = BlockArrayCooMap(bsize,blockids)
  cr = return_cache(k,axs,fx)
  (k,fc,axs,cr)
end

@inline function evaluate!(cache,a::BlockFieldArrayCoo,x::AbstractVector{<:Point})
  k,fc,axs,cr = cache
  fx = evaluate!(fc,a.block,x)
  evaluate!(cr,k,axs,fx)
end

# Gradient

function evaluate!(cache,k::Broadcasting{typeof(∇)},a::BlockFieldArrayCoo)
  g = k(a.block)
  BlockFieldArrayCoo(a.axes,a.blockids,g)
end

function evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::BlockFieldArrayCoo)
  g = k(a.block)
  BlockFieldArrayCoo(a.axes,a.blockids,g)
end

# Transpose

function Base.transpose(a::BlockFieldArrayCoo{T,1} where T)
  r = similar_range(first(axes(a)),1)
  axs = (r,axes(a)...)
  blockids = map(i->(1,i...),a.blockids)
  BlockFieldArrayCoo(axs,blockids,transpose(a.block))
end

# Global optimizations

# Evaluation

function lazy_map(
  ::typeof(evaluate),
  a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},
  x::AbstractArray{<:Point})

  m = a.g.value
  cell_axs, cell_ai = a.f
  cell_aix = lazy_map(evaluate,cell_ai,x)
  lazy_map(BlockArrayCooMap(m.blocksize,m.blockids),cell_axs,cell_aix)
end

function lazy_map(
  ::typeof(evaluate),
  a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},
  cell_x::AbstractArray{<:AbstractVector{<:Point}})

  m = a.g.value
  cell_axs, cell_ai = a.f
  cell_aix = lazy_map(evaluate,cell_ai,cell_x)

  function result_axes_on_point_vector(axs,x)
    pr = similar_range(first(axs),length(x))
    axs = (pr,axs...)
  end

  cell_axs_new = lazy_map(result_axes_on_point_vector,cell_axs,cell_x)

  bsize_new = (1,m.blocksize...)
  bids_new = map(i->(1,i...),m.blockids)

  lazy_map(BlockArrayCooMap(bsize_new,bids_new),cell_axs_new, cell_aix)
end

# Gradient

function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.g.value
  cell_axs, cell_ai = a.f

  cell_gi = lazy_map(k,cell_ai)
  lazy_map(m,cell_axs, cell_gi)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.g.value
  cell_axs, cell_ai = a.f
  cell_gi = lazy_map(k,cell_ai)
  lazy_map(m,cell_axs,cell_gi)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∇)}}})
  b = a.f[1]
  lazy_map(axes,b)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∇∇)}}})
  b = a.f[1]
  lazy_map(axes,b)
end

# Composition

function lazy_map(
  k::Broadcasting{typeof(∘)}, a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}},b::AbstractArray{<:Field})

  m = a.g.value
  cell_axs, cell_ai = a.f
  cell_gi = lazy_map(k,cell_ai,b)
  lazy_map(m,cell_axs,cell_gi)
end

function lazy_map(::typeof(axes),a::LazyArray{<:Fill{Broadcasting{typeof(∘)}}})
  b = a.f[1]
  lazy_map(axes,b)
end

# Transpose

function lazy_map(
  k::typeof(transpose), a::LazyArray{<:Fill{<:BlockFieldArrayCooMap}})

  m = a.g.value
  cell_axs, cell_ai = a.f

  cell_ait = lazy_map(transpose,cell_ai)

  cell_axs_new = lazy_map(cell_axs) do axs
    @check length(axs) == 1
    r = similar_range(first(axs),1)
    (r,axs...)
  end

  bsize_new = (1,m.blocksize...)
  bids_new = map(i->(1,i...),m.blockids)

  lazy_map(BlockFieldArrayCooMap(bsize_new,bids_new),cell_axs_new, cell_ait)
end

# Operations

function _get_axes_and_blocks(f)
  f[1], f[2:end]
end
#
# Unary operations
# Assumption: op is a scaling of a
function lazy_map(
  k::BroadcastingFieldOpMap, a::LazyArray{<:Fill{<:BlockArrayCooMap}})
  m = a.g.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  cell_blocks_new = map(b->lazy_map(k,b),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

# Binary test/field or trial/field

# Assumption: op is a scaling of a
function _op_basis_vs_field(k,a,f)
  @check ndims(eltype(a)) > ndims(eltype(f))
  m = a.g.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  cell_blocks_new = map(b->lazy_map(k,b,f),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{1}}},
  f::AbstractArray{<:Number})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  f::AbstractArray{<:Number})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  f::AbstractArray{<:AbstractVector})
  _op_basis_vs_field(k,a,f)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}},
  f::AbstractArray{<:AbstractVector})
  _op_basis_vs_field(k,a,f)
end

# Assumption: op is a scaling of a
function _op_field_vs_basis(k,f,a)
  @check ndims(eltype(a)) > ndims(eltype(f))
  m = a.g.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  cell_blocks_new = map(b->lazy_map(k,f,b),cell_blocks)
  lazy_map(m,cell_axs,cell_blocks_new...)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:Number},
  a::LazyArray{<:Fill{BlockArrayCooMap{1}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:Number},
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:AbstractVector},
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}})
  _op_field_vs_basis(k,f,a)
end

function lazy_map(
  k::BroadcastingFieldOpMap,
  f::AbstractArray{<:AbstractVector},
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}})
  _op_field_vs_basis(k,f,a)
end

# Binary test/test or trial/trial
# Assumption: op is a linear combination of a and b
function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{N}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{N}}}) where N

  ma = a.g.value
  mb = b.g.value
  @assert blocksize(ma) == blocksize(mb)

  cell_axs, = a.f

  blocks = []
  blockids = eltype(ma.blockids)[]
  for I in eachblockid(ma)
    if is_nonzero_block(ma,I) || is_nonzero_block(mb,I)
      aI = _get_cell_block(ma,a,I)
      bI = _get_cell_block(mb,b,I)
      block = lazy_map(k,aI,bI)
      push!(blocks,block)
      push!(blockids,I.n)
    end
  end

  @assert length(blocks) > 0
  m = BlockArrayCooMap(blocksize(ma),blockids)
  lazy_map(m,cell_axs,blocks...)
end

function _get_cell_block(
  m::BlockArrayCooMap,
  a::AbstractArray{<:BlockArrayCoo{T,N,A}},
  b::Block) where {T,N,A}
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  I = convert(Tuple,b)
  p = m.ptrs[I...]
  if p>0
    return cell_blocks[p]
  else
    return lazy_map(cell_axs) do axs
      laxs = map( local_range, axs, I)
      Arrays._zero_block(A,laxs)
    end
  end
end

# Binary test/trial
# Assumption: op is a product of a and b

#function lazy_map(
#  k::BroadcastingFieldOpMap,
#  a::LazyArray{<:Fill{BlockArrayCooMap{1}}},
#  b::LazyArray{<:Fill{BlockArrayCooMap{2}}})
#
#end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{3}}})

  cell_axs_a, = a.f
  cell_axs_b, = b.f
  cell_axs = lazy_map((a1,a2) -> (a1[1],a1[2],a2[3]),cell_axs_a,cell_axs_b)
  ma = a.g.value
  mb = b.g.value

  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(ma.ptrs,2)
  nfield2 = size(mb.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(ma,I1) && is_nonzero_block(mb,I2)
        aI1 = _get_cell_block(ma,a,I1)
        bI2 = _get_cell_block(mb,b,I2)
        block = lazy_map(k,aI1,bI2)
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end

  @assert length(blocks) > 0
  bs = (ma.blocksize[1],ma.blocksize[2],mb.blocksize[3])
  m = BlockArrayCooMap(bs,blockids)
  lazy_map(m,cell_axs,blocks...)
end

#function lazy_map(
#  k::BroadcastingFieldOpMap,
#  a::LazyArray{<:Fill{BlockArrayCooMap{2}}},
#  b::LazyArray{<:Fill{BlockArrayCooMap{1}}})
#
#end

function lazy_map(
  k::BroadcastingFieldOpMap,
  a::LazyArray{<:Fill{BlockArrayCooMap{3}}},
  b::LazyArray{<:Fill{BlockArrayCooMap{2}}})

  cell_axs_a, = a.f
  cell_axs_b, = b.f
  cell_axs = lazy_map((a1,a2) -> (a2[1],a2[2],a1[3]),cell_axs_a,cell_axs_b)
  ma = a.g.value
  mb = b.g.value

  blocks = []
  blockids = NTuple{3,Int}[]
  nfield1 = size(mb.ptrs,2)
  nfield2 = size(ma.ptrs,3)
  for f1 in 1:nfield1
    I1 = Block(1,f1)
    for f2 in 1:nfield2
      I2 = Block(1,1,f2)
      if is_nonzero_block(mb,I1) && is_nonzero_block(ma,I2)
        bI1 = _get_cell_block(mb,b,I1)
        aI2 = _get_cell_block(ma,a,I2)
        block = lazy_map(k,bI1,aI2)
        push!(blocks,block)
        push!(blockids,(1,f1,f2))
      end
    end
  end

  @assert length(blocks) > 0
  bs = (mb.blocksize[1],mb.blocksize[2],ma.blocksize[3])
  m = BlockArrayCooMap(bs,blockids)
  lazy_map(m,cell_axs,blocks...)
end

# Integration of elem vectors
function lazy_map(k::IntegrationMap,a::LazyArray{<:Fill{BlockArrayCooMap{2}}},w::AbstractArray,j::AbstractArray)
  ma = a.g.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  cell_axs_new = lazy_map(a->(a[2],),cell_axs)
  blocks = map(block->lazy_map(k,block,w,j),cell_blocks)
  blockids = [ (ids[2],) for ids in ma.blockids ]
  m = BlockArrayCooMap(ma.blocksize[2:end],blockids)
  lazy_map(m,cell_axs_new,blocks...)
end

## Integration of elem matrices
function lazy_map(k::IntegrationMap,a::LazyArray{<:Fill{BlockArrayCooMap{3}}},w::AbstractArray,j::AbstractArray)
  ma = a.g.value
  cell_axs, cell_blocks = _get_axes_and_blocks(a.f)
  cell_axs_new = lazy_map(a->(a[2],a[3]),cell_axs)
  blocks = map(block->lazy_map(k,block,w,j),cell_blocks)
  blockids = [ (ids[2],ids[3]) for ids in ma.blockids ]
  m = BlockArrayCooMap(ma.blocksize[2:end],blockids)
  lazy_map(m,cell_axs_new,blocks...)
end
#function lazy_map(k::IntMap,f::VectorOfBlockArrayCoo{T,3} where T,w::AbstractArray,j::AbstractArray)
#  ax = lazy_map(a->(a[2],a[3]),f.axes)
#  blocks = map(block->lazy_map(k,block,w,j),f.blocks)
#  blockids = [ (ids[2], ids[3]) for ids in f.blockids ]
#  VectorOfBlockArrayCoo(blocks,blockids,ax)
#end







## The following masted types are needed to achieve type-stability
## in order to use BlockVectorCoo with arrays of fields
#
#struct MaskedField{F} <: Field
#  field::F
#  mask::Bool
#end
#
#function return_cache(z::MaskedField,x::Point)
#  return_cache(z.field,x)
#end
#
#@inline function evaluate!(cache,z::MaskedField,x::Point) 
#  fx = evaluate!(cache,z.field,x)
#  z.mask ? zero(fx) : fx
#end
#
#testvalue(::Type{MaskedField{F}}) where F = MaskedField(testvalue(F),false)
#
#function return_cache(z::MaskedField,x::AbstractArray{<:Point})
#  return_cache(z.field,x)
#end
#
#function evaluate!(c,z::MaskedField,x::AbstractArray{<:Point})
#  fx = evaluate!(c,z.field,x)
#  if z.mask
#    fill!(fx,zero(eltype(fx)))
#  end
#  fx
#end
#
#@inline gradient(z::MaskedField) = MaskedField(gradient(z.field),z.mask)
#
#struct MaskedFieldArray{T,N,A,X} <: AbstractArray{T,N}
#  axes::X
#  field_array::A
#  mask::Bool
#  function MaskedFieldArray(
#    _axes::NTuple{N},
#    field_array::AbstractArray{F,N},
#    mask::Bool) where {F<:Field,N}
#
#    @check mask || blocks_equal(_axes,axes(field_array)) "Incompatible axes and field_array"
#
#    T = MaskedField{F}
#    A = typeof(field_array)
#    X = typeof(_axes)
#    new{T,N,A,X}(_axes,field_array,mask)
#  end
#end
#
#Base.size(a::MaskedFieldArray) = map(length,Base.axes(a))
#Base.axes(a::MaskedFieldArray) = a.axes
#Base.IndexStyle(::Type{<:MaskedFieldArray}) = IndexCartesian()
#function Base.getindex(a::MaskedFieldArray{T,N},i::Vararg{Integer,N}) where {T,N}
#  if a.mask
#    MaskedField(testitem(a.field_array),a.mask)
#  else
#    MaskedField(a.field_array[i...],a.mask)
#  end
#end
#
##function return_cache(a::MaskedFieldArray,x::Point)
##  cf = return_cache(a.field_array,x)
##  fx = return_value(a.field_array,x)
##  cz = CachedArray(similar(fx,eltype(fx),a.axes))
##  (cf,cz)
##end
##
##@inline function evaluate!(cache,a::MaskedFieldArray,x::Point)
##  cf, cz = cache
##  if a.mask
##    setaxes!(cz,a.axes)
##    r = cz.array
##    fill_entries!(r,zero(eltype(r)))
##    r
##  else
##    evaluate!(cf,a.field_array,x)
##  end
##end
##
##function return_cache(a::MaskedFieldArray,x::AbstractVector{<:Point})
##  cf = return_cache(a.field_array,x)
##  fx = return_value(a.field_array,x)
##  rx = similar_range(first(a.axes),length(x))
##  shape = (,)
##  cz = CachedArray(similar(fx,eltype(fx),a.axes))
##  (cf,cz)
##end
##
##@inline function evaluate!(cache,a::MaskedFieldArray,x::AbstractVector{<:Point})
##  cf, cz = cache
##  if a.mask
##    setaxes!(cz,a.axes)
##    r = cz.array
##    fill_entries!(r,zero(eltype(r)))
##    r
##  else
##    evaluate!(cf,a.field_array,x)
##  end
##end

#function _x_range(x:AbstractVector{<:Point},ran::Tuple{Base.OneTo})
#  Base.OneTo(length(x))
#end
#
#function _x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#  blockedrange([length(x)])
#end
#
#function _x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#_x_range(x:AbstractVector{<:Point},ran::Tuple{BlockedUnitRange})
#  blockedrange([length(x)])
#end
#  blockedrange([length(x)])
#end
#
#function _axes_with_x(x:AbstractVector{<:Point},ran::NTuple{N,<:BlockedUnitRange} where N)
#  np = length(x)
#  (blockedrange([np]),ran...)
#end
#
#function _new_axes(x,ran::NTuple{N,<:MultiLevelBlockedUnitRange} where N)
#  np = length(x)
#  a = _new_axes(x,first())
#  r = blockedrange([np])
#  (blockedrange([r]),ran...)
#end


#
#@inline function evaluate!(cache,a::BlockFieldArray,x::Point)
#  fc, blockids, cr = cache
#  fx = evaluate!(fc,a.field_array,x)
#  blockids[1] = a.blockid
#  evaluate!(cr,BlockVectorCoo,a.axes,blockids,fx)
#end


#struct BlockFieldArray{T,N,A,X} <: AbstractBlockArray{T,N}
#  axes::X
#  blockid::NTuple{N,Int}
#  field_array::A
#  function BlockFieldArray(
#    axes::NTuple{N},
#    blockid::NTuple{N,Int},
#    field_array::AbstractArray{T,N}) where {T<:Field,N}
#
#    @check begin 
#      msg = "The given field_array and axes are incompatible."
#      blocks_equal(axes(field_array),map(local_range,axes,blockid)) msg
#    end
#
#    A = typeof(field_array)
#    X = typeof(axes)
#    new{T,N,A,X}(axes,blockid,field_array)
#  end
#end
#
#Base.size(a::BlockFieldArray) = map(length,Base.axes(a))
#Base.axes(a::BlockFieldArray) = a.axes
#Base.IndexStyle(::Type{<:BlockFieldArray}) = IndexCartesian()
#Base.getindex(a::BlockFieldArray{T,N},i::Vararg{Integer,N}) where {T,N} = @notimplemented
#Base.setindex!(a::BlockFieldArray{T,N},v,i::Vararg{Integer,N}) where {T,N} = @notimplemented
#
#is_zero_block(a::BlockFieldArray{T,N},i::Vararg{Integer,N}) where {T,N} = i!=a.blockid
#BlockArrays.eachblock(a::BlockFieldArray) = ( a[I]  for I in eachblockid(a) )
#function BlockArrays.getblock(a::BlockFieldArray{T,N}, block::Vararg{Integer, N}) where {T,N}
#  if block == a.blockid
#    a.field_array
#  else
#    ai = testitem(a.field_array)
#    laxes = 
#    fill(zero(a),)
#
#  end
#end
#
#function return_cache(a::BlockFieldArray,x::Point)
#  fc = return_cache(a.field_array,x)
#  fx = return_value(a.field_array,x)
#  blockids = [a.blockid]
#  cr = return_cache(BlockVectorCoo,a.axes,blockids,fx)
#  (fc,blockids,cr)
#end
#
#@inline function evaluate!(cache,a::BlockFieldArray,x::Point)
#  fc, blockids, cr = cache
#  fx = evaluate!(fc,a.field_array,x)
#  blockids[1] = a.blockid
#  evaluate!(cr,BlockVectorCoo,a.axes,blockids,fx)
#end
