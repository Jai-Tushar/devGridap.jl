
"""
    struct QGradLegendrePolynomialBasis{...} <: AbstractVector{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Nedelec reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QGradLegendrePolynomialBasis{D,T} <: AbstractVector{Legendre}
  order::Int
  terms::CartesianIndices{D}
  perms::Matrix{Int}
  function QGradLegendrePolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    new{D,T}(order,terms,perms)
  end
end

Base.size(a::QGradLegendrePolynomialBasis) = (_ndofs_qgrad_leg(a),)
Base.getindex(a::QGradLegendrePolynomialBasis,i::Integer) = Legendre()
Base.IndexStyle(::QGradLegendrePolynomialBasis) = IndexLinear()

"""
    QGradLegendrePolynomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QGradLegendrePolynomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the curl of the  functions in this basis
is in the Q space of degree `order`.
"""
function QGradLegendrePolynomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order + 1
  _t = tfill(_order+1,Val{D-1}())
  t = (_order,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QGradLegendrePolynomialBasis(T,order,terms,perms)
end

"""
    num_terms(f::QGradLegendrePolynomialBasis{D,T}) where {D,T}
"""
num_terms(f::QGradLegendrePolynomialBasis{D,T}) where {D,T} = length(f.terms)*D

get_order(f::QGradLegendrePolynomialBasis) = f.order

return_type(::QGradLegendrePolynomialBasis{D,T}) where {D,T} = VectorValue{D,T}

function return_cache(f::QGradLegendrePolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad_leg(f)
  n = 1 + f.order+1
  V = VectorValue{D,T}
  r = CachedArray(zeros(V,(np,ndof)))
  v = CachedArray(zeros(V,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  (r, v, c)
end

function evaluate!(cache,f::QGradLegendrePolynomialBasis{D,T},x::AbstractVector{<:Point}) where {D,T}
  r, v, c = cache
  np = length(x)
  ndof = _ndofs_qgrad_leg(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  for i in 1:np
    @inbounds xi = x[i]
      _evaluate_nd_qgrad_leg!(v,xi,f.order+1,f.terms,f.perms,c)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

function return_cache(
  fg::FieldGradientArray{1,QGradLegendrePolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  @check D == length(eltype(x)) "Incorrect number of point components"
  np = length(x)
  ndof = _ndofs_qgrad_leg(f)
  n = 1 + f.order+1
  xi = testitem(x)
  V = VectorValue{D,T}
  G = gradient_type(V,xi)
  r = CachedArray(zeros(G,(np,ndof)))
  v = CachedArray(zeros(G,(ndof,)))
  c = CachedArray(zeros(T,(D,n)))
  g = CachedArray(zeros(T,(D,n)))
  (r, v, c, g)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{1,QGradLegendrePolynomialBasis{D,T}},
  x::AbstractVector{<:Point}) where {D,T}

  f = fg.fa
  r, v, c, g = cache
  np = length(x)
  ndof = _ndofs_qgrad_leg(f)
  n = 1 + f.order+1
  setsize!(r,(np,ndof))
  setsize!(v,(ndof,))
  setsize!(c,(D,n))
  setsize!(g,(D,n))
  V = VectorValue{D,T}
  for i in 1:np
    @inbounds xi = x[i]
    _gradient_nd_qgrad_leg!(v,xi,f.order+1,f.terms,f.perms,c,g,V)
    for j in 1:ndof
      @inbounds r[i,j] = v[j]
    end
  end
  r.array
end

# Helpers

_ndofs_qgrad_leg(f::QGradLegendrePolynomialBasis{D}) where D = D*(length(f.terms))

function _evaluate_nd_qgrad_leg!(
  v::AbstractVector{V},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T}) where {V,T,D}

  dim = D
  for d in 1:dim
    _evaluate_1d!(Legendre,order,c,x,d)
  end

  o = one(T)
  k = 1
  m = zero(Mutable(V))
  js = eachindex(m)
  z = zero(T)

  for ci in terms

    for j in js

      @inbounds for i in js
        m[i] = z
      end

      s = o
      @inbounds for d in 1:dim
        s *= c[d,ci[perms[d,j]]]
      end

      m[j] = s
      v[k] = m
      k += 1

    end

  end

end

function _gradient_nd_qgrad_leg!(
  v::AbstractVector{G},
  x,
  order,
  terms::CartesianIndices{D},
  perms::Matrix{Int},
  c::AbstractMatrix{T},
  g::AbstractMatrix{T},
  ::Type{V}) where {G,T,D,V}

  dim = D
  for d in 1:dim
    _derivatives_1d!(Legendre,order,(c,g),x,d)
  end

  z = zero(Mutable(V))
  m = zero(Mutable(G))
  js = eachindex(z)
  mjs = eachindex(m)
  o = one(T)
  zi = zero(T)
  k = 1

  for ci in terms

    for j in js

      s = z
      for i in js
        s[i] = o
      end

      for q in 1:dim
        for d in 1:dim
          if d != q
            @inbounds s[q] *= c[d,ci[perms[d,j]]]
          else
            @inbounds s[q] *= g[d,ci[perms[d,j]]]
          end
        end
      end

      @inbounds for i in mjs
        m[i] = zi
      end

      for i in js
        @inbounds m[i,j] = s[i]
      end
        @inbounds v[k] = m
      k += 1

    end

  end

end

############################################################################################

"""
    struct QCurlGradLegendrePolynomialBasis{...} <: AbstractArray{Monomial}

This type implements a multivariate vector-valued polynomial basis
spanning the space needed for Raviart-Thomas reference elements on n-cubes.
The type parameters and fields of this `struct` are not public.
This type fully implements the [`Field`](@ref) interface, with up to first order
derivatives.
"""
struct QCurlGradLegendrePolynomialBasis{D,T} <: AbstractVector{Legendre}
  qgrad::QGradLegendrePolynomialBasis{D,T}
  function QCurlGradLegendrePolynomialBasis(::Type{T},order::Int,terms::CartesianIndices{D},perms::Matrix{Int}) where {D,T}
    qgrad = QGradLegendrePolynomialBasis(T,order,terms,perms)
    new{D,T}(qgrad)
  end
end

Base.size(a::QCurlGradLegendrePolynomialBasis) = (length(a.qgrad),)
Base.getindex(a::QCurlGradLegendrePolynomialBasis,i::Integer) = Legendre()
Base.IndexStyle(::QCurlGradLegendrePolynomialBasis) = IndexLinear()

"""
    QCurlGradLegendrePolynomialBasis{D}(::Type{T},order::Int) where {D,T}

Returns a `QCurlGradLegendrePolynomialBasis` object. `D` is the dimension
of the coordinate space and `T` is the type of the components in the vector-value.
The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function QCurlGradLegendrePolynomialBasis{D}(::Type{T},order::Int) where {D,T}
  @check T<:Real "T needs to be <:Real since represents the type of the components of the vector value"
  _order = order+1
  _t = tfill(_order,Val{D-1}())
  t = (_order+1,_t...)
  terms = CartesianIndices(t)
  perms = _prepare_perms(D)
  QCurlGradLegendrePolynomialBasis(T,order,terms,perms)
end

return_type(::QCurlGradLegendrePolynomialBasis{D,T}) where {D,T} = VectorValue{D,T}

function return_cache(f::QCurlGradLegendrePolynomialBasis,x::AbstractVector{<:Point})
  return_cache(f.qgrad,x)
end

function evaluate!(cache,f::QCurlGradLegendrePolynomialBasis,x::AbstractVector{<:Point})
  evaluate!(cache,f.qgrad,x)
end

function return_cache(
  fg::FieldGradientArray{N,<:QCurlGradLegendrePolynomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  return_cache(FieldGradientArray{N}(f.qgrad),x)
end

function evaluate!(
  cache,
  fg::FieldGradientArray{N,<:QCurlGradLegendrePolynomialBasis},
  x::AbstractVector{<:Point}) where N

  f = fg.fa
  evaluate!(cache,FieldGradientArray{N}(f.qgrad),x)
end

"""
    num_terms(f::QCurlGradLegendrePolynomialBasis{D,T}) where {D,T}
"""
num_terms(f::QCurlGradLegendrePolynomialBasis{D,T}) where {D,T} = length(f.qgrad.terms)*D

get_order(f::QCurlGradLegendrePolynomialBasis{D,T}) where {D,T} = get_order(f.qgrad)

############################################################################################
