
"""
    abstract type Pushforward <: Map end

Represents a pushforward map F*, defined as 
  F* : V̂ -> V
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.
"""
abstract type Pushforward <: Map end

abstract type PushforwardRefFE <: ReferenceFEName end

Pushforward(::Type{<:PushforwardRefFE}) = @abstractmethod
Pushforward(name::PushforwardRefFE) = Pushforward(typeof(name))

function Arrays.lazy_map(
  k::Pushforward, ref_cell_fields::AbstractArray, pf_args::AbstractArray...
)
  lazy_map(Broadcasting(Operation(k)), ref_cell_fields, pf_args...)
end

function Arrays.evaluate!(
  cache, k::Pushforward, v_ref::Number, args...
)
  @abstractmethod
end

function evaluate!(
  cache, k::Pushforward, f_ref::AbstractVector{<:Field}, args...
)
  Broadcasting(Operation(k))(f_ref,args...)
end

function Arrays.lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{<:Pushforward}}}}
)
  cell_ref_fields, args... = a.args
  cell_ref_gradient = lazy_map(Broadcasting(∇),cell_ref_fields)
  return lazy_map(a.maps.value,cell_ref_gradient,args...)
end

function Arrays.evaluate!(
  cache,
  ::Broadcasting{typeof(gradient)},
  a::Fields.BroadcastOpFieldArray{<:Pushforward}
)
  v, pf_args... = a.args
  grad_v = Broadcasting(∇)(v)
  Broadcasting(Operation(a.op))(grad_v,pf_args...)
end

# InversePushforward

"""
    const InversePushforward{PF} = InverseMap{PF} where PF <: Pushforward

Represents the inverse of a pushforward map F*, defined as 
  (F*)^-1 : V -> V̂
where 
  - V̂ is a function space on the reference cell K̂ and 
  - V is a function space on the physical cell K.

# Note: 

Given a pushforward F*, we provide a default implementation for the inverse (F*)^-1 that 
is not optimal (and I think wrong for nonlinear geometries???). 

For better performance, the user should overload the following methods for each 
specific pushforward type:  

- `return_cache(k::InversePushforward{PF}, v_phys::Number, args...)`
- `evaluate!(cache, k::InversePushforward{PF}, v_phys::Number, args...)`

"""
const InversePushforward{PF} = InverseMap{PF} where PF <: Pushforward

function Arrays.lazy_map(
  k::InversePushforward, phys_cell_fields::AbstractArray, pf_args::AbstractArray...
)
  lazy_map(Broadcasting(Operation(k)), phys_cell_fields, pf_args...)
end

function Arrays.return_cache(
  k::InversePushforward, v_phys::Number, args...
)
  mock_basis(::VectorValue{D,T}) where {D,T} = one(TensorValue{D,D,T})
  v_ref_basis = mock_basis(v_phys)
  pf_cache = return_cache(inverse_map(k),v_ref_basis,args...)
  return v_ref_basis, pf_cache
end

function Arrays.evaluate!(
  cache, k::InversePushforward, v_phys::Number, args...
)
  v_ref_basis, pf_cache = cache
  change = evaluate!(pf_cache,inverse_map(k),v_ref_basis,args...)
  return v_phys⋅inv(change)
end

function evaluate!(
  cache, k::InversePushforward, f_phys::AbstractVector{<:Field}, args...
)
  Broadcasting(Operation(k))(f_phys,args...)
end

function evaluate!(
  cache, k::InversePushforward, f_phys::Field, args...
)
  Operation(k)(f_phys,args...)
end

# Pullback

"""
    struct Pullback{PF <: Pushforward} <: Map end

Represents a pullback map F**, defined as 
  F** : V* -> V̂*
where 
  - V̂* is a dof space on the reference cell K̂ and 
  - V* is a dof space on the physical cell K.
Its action on physical dofs σ : V -> R is defined in terms of the pushforward map F* as
  ̂σ = F**(σ) := σ∘F* : V̂ -> R
"""
struct Pullback{PF <: Pushforward} <: Map
  pushforward::PF
end

function Arrays.lazy_map(
  ::typeof(evaluate),k::LazyArray{<:Fill{<:Pullback}},ref_cell_fields::AbstractArray
)
  pb = k.maps.value
  phys_cell_dofs, pf_args... = k.args
  phys_cell_fields = lazy_map(pb.pushforward,ref_cell_fields,pf_args...)
  return lazy_map(evaluate,phys_cell_dofs,phys_cell_fields)
end

function evaluate!(
  cache, k::Pullback, σ_phys::AbstractVector{<:Dof}, args...
)
  return MappedDofBasis(k.pushforward,σ_phys,args...)
end

# InversePullback

"""
    struct InversePullback{PF <: Pushforward} <: Map end

Represents the inverse of the pullback map F**, defined as 
  (F**)^-1 : V̂* -> V*
where 
  - V̂* is a dof space on the reference cell K̂ and 
  - V* is a dof space on the physical cell K.
Its action on reference dofs ̂σ : V̂ -> R is defined in terms of the pushforward map F* as
  σ = (F**)^-1(̂σ) := ̂σ∘(F*)^-1 : V -> R
"""
const InversePullback{PB} = InverseMap{PB} where PB <: Pullback

function Arrays.lazy_map(
  ::typeof(evaluate), k::LazyArray{<:Fill{<:InversePullback}}, phys_cell_fields::AbstractArray
)
  pb = inverse_map(k.maps.value)
  ref_cell_dofs, pf_args... = k.args
  ref_cell_fields = lazy_map(inverse_map(pb.pushforward), phys_cell_fields, pf_args...)
  return lazy_map(evaluate,ref_cell_dofs,ref_cell_fields)
end

function evaluate!(
  cache, k::InversePullback, σ_ref::AbstractVector{<:Dof}, args...
)
  pb = inverse_map(k)
  return MappedDofBasis(inverse_map(pb.pushforward),σ_ref,args...)
end

# ContraVariantPiolaMap

struct ContraVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::ContraVariantPiolaMap, v_ref::Number, Jt::Number
)
  idetJ = 1. / meas(Jt)
  return v_ref ⋅ (idetJ * Jt)
end

function return_cache(
  ::InversePushforward{ContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  nothing
end

function evaluate!(
  cache, ::InversePushforward{ContraVariantPiolaMap}, v_phys::Number, Jt::Number
)
  detJ = meas(Jt)
  return v_phys ⋅ (detJ * pinvJt(Jt))
end

# TODO: Should this be here? Probably not...

function Fields.DIV(f::LazyArray{<:Fill})
  df = Fields.DIV(f.args[1])
  k  = f.maps.value
  lazy_map(k,df)
end

function Fields.DIV(a::LazyArray{<:Fill{typeof(linear_combination)}})
  i_to_basis  = Fields.DIV(a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

function Fields.DIV(f::LazyArray{<:Fill{Broadcasting{Operation{ContraVariantPiolaMap}}}})
  ϕrgₖ = f.args[1]
  return lazy_map(Broadcasting(divergence),ϕrgₖ)
end

function Fields.DIV(f::Fill{<:Fields.BroadcastOpFieldArray{ContraVariantPiolaMap}})
  ϕrgₖ = f.value.args[1]
  return Fill(Broadcasting(divergence)(ϕrgₖ),length(f))
end

# CoVariantPiolaMap

struct CoVariantPiolaMap <: Pushforward end

function evaluate!(
  cache, ::CoVariantPiolaMap, v_ref::Number, Jt::Number
)
  # we right-multiply to compute the gradient correctly
  return v_ref ⋅ transpose(inv(Jt))
end

function return_cache(
  ::InversePushforward{CoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return nothing
end

function evaluate!(
  cache, ::InversePushforward{CoVariantPiolaMap}, v_phys::Number, Jt::Number
)
  return v_phys ⋅ transpose(Jt)
end
