
struct DivConformity <: Conformity end
abstract type DivConforming <: ReferenceFEName end

# RaviartThomas

struct RaviartThomas <: DivConforming end
const raviart_thomas = RaviartThomas()

"""
    RaviartThomasRefFE(::Type{et},p::Polytope,order::Integer) where et

The `order` argument has the following meaning: the divergence of the  functions in this basis
is in the Q space of degree `order`.
"""
function RaviartThomasRefFE(
  ::Type{T},p::Polytope{D},order::Integer
) where {T,D}

  if is_n_cube(p)
    prebasis = QCurlGradJacobiPolynomialBasis{D}(T,order)          # Prebasis
    cb = QGradJacobiPolynomialBasis{D}(T,order-1)                  # Cell basis
    fb = JacobiPolynomialBasis{D-1}(T,order,Polynomials._q_filter) # Face basis
  elseif is_simplex(p)
    prebasis = PCurlGradMonomialBasis{D}(T,order)                                 # Prebasis
    cb = JacobiPolynomialBasis{D}(VectorValue{D,T},order-1,Polynomials._p_filter) # Cell basis
    fb = JacobiPolynomialBasis{D-1}(T,order,Polynomials._p_filter)                # Face basis
  else
    @notimplemented "Raviart-Thomas Reference FE only available for cubes and simplices"
  end

  function cmom(φ,μ,ds) # Cell moment function: σ_K(φ,μ) = ∫(φ·μ)dK
    Broadcasting(Operation(⋅))(φ,μ)
  end
  function fmom(φ,μ,ds) # Face moment function : σ_F(φ,μ) = ∫((φ·n)*μ)dF
    n = get_facet_normal(ds)
    φn = Broadcasting(Operation(⋅))(φ,n)
    Broadcasting(Operation(*))(φn,μ)
  end

  moments = Tuple[
    (get_dimrange(p,D-1),fmom,fb), # Face moments
  ]
  if (order > 0)
    push!(moments,(get_dimrange(p,D),cmom,cb)) # Cell moments
  end

  return MomentBasedReferenceFE(RaviartThomas(),p,prebasis,moments,DivConformity())
end

function ReferenceFE(p::Polytope,::RaviartThomas,order;kwargs...)
  RaviartThomasRefFE(Float64,p,order;kwargs...)
end

function ReferenceFE(p::Polytope,::RaviartThomas,::Type{T},order;kwargs...) where T
  RaviartThomasRefFE(T,p,order;kwargs...)
end

function Conformity(reffe::GenericRefFE{RaviartThomas},sym::Symbol)
  hdiv = (:Hdiv,:HDiv)
  if sym == :L2
    L2Conformity()
  elseif sym in hdiv
    DivConformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a Raviart Thomas reference FE.

    Possible values of conformity for this reference fe are $((:L2, hdiv...)).
    """
  end
end

function get_face_own_dofs(reffe::GenericRefFE{RaviartThomas}, conf::DivConformity)
  get_face_dofs(reffe)
end

# TODO: Please remove me
function JacobiBasis(::Type{T},p::Polytope,orders) where T
  compute_jacobi_basis(T,p,orders)
end
function JacobiBasis(::Type{T},p::Polytope{D},order::Int) where {D,T}
  orders = tfill(order,Val{D}())
  JacobiBasis(T,p,orders)
end
function compute_jacobi_basis(::Type{T},p::ExtrusionPolytope{D},orders) where {D,T}
  extrusion = Tuple(p.extrusion)
  terms = _monomial_terms(extrusion,orders)
  JacobiPolynomialBasis{D}(T,orders,terms)
end
