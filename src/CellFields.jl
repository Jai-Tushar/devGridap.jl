export CellField
export evaluate, compose

"""
Abstract type that represents a cell-wise field, where
T stands for the type that represents the field at a point
(e.g., scalar, vector, tensor) and D stands for the space
dimension
"""
const CellField{D,T} = EvaluableCellArray{D,T,1} where {D,T}

evaluate(::CellField{D,T} where {D,T},::CellPoints{D} where D)::CellFieldValues{T} = @abstractmethod

"""
Returns another CellField object that represents the gradient.
TG has a rank one order greater than the one of T
"""
gradient(::CellField{D,T} where {D,T})::CellField{D,TG} = @abstractmethod

"""
Composes a lambda function with a `CellField`
`f` has to be overloaded with 2 methods one that returns the type of the
result, and another that returns the result
"""
function compose(f,g::CellField{D,S}) where {D,S}
  T = f(S)
  CellFieldFromComposeWithLambda{D,S,T}(f,g)
end

