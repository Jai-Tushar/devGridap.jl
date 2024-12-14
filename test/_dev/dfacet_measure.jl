using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Helpers

using FillArrays


p = TRI
D = num_dims(p)

# T = Float64
T = VectorValue{D,Float64}


cr_reffe = CRRefFE(T,p,1)

cr_dofs = get_dof_basis(cr_reffe)
cr_prebasis  = get_prebasis(cr_reffe)
cr_shapefuns = get_shapefuns(cr_reffe) 

M = evaluate(cr_dofs, cr_shapefuns)

partition = (0,1,0,1)
cells = (2,2)
model = simplexify(CartesianDiscreteModel(partition, cells))

V = FESpace(model,cr_reffe)
get_cell_dof_ids(V)


# 
# p = TET
p = QUAD
face = 2



function get_dfacet_measure(p::Polytope{2}, face::Int) 

    measures = Float64[]

    if p == QUAD 
        perm = [1,2,4,3]
    elseif p == TRI
        perm = [1,2,3]
    end

    dim = get_dimranges(p)[face+1]

    if face == 0
        push!(measures, 0.0)
    elseif face == 1
        face_ents = get_face_coordinates(p)[dim]
        for entity in face_ents
            p1, p2 = entity
            push!(measures, norm(p2-p1))
        end
    elseif face == 2
        # Shoelace / Gauss area algo
        face_ents = get_face_coordinates(p)[dim]
        face_ents = map(Reindex(face_ents...),perm)
        shift = circshift(face_ents, -1)
    
        sum1 = map(face_ents, shift) do x1, x2
            x1[1] * x2[2]  
        end
        sum2 = map(face_ents, shift) do x1, x2
            x1[2] * x2[1]  
        end
        area = 0.5 * abs(sum(sum1)-sum(sum2))
        push!(measures, area)
    end
    return measures
end

function get_dfacet_measure(p::Polytope{3}, face::Int) 

    measures = Float64[]

    if p == HEX
        perm = [1,2,4,3]
    elseif p == TET
        perm = [1,2,3]
    end

    dim = get_dimranges(p)[face+1]

    if face == 0
        push!(measures, 0.0)
    elseif face == 1
        face_ents = get_face_coordinates(p)[dim]
        for entity in face_ents
            p1, p2 = entity
            push!(measures, norm(p2-p1))
        end
    elseif face == 2
        # # Shoelace / Gauss area algo
        # face_ents = get_face_coordinates(p)[dim]
        # face_ents = map(Reindex(face_ents...),perm)
        # shift = circshift(face_ents, -1)
    
        # sum1 = map(face_ents, shift) do x1, x2
        #     x1[1] * x2[2]  
        # end
        # sum2 = map(face_ents, shift) do x1, x2
        #     x1[2] * x2[1]  
        # end
        # area = 0.5 * abs(sum(sum1)-sum(sum2))
        # push!(measures, area)
    elseif face == 3
        @notimplemented "not implemented yet"
    end
    return measures
end

# p = TET
# dim = get_dimranges(p)[face+1]
# face_ents = get_face_coordinates(p)[dim]

# face_ents = map(face_ents) do f
#     map(Reindex(f),perm)
# end

# shift = map(face_ents) do f
#     circshift(f, -1)
# end

# sum1 = map(face_ents, shift) do face_ent, sh
#     map(face_ent, sh) do x1, x2
#         x1[1] * x2[2]  
#     end
# end

# sum2 = map(face_ents, shift) do face_ent, sh
#     map(face_ent, sh) do x1, x2
#         x1[2] * x2[1]  
#     end
# end

# facet_measure = get_dfacet_measure(p,1)