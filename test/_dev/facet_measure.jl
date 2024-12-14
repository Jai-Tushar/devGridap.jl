using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Helpers

using FillArrays
using LinearAlgebra


p = TRI
dim = get_dimranges(p)[3]
perm = [1,2,3]
face_ents = get_face_coordinates(p)[dim]
face_ents = map(Reindex(face_ents...),perm)

# T = Float64
T = VectorValue{D,Float64}


# cr_reffe = CRRefFE(T,p,1)

# cr_dofs = get_dof_basis(cr_reffe)
# cr_prebasis  = get_prebasis(cr_reffe)
# cr_shapefuns = get_shapefuns(cr_reffe) 

# M = evaluate(cr_dofs, cr_shapefuns)

partition = (0,1,0,1)
cells = (2,2)
model = simplexify(CartesianDiscreteModel(partition, cells))

V = FESpace(model,cr_reffe)
get_cell_dof_ids(V)


# 
# p = TET
p = QUAD
face = 2



function get_facet_measure(p::Polytope{2}, face::Int) 

    measures = Float64[]

    if p == QUAD 
        perm = [1,2,4,3]
    elseif p == TRI
        perm = [1,2,3]
    end

    dim = get_dimranges(p)[face+1]

    if face == 0
        face_ents = get_face_coordinates(p)[dim]
        for entity in face_ents
            push!(measures, 0.0)
        end
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

function get_face_measure(p::Polytope{3}, face::Int) 

    measures = Float64[]

    dim = get_dimranges(p)[face+1]

    if face == 0
        face_ents = get_face_coordinates(p)[dim]
        for entity in face_ents
            push!(measures, 0.0)
        end
    elseif face == 1
        face_ents = get_face_coordinates(p)[dim]
        for entity in face_ents
            p1, p2 = entity
            push!(measures, norm(p2-p1))
        end
    elseif face == 2
        if p == HEX 
            perm = [1,2,4,3]
        elseif p == TET
            perm = [1,2,3]
        end

        @notimplemented "not implemented yet"

    elseif face == 3
        @notimplemented "not implemented yet"
    end
    return measures
end

# p = TET
# dim = get_dimranges(p)[2+1]
# face_ents = get_face_coordinates(p)[dim]
# # perm = [1,2,4,3]
# perm = [1,2,3]
# face_ents = map(face_ents) do x
#     map(Reindex(x),perm)
# end

# area_vec = []
# for i in 1:length(face_ents)
#     push!(area_vec,[])
#     for j in 1:length(face_ents[i])
#        v1 = face_ents[i][j]
#        v2 = face_ents[i][mod(j,length(face_ents[i]))+1]
#        vec1 = v1 - face_ents[i][1]
#        push!(area_vec[i],cross(v1,v2))
#     end
# end

# 0.5 * norm(sum(area_vec[2]))

# v1 = face_ents[1][1]
# v2 = face_ents[1][2]
# push!(area_vec,cross(v1,v2))