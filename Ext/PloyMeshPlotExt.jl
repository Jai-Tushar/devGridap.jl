module PolyMeshPlotExt

using Gridap.Arrays, Gridap.Fields
using Meshes, MeshViz
import CairoMakie as Mke

function plot_polymesh(nodes::Vector{VectorValue{2, Float64}}, connectivity::Table{Int64, Vector{Int64}, Vector{Int64}})
    elem = connect.(map(Tuple,connectivity))
    verts = map(v -> Point2(v...), nodes)
    Polymesh = SimpleMesh(verts, elem)
    viz(Polymesh, showfacets=true)
end

export plot_polymesh

end