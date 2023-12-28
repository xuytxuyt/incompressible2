using ApproxOperator
include("input_gmsh.jl")

ndiv = 3
elements,nodes,elms = import_curved_beam("./msh/bar_"*string(ndiv)*".msh");