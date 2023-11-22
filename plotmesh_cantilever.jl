using ApproxOperator,CairoMakie

lwb = 2.0;lwm = 0.5;mso = 7;msx = 7;ppu = 2.5;α = 0.7;
filename1 = "./msh/cantilever_20.msh"
filename2 = "./msh/cantilever_bubble_1600.msh"
savename = "./png/1.png"
elms,~= ApproxOperator.importmsh(filename1)
elms_p,~ = ApproxOperator.importmsh(filename2)

x = elms["Ω"][1].x
y = elms["Ω"][1].y
z = elms["Ω"][1].z
xᵖ = elms_p["Ω"][1].x
yᵖ = elms_p["Ω"][1].y
zᵖ = elms_p["Ω"][1].z

if occursin("quad",filename1)
    index = [1,2,3,4,1]
else
    index = [1,2,3,1]
end

f = Figure(backgroundcolor = :transparent)
ax = Axis(f[1,1],aspect = DataAspect(),backgroundcolor = :transparent)
hidespines!(ax)
hidedecorations!(ax)
L = 48.
b = 12.
lines!([0.0,L,L,0.0,0.0],[-b/2,-b/2,b/2,b/2,-b/2], linewidth = lwb, color = :black)

for elm in elms["Ω"]
    id = [i for i in elm.i]
    lines!(x[id[index]],y[id[index]], linewidth = lwm, color = :black)
end
scatter!(x,y,marker = :circle, markersize = mso, color = :black)

# for elm in elms_p["Ω"]
#     id = [i for i in elm.i]
#     lines!(xᵖ[id[[1,2,3,1]]],yᵖ[id[[1,2,3,1]]], linewidth = lwm, color = :blue)
# end
scatter!(xᵖ,yᵖ,marker = :xcross, markersize = msx, color = (:blue, α))
save(savename,f,px_per_unit = ppu)
f