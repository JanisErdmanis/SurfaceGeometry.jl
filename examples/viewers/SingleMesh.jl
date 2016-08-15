using ThreeJS
import ThreeJS
using Compat

function ConvertToTuple(p,t)
    points = Array(Tuple{Float64,Float64,Float64},size(p,2))
    for i in 1:size(p,2)
        points[i] = tuple(p[:,i]...)
    end
    
    faces = Array(Tuple{Int64,Int64,Int64},size(t,2))
    for i in 1:size(t,2)
        faces[i] = tuple(t[:,i]...)
    end

    return points,faces
end


function main(window)
    
    push!(window.assets,("ThreeJS","threejs"))
    push!(window.assets,"widgets")

    convert(Int,round(3.6))

    vbox(         
         vbox(
          outerdiv() <<
          (
           initscene() <<
           [
            # For presentations
             ThreeJS.mesh(0.0, 0.0, 0.0) << [ geometry(ConvertToTuple(points,faces)...), material(Dict(:color=>"white",:kind=>"basic")) ],
             ThreeJS.mesh(0.0, 0.0, 0.0) << [ geometry(ConvertToTuple(points,faces)...), material(Dict(:color=>"black",:wireframe=>true,:kind=>"basic")) ],
            #ThreeJS.mesh(0.0, 0.0, 0.0) << [ geometry(ConvertToTuple(points,faces)...), material(Dict(:color=>"red",:kind=>"normal")) ],    
    #            pointlight(0.,10.,0.),
            camera(10.0, 0., 0.),
           ]
          )
          )
    )
end
