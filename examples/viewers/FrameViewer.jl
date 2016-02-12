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

    geom(frame) = geometry(ConvertToTuple(frames[frame][1],frames[frame][2])...)
    
    ff = Signal(1)
    #ThreeJS.title(2,"Viewer for time dependant calculation"),
    vbox(         
                  #"Viewer for time dependant calculation $fname",
         hbox(
            vbox(
                 #"data",
                 hbox("frame",slider(1:size(frames,1)) >>> ff),
                 ),
             ),
         vskip(2em),

         map(ff) do ff
         vbox(
             """
             Time is $(round(t[ff],2));          
             Volume is $(round(volume(frames[ff]...),4));
             Elipsoid parameters $([round(i,2) for i in ellipsoid_parameters(frames[ff][1])])
             """,
          outerdiv() <<
          (
           initscene() <<
           [
         ### Pozikridis veids ka attlot meshu
            ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:color=>"white", :kind=>"basic")) ], 
            ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:kind=>"basic",:color=>"black",:wireframe=>true,:wireframeLinewidth=>2)) ], # :wireframe=>true ],
            camera(10.0, 0., 0.),
           ]
          )
          )

         end
    )
        
end
