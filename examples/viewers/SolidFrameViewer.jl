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

function elipspar(p)
    elpar = (0,0,0)
    try
        elpar = ellipsoid_parameters(p)
    catch
        elpar = (0.,0.,0.)
    end

    return [round(i,2) for i in elpar]
end

function main(window)
    
    push!(window.assets,("ThreeJS","threejs"))
    push!(window.assets,"widgets")

    convert(Int,round(3.6))

    geom(frame) = geometry(ConvertToTuple(memory[frame][2],memory[frame][3])...)
    
    ff = Signal(1)
    #ThreeJS.title(2,"Viewer for time dependant calculation"),
    vbox(         
                  #"Viewer for time dependant calculation $fname",
         hbox(
            vbox(
                 #"data",
                 hbox("frame",slider(1:size(memory,1)) >>> ff),
                 ),
             ),
         vskip(2em),

         map(ff) do ff
         vbox(
              """
              Time is $(try round(memory[ff][1],2) catch NaN end);          
              Volume is $(try round(volume(memory[ff][2],memory[ff][3]),4) catch NaN end);
              Elipsoid parameters $(try [round(i,2) for i in FitEllipsoid(memory[ff][2])] catch NaN end)
             """,
          outerdiv() <<
          (
           initscene() <<
           [
         ### Pozikridis veids ka attlot meshu
            ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:color=>"red",:kind=>"normal")) ], 
            #ThreeJS.mesh(0.0, 0.0, 0.0) << [ geom(ff), material(Dict(:kind=>"basic",:color=>"black",:wireframe=>true,:wireframeLinewidth=>2)) ], # :wireframe=>true ],
            pointlight(0.,10.,0.),
            pointlight(0.,-10.,0.),
            pointlight(0.,0.,10.),
             pointlight(0.,0.,-10.),
            camera(10.0, 0., 0.),
           ]
          )
          )

         end
    )
        
end
