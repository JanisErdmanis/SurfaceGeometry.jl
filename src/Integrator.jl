### Kautkas līdzīgs jāpamēģina uzrakstīt katrai virsotnei
abstract TimeStepers

type EilerStep <: TimeStepers
    t
    x
end

function step!(s::EilerStep,velocity,h)
    s.x += h*velocity
    s.t += h
end


type AdamsStep <: TimeStepers
    t::Float64                  # This variable will be shared, we can make it asan age variable
    x
    lasth::Float64
    lastv
    stepcount::Int

    function AdamsStep(t,x)
        new(t,x,0,similar(x),0)
    end
end

function step!(s::AdamsStep,velocity,h)
    if s.stepcount==0
        s.x += velocity*h
    else
        #y2 = s.y + 0.5*h/s.h*((2*s.h + h)*velocity - h*s.velocity)
        s.x += 0.5*h/s.lasth*( (2*s.lasth + h)*velocity - h*s.lastv)
    end
    s.t += h
    s.lastv = velocity
    s.stepcount += 1
    s.lasth = h
end

function step!{T <: TimeStepers}(steppers::Array{T,1},velocity::Array{Float64,2},h::Float64)
    points = Array(Float64,3,size(steppers)...)
    for xkey in 1:size(steppers,1)
        #s = steppers[xkey]
        step!(steppers[xkey],velocity[:,xkey],h)
        points[:,xkey] = steppers[xkey].x[:]
    end

    return points
end


### This should be usefull with use of linspace and some adaptive stepping technique
function multistep!(s,f,xf)

    stepsize(x,y,velocity) = 0.01
    
    iterations = 0
    while true
        iterations+=1
        velocity = f(s.x,s.y)
        h = stepsize(s.x,s.y,velocity)
        if (s.x+h)<xf
            step!(s,velocity,h)
            #println("x=$(s.x) y=$(s.y)")
        else
            h = xf-s.x
            step!(s,velocity,h)
            #println("x=$(s.x) y=$(s.y)")
            break
        end
    end
    return iterations
end


### Estimation of proper step
function properstep(cmsh,velocity)
    points = cmsh.points
    time = Inf
    
    for vi in 1:size(cmsh.points,2)
        face = cmsh.vfaces[vi]
        tri = cmsh.faces[:,face]
        w = tri .== vi
        
        p0 = points[:,tri[w]...]
        p1 = points[:,tri[w[[3,1,2]]]...]

        veloc = norm(velocity[:,vi])
        timei, = abs(norm(p0-p1)/veloc) 
        if timei<time
            time = timei
        end
    end
    return time
end
