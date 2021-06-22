
using ShunnHamQuadrature

using FastGaussQuadrature
using StaticArrays
using Combinatorics
using Plots

function _legendre(n,a,b)
    x, w = FastGaussQuadrature.gausslegendre(n)
    w .*= (b-a)/2
    x = (x.+1)/2*(b-a).+a
    collect(zip(x,w))
end

function  duffy2D( rule)
   
    qps = _legendre(rule, 0.0, 1.0)

    QL = rule

    QP = QL*QL
    uP = zeros(2,QP)
    wP = zeros(QP)
    n = 1
    uP = Vector{SVector{2,Float64}}()
  
    for (x1,w1) in qps 
        for (x2,w2) in qps
            push!(uP,[x1,x2*x1])
                         
            wP[n]   = w1*w2*x1
            n += 1
        end
    end

    return uP,wP
end

function  duffy3D(rule)
   
    qps = _legendre(rule, 0.0, 1.0)

    QL = rule

    QP = QL*QL*QL
  
    wP = zeros(QP)
    n = 1
    uP = Vector{SVector{3,Float64}}()
  
    for (x1,w1) in qps 
        for (x2,w2) in qps
            for (x3,w3) in qps
                push!(uP,[x1,x2*x1,x3*x2*x1])
                            
                wP[n]   = w1*w2*w3*x1^2*x2
                n += 1
            end
        end
    end

    return uP,wP
end


function  duffy4D(rule)
   
    qps = _legendre(rule, 0.0, 1.0)

    QL = rule

    QP = QL*QL*QL*QL
  
    wP = zeros(QP)
    n = 1
    uP = Vector{SVector{4,Float64}}()
  
    for (x1,w1) in qps 
        for (x2,w2) in qps
            for (x3,w3) in qps
                for (x4,w4) in qps
                    push!(uP,[x1,x2*x1,x3*x2*x1,x4*x3*x2*x1])
                                
                    wP[n]   = w1*w2*w3*w4*x1^3*x2^2*x3
                    n += 1
                end
            end
        end
    end

    return uP,wP
end

function  duffy5D(rule)
   
    qps = _legendre(rule, 0.0, 1.0)

    QL = rule

    QP = QL*QL*QL*QL*QL
  
    wP = zeros(QP)
    n = 1
    uP = Vector{SVector{5,Float64}}()
  
    for (x1,w1) in qps 
        for (x2,w2) in qps
            for (x3,w3) in qps
                for (x4,w4) in qps
                    for (x5,w5) in qps
                        push!(uP,[x1,x2*x1,x3*x2*x1,x4*x3*x2*x1,x5*x4*x3*x2*x1])
                                    
                        wP[n]   = w1*w2*w3*w4*w5*x1^4*x2^3*x3^2*x4
                        n += 1
                    end
                end
            end
        end
    end

    return uP,wP
end





function test_one()
    return 1
end

function test_poly_3d(x,m)
    sum = 0.0
    denom = 0.0
    for i in 0:m
        for j in 0:m-i
            for k in 0:m-i-j
                sum += (i+1)*(j+1)*(k+1)*x[1]^i*x[2]^j*x[3]^k
            end
        end
        for j in 0:i+1
            denom += j
        end
    end
    return sum/denom
end

#evaluate polynom
function eval_poly(x,e)
    p = 1.0
   
    for (i,y) in enumerate(x)
        
        p*=y^e[i]
    end
    return p
end

function test_poly_2d(x,m)
    sum = 0.0
    denom = (m+1)*(m+2)
    for i in 0:m
        for j in 0:m-i
            sum += 2*(i+1)*(j+1)*x[1]^i*x[2]^j
        end
    end
    return sum/denom
end

#Test Order
exact=0.5

#vertices = [[0 0], [1 0], [1 1]]
#vertices = [[0 0 0], [1 0 0], [1 1 0], [1 1 1]]
#vertices = [[0 0 0 0], [1 0 0 0], [1 1 0 0], [1 1 1 0], [1 1 1 1]]
vertices = [[0 0 0 0 0], [1 0 0 0 0], [1 1 0 0 0], [1 1 1 0 0], [1 1 1 1 0], [1 1 1 1 1]]

(x,w) = duffy5D(12)
qps_high = zip(x,w)
println(sum(w))

order = []
err_ref =[]
err_new =[]
for k in 5:5

    

    #sh_ref = shunnham5D_ref(k)
    #(x,w) = get_pts_wts(sh_ref,vertices)
    #qps_ref = zip(x,w)

    #println(sum(w))
    

    #sh_new = shunnham3D_alt(k)
    #(x,w) = get_pts_wts(sh_new,vertices)
    #qps_new = zip(x,w)

    (x,w) = duffy5D(k)
    qps_new = zip(x,w)

    println(sum(w))

    for i in 0:8

        exponents = collect(multiexponents(5,i))
        println(length(exponents))
        for e in exponents
            function integrand(x)
                eval_poly(x,e)
            end

            #Refrence
            exact  =  sum(w*integrand(x) for (x,w) in qps_high)

            #int_ref = sum(w*integrand(x) for (x,w) in qps_ref)
            int_new = sum(w*integrand(x) for (x,w) in qps_new)

            println(i," ",e)
            println(exact)
            #println(int_ref)
            println(int_new)

            push!(order,i)
            #push!(err_ref,abs(int_ref-exact))
            push!(err_new,abs(int_new-exact))
        end
      
    end

end


plot(size=(800,500),yaxis=:log)


#scatter!(order,err_ref,label="Ref",markershape=:x)
scatter!(order,err_new,label="New",markershape=:o)


plot!(xlims=(-1,13),xticks=0:1:12,ylims=(1e-17,10),ylabel="Max. rel. error",xlabel="Max. degree of polynom",legend=:outerright, title="Order of Shunn-Ham Quadrature")

