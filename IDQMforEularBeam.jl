

using Pkg
Pkg.add("LinearAlgebra")
Pkg.add("CairoMakie")
using CairoMakie
using LinearAlgebra

struct Beam1{T}
    p0::T #kN/m #
    L::T  #m
    EI::T #kNm2
end

# Mesh scheme
function glcPointsGlobal(left::Number,right::Number,N::Int)
    glcPointsLocal= [0.5(1-cos((i-1)pi/(N-1))) for i in 1:N]
    (right-left).*glcPointsLocal.+left
end
# Calculate the DQ matrix
function CofDQM(x::Vector,order::Int)
    #x: vector of sampling points;
    #order: largest derivative order;
    #return: matrix with dimension of (x,x,m)

    n = length(x)
    C = fill(0.0,(n,n,order))
    (n>order) || error("need more sampling points")

    # 1st order
    
    for i in 1:n, j in 1:n
        if i != j
            C[i,j,1] = prod([((m!=i && m!=j) ? x[i]-x[m] : 1) for  m in 1:n])/prod([((m!=j) ? x[j]-x[m] : 1) for  m in 1:n])

            # m_range_denominater = filter!(xx->xx!=i,[m for m in 1:n])
            # m_range_numerator = filter!(xx->xx!=j,m_range_denominater)
            #C[i,j,1] = numstring(filter(xx->xx!=0,[x[i]-x[m] for m in m_range_numerator]))/numstring(filter(xx->xx!=0,[x[j]-x[m] for m in m_range_denominater]))
            # C[i,j,1] = pop!(cumprod(filter(xx->xx!=0,[x[i]-x[m] for m in m_range_numerator])))/pop!(cumprod(filter(xx->xx!=0,[x[j]-x[m] for m in 1:n])))
        else

            # C[i,i,1] = sum(filter(xx->xx!=Inf,[1/(x[i]-x[m]) for m in 1:n]))
            C[i,i,1] = sum([((m!=i) ? 1/(x[i]-x[m]) : 0) for m in 1:n]) #这个慢了4个微秒

        
        end
    end


    # higher order
    for i in 1:n, j in 1:n, k in 2:order
        if i != j
            C[i,j,k] = k*(C[i,i,k-1]*C[i,j,1] - C[i,j,k-1]/(x[i]-x[j]))
        end
    end
    for i in 1:n, j in 1:n, k in 2:order
        if i == j
            C[i,j,k] = -sum([((m!=i) ? C[i,m,k] : 0) for m in 1:n])
        end
    end
    C
end

#Plot figures, Played with mutipleDispatch
function Plot4(xGlobal,w,w_corrected,wReal)
    fig1 = Figure(backgroundcolor = RGBf(98, 98, 98),
    resolution = (1000, 300))


    Canvaskappa11 = Axis(fig1[1,1],title="Λ^{2}")
    Canvaskappa21 = Axis(fig1[2,1],title="residual")

    Canvaskappa12 = Axis(fig1[1,2],title="Λ^{1}")
    Canvaskappa22 = Axis(fig1[2,2],title="residual")



    scatter!(Canvaskappa11,[xGlobal w])
    lines!(Canvaskappa11,[xGlobal wReal.(xGlobal)],color = "black")

    scatter!(Canvaskappa21,[xGlobal wReal.(xGlobal)-w],color = "pink")
    scatter!(Canvaskappa12,[xGlobal w_corrected],color = "yellow")
    lines!(Canvaskappa12,[xGlobal wReal.(xGlobal)])

    scatter!(Canvaskappa22,[xGlobal wReal.(xGlobal)-w_corrected],color = "pink") 
    fig1
end

function Plot4(xGlobal,w,w_corrected)
    fig1 = Figure(backgroundcolor = RGBf(98, 98, 98),
    resolution = (1000, 700))


    Canvaskappa11 = Axis(fig1[1,1],title="w")
    Canvaskappa21 = Axis(fig1[2,1],title="残差")

    Canvaskappa12 = Axis(fig1[1,2],title="w_corrected")
    Canvaskappa22 = Axis(fig1[2,2],title="残差")



    # scatter!(Canvaskappa11,[xGlobal w],color = "gray")
    # lines!(Canvaskappa11,[xGlobal wReal.(xGlobal)],color = "black")

    # scatter!(Canvaskappa21,[xGlobal wReal.(xGlobal)-w],color = "pink")
    scatter!(Canvaskappa12,[xGlobal w_corrected],color = "yellow")
    # lines!(Canvaskappa12,[xGlobal wReal.(xGlobal)])

    # scatter!(Canvaskappa22,[xGlobal wReal.(xGlobal)-w_corrected],color = "yellow") 
    scatter!(Canvaskappa22,[xGlobal w-w_corrected],color = "pink") 
    fig1
end

function Plot4(xGlobal,w_corrected,wReal)
    fig1 = Figure(backgroundcolor = RGBf(98, 98, 98),
    resolution = (1000, 300))




    Canvaskappa12 = Axis(fig1[1,2],title="w_corrected")
    Canvaskappa22 = Axis(fig1[2,2],title="残差")




    scatter!(Canvaskappa12,[xGlobal w_corrected],color = "yellow")
    lines!(Canvaskappa12,[xGlobal wReal.(xGlobal)])

    scatter!(Canvaskappa22,[xGlobal wReal.(xGlobal)-w_corrected],color = "yellow") 
    fig1
end

# beam info
inputBeam = Beam1{Float64}(2,1,20)

#sampling points
N=60;


# A and Λ
xGlobal = glcPointsGlobal(0.,inputBeam.L,N)
xInner = map(x->x-inputBeam.L/2.,xGlobal) # 
DQm2 = CofDQM(xGlobal,2)[:,:,2] #A^2
DQm1 = CofDQM(xGlobal,2)[:,:,1] #A
iDQm1 = inv(DQm1)   # Λ
iDQm2 = inv(DQm2)   # Λ^2 



fig1 = Figure(backgroundcolor = RGBf(98, 98, 98),
    resolution = (2400, 700))

    labels = ["Theoretical results", "IDQM numerical results", "Residual"]
#for ∂2 B.C. hinged-hinged

if 1 == 1
    wReal(x) = x<=0.5 ? -3.125e-3x+4.1667e-3x^3 : 1.0417e-3-9.3750e-3x+1.25e-2x^2-4.1667e-3x^3
    ∂2w(x) = 0<x<=0.5inputBeam.L ? -(0.5x)/inputBeam.EI : -(0.5-0.5x)/inputBeam.EI 
    yGlobal = ∂2w.(xGlobal)
    InvInt0 = iDQm1*yGlobal # ϕ 
    InvInt1 = InvInt0 .+ (0-InvInt0[1]) # ϕ with B.C.
    InvInt00 = iDQm1*InvInt1
    InvInt11 = InvInt00 .+ [(0-InvInt00[end])+(0-InvInt00[1])]/2
    InvInt11 = -InvInt11
    Inv2Int0 = iDQm2*yGlobal
    yInv = -(Inv2Int0 .+ (0-Inv2Int0[1]) + ((0-Inv2Int0[end])-(0-Inv2Int0[1]))/inputBeam.L.*xGlobal )
    # Plot4(xGlobal, yInv,wReal)
    canticanvas11 = Axis(fig1[1,1],title="Hinged-Hinged")
    canticanvas12 = Axis(fig1[2,1],title="N = $N")
    scatter!(canticanvas11,[xGlobal yInv],color = "yellow",markersize = 8,label = labels[2])
    lines!(canticanvas11,[xGlobal wReal.(xGlobal)],color = "gray",label = labels[1])
    scatter!(canticanvas12,[xGlobal wReal.(xGlobal)-yInv],color = "pink",markersize = 6,label = labels[3]) 
    leg = Legend(fig1[1, 4],canticanvas11)
    leg = Legend(fig1[2, 4],canticanvas12)
end
#for ∂3 B.C. Cantilever beam  
if 1 == 1
    wReal(x) = -2.5e-2x^2+8.3333e-3x^3
    ∂2w(x) = x<=1 ? -(-1+1.0x)/inputBeam.EI : 0.0
    yGlobal = ∂2w.(xGlobal)
    ∂3w = DQm1*yGlobal
    ∂3w[end] = 0
    InvInt0 = iDQm1*∂3w
    # InvInt1 = InvInt0.+[(yGlobal[end]-InvInt0[end])+(yGlobal[1]-InvInt0[1])]/2
    InvInt1 = InvInt0.+(yGlobal[1]-InvInt0[1])
    InvInt00 = iDQm1*InvInt1
    InvInt11 = InvInt00.+(0-InvInt00[1])
    InvInt000 = iDQm1*InvInt11
    InvInt111 = -(InvInt000.+[0-InvInt000[1]])
    yInv = iDQm2*yGlobal
    # Plot4(xGlobal, InvInt111,wReal)
    canticanvas11 = Axis(fig1[1,2],title="Fixed-Free")
    canticanvas12 = Axis(fig1[2,2],title="N = $N")
    scatter!(canticanvas11,[xGlobal InvInt111],color = "yellow",markersize = 8)
    lines!(canticanvas11,[xGlobal wReal.(xGlobal)],color = "gray")
    scatter!(canticanvas12,[xGlobal wReal.(xGlobal)-InvInt111],color = "pink",markersize = 6) 
end

#for ∂2 B.C. Fix-hinged

if 1 == 1
    wReal(x) = -3.125e-3x^2+5.2083e-3x^3-2.0833e-3x^4
    ∂2w(x) = 0.0<=x<=inputBeam.L ? -(-0.125+0.625x-0.5x^2)/inputBeam.EI : 0.0
    yGlobal = ∂2w.(xGlobal)
    InvInt0 = iDQm1*yGlobal # ϕ 
    InvInt1 = InvInt0 .+ (0-InvInt0[1]) # ϕ with B.C.
    InvInt00 = iDQm1*InvInt1 # w_notail
    InvInt11 = InvInt00 .+ (0-InvInt00[1]) # w
    InvInt11 = -InvInt11
    yInv = -(iDQm2*yGlobal .+ (0-(iDQm2*yGlobal)[1]) + ((0-(iDQm2*yGlobal)[end])-(0-(iDQm2*yGlobal)[1]))/inputBeam.L.*xGlobal)
    FixHinge = [wReal.(xGlobal) InvInt11]
    # Plot4(xGlobal, InvInt11,wReal)
    canticanvas11 = Axis(fig1[1,3],title="Fix-Hinged")
    canticanvas12 = Axis(fig1[2,3],title="N = $N")
    scatter!(canticanvas11,[xGlobal InvInt11],color = "yellow",markersize = 8)
    lines!(canticanvas11,[xGlobal wReal.(xGlobal)],color = "gray")
    scatter!(canticanvas12,[xGlobal wReal.(xGlobal)-InvInt11],color = "pink",markersize = 6) 

    fig2 = Plot4(xGlobal,yInv,InvInt11,wReal)
end



fig1  # Figure 1 in pdf
# fig2  # Figure 2 in pdf



