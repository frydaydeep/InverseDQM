struct Beam1{T}
    p0::T #kN/m #
    L::T  #m
    EI::T #kNm2
end
include("PkgInvDQM.jl")
using CairoMakie



inputBeam = Beam1{Float64}(2,1,20)







N=60;



xGlobal = glcPointsGlobal(0.,inputBeam.L,N)

xInner = map(x->x-inputBeam.L/2.,xGlobal)

DQm2 = CofDQM(xGlobal,2)[:,:,2]
DQm1 = CofDQM(xGlobal,2)[:,:,1]
iDQm1 = inv(DQm1)
iDQm2 = inv(DQm2)

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
    InvInt00 = iDQm1*InvInt1
    InvInt11 = InvInt00 .+ (0-InvInt00[1])
    # InvInt11 = InvInt00 .+ [(0-InvInt00[end])+(0-InvInt00[1])]/2
    InvInt11 = -InvInt11
    # yInv = -(iDQm2*yGlobal.+0.5(0-(iDQm2*yGlobal)[1]+0-(iDQm2*yGlobal)[end])+((0-(iDQm2*yGlobal)[end])-(0-(iDQm2*yGlobal)[1])).*yGlobal+DQm1*(0.5(0-(iDQm2*yGlobal)[1]+0-(iDQm2*yGlobal)[end])*ones(N)))
    # yInv = -(iDQm2*yGlobal.+0.5(0-(iDQm2*yGlobal)[1]+0-(iDQm2*yGlobal)[end])+((0-(iDQm2*yGlobal)[end])-(0-(iDQm2*yGlobal)[1])).*yGlobal)
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



fig1
fig2

# Plot4(xGlobal,yGlobal,InvInt1)

