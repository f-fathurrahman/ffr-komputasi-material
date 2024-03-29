function calc_weight(typ, param, di, dmI)

    r  = abs(di)/dmI

    if di >= 0.0
       drdx = 1.0/dmI
    else
       drdx = -1.0/dmI
    end

    if typ == "gauss"
        w, dwdr, dwdrr = weight_gauss(param, r)
    elseif typ == "cubic"
        w, dwdr, dwdrr = weight_cubic(r)
    end

    #% EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r
    #if (type == 'GAUSS')
    #   [w,dwdr,dwdrr] = Gauss(para,r);
    #elseif (type == 'CUBIC')
    #   [w,dwdr,dwdrr] = Cubic(r);
    #elseif (type == 'SPLIN')
    #   [w,dwdr,dwdrr] = Spline(r);
    #elseif (type == 'power')
    #   [w,dwdr,dwdrr] = power_function(para,r); 
    #elseif (type == 'CSRBF')
    #   [w,dwdr,dwdrr] = CSRBF2(r);
    #else
    #   error('Invalid type of weight function.');
    #end

    dwdx  = dwdr * drdx
    dwdxx = dwdrr * drdx * drdx

    return w, dwdx, dwdxx

end

function weight_gauss(β, r)
    if r > 1.0
       w     = 0.0
       dwdr  = 0.0
       dwdrr = 0.0
    else
       b2 = β*β
       r2 = r*r
       eb2 = exp(-b2)
       w     = (exp(-b2*r2) - eb2) / (1.0 - eb2)
       dwdr  = -2*b2*r*exp(-b2*r2) / (1.0 - eb2)
       dwdrr = -2*b2*exp(-b2*r2)*(1-2*b2*r2) / (1.0 - eb2)
    end
    return w, dwdr, dwdrr
end



function weight_cubic(r)
    if r > 1.0
       w     = 0.0
       dwdr  = 0.0
       dwdrr = 0.0
    else
       w     = 1-6*r^2+8*r^3-3*r^4
       dwdr  = -12*r+24*r^2-12*r^3
       dwdrr = -12+48*r-36*r^2
    end
    return w, dwdr, dwdrr
end

#=
function [w,dwdr,dwdrr] = power_function(arfa,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0.0;
else
   a2 = arfa*arfa;
   r2 = r*r;
   w     = exp(-r2/a2);
   dwdr  = (-2*r/a2)*exp(-r2/a2);
   dwdrr  = (-2/a2+(-2*r/a2).^2)*exp(-r2/a2);
end

function [w,dwdr,dwdrr] = Spline(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0.0;
elseif (r<=0.5)
   w     = 2/3 - 4*r^2 + 4*r^3;
   dwdr  = -8*r + 12*r^2;
   dwdrr = -8 + 24*r;
else
   w     = 4/3 - 4*r + 4*r^2 - 4*r^3/3;
   dwdr  = -4 + 8*r -4*r^2;
   dwdrr = 8 - 8*r;
end

function [w,dwdr,dwdrr] = CSRBF2(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   dwdrr = 0.0;
else
   w     = (1-r)^6*(6+36*r+82*r^2+72*r^3+30*r^4+5*r^5);
   dwdr  = 11*r*(r+2)*(5*r^3+15*r^2+18*r+4)*(r-1)^5;
   dwdrr = 22*(25*r^5+100*r^4+142*r^3+68*r^2-16*r-4)*(r-1)^4;
end
=#


function MLS1d_shape( m_order, Nnodes, xi, Npoints, x, dmI, wtype, param )

    wi   = zeros(1,Nnodes)  # Weight funciton
    dwi  = zeros(1,Nnodes)
    ddwi = zeros(1,Nnodes)

    ϕ   = zeros(Npoints, Nnodes)
    dϕ  = zeros(Npoints, Nnodes)
    d2ϕ = zeros(Npoints, Nnodes)


    p   = zeros(m_order, Nnodes)
    B   = zeros(m_order, Nnodes)
    DB  = zeros(m_order, Nnodes)
    DDB = zeros(m_order, Nnodes)

    A   = zeros(m_order, m_order)
    DA  = zeros(m_order, m_order)
    DDA = zeros(m_order, m_order)

    # LOOP OVER ALL EVALUATION POINTS TO CALCULATE VALUE OF SHAPE FUNCTION Fi(X)
    for j in 1:Npoints

        # DETERMINE WEIGHT FUNCTIONS AND THEIR DERIVATIVES AT EVERY NODE
        for i in 1:Nnodes
            di = x[j] - xi[i]
            wi[i], dwi[i], ddwi[i] = calc_weight(wtype, param, di, dmI[i])
        end
   
        # EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
        if m_order == 1  # Shepard function
        
            p[1,:] = ones(1,Nnodes)
            
            px   = [1]
            dpx  = [0]
            ddpx = [0]
        
        elseif m_order == 2

            p[1,:] = ones(1,Nnodes)
            p[2,:] = xi[:]
            
            px   = [1, x[j]]
            dpx  = [0, 1]
            ddpx = [0, 0]
            
        elseif m_order == 3

            p[1,:] = ones(1,Nnodes)
            p[2,:] = xi[:]
            p[3,:] = xi.^2
            
            px   = [ 1, x[j], x[j]^2 ]
            dpx  = [ 0, 1, 2*x[j] ]
            ddpx = [ 0, 0, 2 ]
        
        else
           
           error("Invalid order of basis")
        
        end

        for i in 1:Nnodes
            for m in 1:m_order
                B[m,i] = p[m,i] * wi[i]
                DB[m,i] = p[m,i] * dwi[i]
                DDB[m,i] = p[m,i] * ddwi[i]
            end
        end

   
        # EVALUATE MATRICES A AND ITS DERIVATIVES
        # reset
        A[:,:]   .= 0.0
        DA[:,:]  .= 0.0
        DDA[:,:] .= 0.0
        for i in 1:Nnodes
            pp = p[:,i] * p[:,i]'
            A[:,:]   .= A[:,:]   .+ wi[i] * pp
            DA[:,:]  .= DA[:,:]  .+ dwi[i] * pp
            DDA[:,:] .= DDA[:,:] .+ ddwi[i] * pp
        end

        AInv = inv(A)

        rx  = AInv * px

        ϕ[j,:] = rx' * B   # shape function

        drx  = AInv * (dpx - DA * rx)
        dϕ[j,:] = drx' * B + rx' * DB
   
        ddrx  = AInv * (ddpx - 2 * DA * drx - DDA * rx)
        d2ϕ[j,:] = ddrx' * B + 2 * drx' * DB + rx' * DDB     # second order derivatives of shape function

    end

    return ϕ, dϕ, d2ϕ

end