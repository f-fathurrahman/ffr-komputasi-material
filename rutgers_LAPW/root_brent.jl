function root_brent(f, xl::Float64, xu::Float64; TOL=1.0e-10, NiterMax=100)

    a = xl
    b = xu
    fa = f(a)
    fb = f(b)

    c = a
    fc = fa
    d = b - c
    e = d

    b_old = 0.0

    for i in 1:NiterMax

        b_old = b

        if abs(fb) <= TOL
            return b
        end

        if fa*fb > 0.0
            a = c
            fa = fc
            d = b - c
            e = d
        end

        if abs(fa) < abs(fb)
            c = b
            b = a
            a = c
            #
            fc = fb
            fb = fa
            fa = fc
        end

        m = 0.5*(a - b)
        if (abs(m) <= TOL) || (abs(fb) <= TOL)
            break
        end

        if (abs(e) >= TOL) && (abs(fc) > abs(fb))
            s = fb/fc
            if abs(a - c) < eps()
                p = 2*m*s
                q = 1 - s
            else
                q = fc/fa
                r = fb/fa
                p = s*( 2*m*q*(q-r) - (b-c)*(r-1) )
                q = (q-1)*(r-1)*(s-1)
            end
            #
            if p > 0.0
                q = -q
            else
                p = -p
            end
            #
            if ( 2*p < (3*m*q - abs(TOL*q)) ) && (p < abs(0.5*e*q))
                e = d
                d = p/q
            else
                d = m
                e = m
            end
        else
            d = m
            e = m
        end
            
        c = b
        fc = fb

        if abs(d) > TOL
            b = b + d
        else
            b = b - sign(b-a)*TOL
        end
    
        fb = f(b)

        @printf("brent: %5d %18.10f %15.5e\n", i, b, abs(fb))
    end

    return b
end
