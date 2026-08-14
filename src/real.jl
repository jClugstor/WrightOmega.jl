# The Wright omega function, omega(x): the unique real solution of
#
#     omega + ln(omega) - x = 0
#
# It is the same object as W0(exp(x)), but evaluated WITHOUT ever forming exp(x),
# which is the entire reason it is worth having: the bypass form W0(exp(x))
# overflows for x > 709.78 and underflows for x < -745.1, and the failure is
# permanent because Inf and NaN do not recover inside a Newton iteration.
#
# Method: Fukushima (2020), "Fast computation of Wright omega function by piecewise
# minimax rational function approximation", J. Comput. Appl. Math.  A conditional
# switch between
#
#   (i)   a rational function of z = exp(x) for x below x_0, where W0 is analytic
#         at z = 0 and so is efficiently approximated in z rather than x;
#   (ii)  twelve piecewise minimax rational functions R_1..R_12 in x;
#   (iii) L(x) = x*(1 - ln(x)/(x+1)) for very large x, which is one Newton
#         correction applied to the asymptote omega(x) -> x.
#
# Coefficients are Fukushima's Tables 2 and 11-23 (double precision)


# ---------------------------------------------------------------- region bounds
# Fukushima Table 2. R_k is used on [X[k-1], X[k]); R0(exp(x)) below X0; L(x) above X12.
const WO_X0  = -1.7969970703125
const WO_X1  =  0.0
const WO_X2  =  2.670658111572265625
const WO_X3  =  9.6968069076538085938
const WO_X4  =  3.8007662773132324219e+1
const WO_X5  =  2.2405800342559814453e+2
const WO_X6  =  1.3188293657302856445e+3
const WO_X7  =  9.9018982133865356445e+3
const WO_X8  =  1.0472998977375030518e+5
const WO_X9  =  1.7873714379758834839e+6
const WO_X10 =  6.4666226364722251892e+7
const WO_X11 =  9.6537874592561578751e+9
const WO_X12 =  3.7284375517631785583e+11

# ---------------------------------------------------------------- coefficients
# Each pair is (numerator, denominator); the denominator's leading 1 is explicit so
# that `evalpoly` can be used directly on both.

# Table 11: R0(z), type (5,5), argument z = exp(x)
const WO_P0 = (0.0,
               9.9999999999999989310e-1,
               5.0745983953466943659,
               8.2606009602900812771,
               4.7024717601199799556,
               6.5482227387481846600e-1)
const WO_Q0 = (1.0,
               6.0745983953465398344,
               1.2835199355673150217e+1,
               1.1092440186086891869e+1,
               3.4850593056113482442,
               2.3505959341677892216e-1)

# Table 12: R1(x), type (7,7)
const WO_P1 = (5.6714329040978393597e-1,
               4.0127273991556628416e-1,
               1.6948622586211161962e-1,
               4.6665932688016784538e-2,
               8.5703273322828587102e-3,
               1.0434401063688186806e-3,
               7.7800083733980925540e-5,
               2.7204772831680988214e-6)
const WO_Q1 = (1.0,
               6.9429514456939747074e-2,
               1.2462836200029189075e-1,
              -3.8952406229635946107e-3,
               4.4555853359591692983e-3,
              -4.1123146046694083146e-4,
               5.9020981275047793275e-5,
              -3.5492274716282929464e-6)

# Table 13: R2(x), type (7,7)
const WO_P2 = (5.6714329040978381003e-1,
               7.7356892277028762975e-1,
               5.1832380984628428402e-1,
               2.1755702154802737458e-1,
               6.1224277815768787167e-2,
               1.1518883525838608343e-2,
               1.3449706921997516673e-3,
               7.5896892739201294405e-5)
const WO_Q2 = (1.0,
               7.2587064520843295454e-1,
               3.2082921273199275229e-1,
               8.6948809834337515365e-2,
               1.5394071575832887833e-2,
               1.6360582874079300062e-3,
               7.8307305162627353776e-5,
              -3.1442237428129457031e-8)

# Table 14: R3(x), type (7,7)
const WO_P3 = (5.6713501343190366161e-1,
               7.8382455352242915730e-1,
               5.1756353812338515893e-1,
               2.0726853039949915898e-1,
               5.3010998762456477240e-2,
               8.3925496655250859740e-3,
               6.6092716371736707786e-4,
               1.3614908670113335271e-5)
const WO_Q3 = (1.0,
               7.4389593435342377573e-1,
               3.0809335954228387897e-1,
               7.4470884949693883215e-2,
               1.0667765269781492069e-2,
               7.2813578520727560898e-4,
               1.3735475456227982957e-5,
              -3.9617637001432374240e-10)

# Table 15: R4(x), type (7,7)
const WO_P4 = (5.6401353907858227635e-1,
               7.1556555400973446955e-1,
               4.0606219668682194772e-1,
               1.3098742100043090882e-1,
               2.3792021846516585517e-2,
               1.4712087559790111411e-3,
               2.8262095549717529699e-5,
               1.3176541606714501011e-7)
const WO_Q4 = (1.0,
               6.1723923788369990709e-1,
               1.9472671474858151447e-1,
               2.9166374057621022340e-2,
               1.6049364440594188061e-3,
               2.9097481892430916950e-5,
               1.3204819260704942813e-7,
              -2.2759460679274951500e-13)

# Table 16: R5(x), type (8,6) -- note the asymmetric degrees
const WO_P5 = (7.4904304655011303763e-1,
               4.1475765300820204606e-1,
               2.6229090010420327086e-1,
               2.2231827737742717397e-2,
               5.5769724215842159151e-4,
               4.8929940562358938904e-6,
               1.4576036364817152103e-8,
               1.1211845148576067899e-11,
               5.8617957261146552764e-19)
const WO_Q5 = (1.0,
               3.3914957488601671439e-1,
               2.4691496687221735503e-2,
               5.8391033761937493228e-4,
               4.9870990455400355673e-6,
               1.4666816425760107778e-8,
               1.1216014296054046261e-11)

# Table 17: R6(x), type (7,7)
const WO_P6 = (-5.2847240530691242510e-1,
                6.6439555408432191153e-1,
                9.5734165623936586970e-2,
                1.7984402717797458483e-3,
                9.5608645967684762747e-6,
                1.7130211053307026173e-8,
                1.0127406535507852349e-11,
                1.5119473539633986399e-15)
const WO_Q6 = (1.0,
               1.0436793508415000619e-1,
               1.8544824235611654910e-3,
               9.6789573725809544085e-6,
               1.7211790608220386441e-8,
               1.0142134330206435119e-11,
               1.5120546649245173808e-15,
              -2.8371714954759691899e-24)

# Table 18: R7(x), type (7,7)
const WO_P7 = (-2.7265755587675828255,
                9.1886110507404644509e-1,
                1.5437556964724977333e-2,
                4.2379757069402316910e-5,
                3.3788408654334341745e-8,
                9.0995247456328828154e-12,
                8.0541623626962845941e-16,
                1.7892200405606258679e-20)
const WO_Q7 = (1.0,
               1.5717280527236316462e-2,
               4.2640053213855619617e-5,
               3.3868168872656516608e-8,
               9.1075366623274182120e-12,
               8.0562465016652583741e-16,
               1.7892388148342752402e-20,
              -7.2802540742878224330e-31)

# Table 19: R8(x), type (7,7)
const WO_P8 = (-4.9758448332804315705,
                9.8672919773901375858e-1,
                1.8048502584507968137e-3,
                5.6727618232662977651e-7,
                5.1693477737117469717e-11,
                1.5763106040543895563e-15,
                1.5615457093141147610e-20,
                3.8348282655423816844e-26)
const WO_Q8 = (1.0,
               1.8098048079446830033e-3,
               5.6778598593151384010e-7,
               5.1710727052210149252e-11,
               1.5765001969393981202e-15,
               1.5615988471610678184e-20,
               3.8348326592462094894e-26,
              -1.8376700258834952892e-38)

# Table 20: R9(x), type (7,7)
const WO_P9 = (-7.5064713916128040209,
                9.9856667809261792330e-1,
                1.4527404521491427728e-4,
                3.6058443746980183631e-9,
                2.5487488372192215976e-14,
                5.9071573396622281840e-20,
                4.3540463254359589235e-26,
                7.7870710830579392924e-33)
const WO_Q9 = (1.0,
               1.4531469783952414582e-4,
               3.6061610745034193796e-9,
               2.5488287561466743981e-14,
               5.9072215616613797183e-20,
               4.3540591638996294464e-26,
               7.7870717193388587180e-33,
              -1.8603410329812104693e-47)

# Table 21: R10(x), type (7,7)
const WO_P10 = (-1.0553133866290087103e+1,
                 9.9991182205756461588e-1,
                 6.8147130642100950299e-6,
                 7.6772315747426585717e-12,
                 2.3749238044176130427e-18,
                 2.3156775350327431915e-25,
                 6.8906814806601365039e-33,
                 4.7733356492001687570e-41)
const WO_Q10 = (1.0,
                6.8148234879449250418e-6,
                7.6772685608307803918e-12,
                2.3749276757415131404e-18,
                2.3156787740031116785e-25,
                6.8906824237918392311e-33,
                4.7733356637141783505e-41,
               -1.5212714240163636517e-58)

# Table 22: R11(x), type (7,7)
const WO_P11 = (-1.4453233847636304145e+1,
                 9.9999773787655695790e-1,
                 1.3384094485652721805e-7,
                 2.7649127536426257913e-15,
                 1.4462424184214411200e-23,
                 2.1739442557543959836e-32,
                 9.0218885020989014506e-42,
                 7.8687375038593681713e-52)
const WO_Q11 = (1.0,
                1.3384099578765282610e-7,
                2.7649130385236672771e-15,
                1.4462424639425969762e-23,
                2.1739442758770299520e-32,
                9.0218885211281218582e-42,
                7.8687375041325222089e-52,
               -2.9925825575484892462e-73)

# Table 23: R12(x), type (5,4)
const WO_P12 = (-1.9911981551589239565e+1,
                 9.9999998802702964239e-1,
                 5.3193290756933022976e-10,
                 3.6315161298281125840e-20,
                 4.6372938690029442720e-31,
                 9.2610561378209620950e-43)
const WO_Q12 = (1.0,
                5.3193290845021667622e-10,
                3.6315161310308030589e-20,
                4.6372938692633433675e-31,
                9.2610561378259071409e-43)

@inline _wo_rat(x, P, Q) = evalpoly(x, P) / evalpoly(x, Q)

"""
    wrightomega(x::Real) -> Real

The Wright omega function: the unique real solution `w > 0` of `w + log(w) == x`.

Equal to `lambertw(exp(x))` but evaluated without forming `exp(x)`, so it is finite and
accurate over the whole floating-point range of `x`. Accurate to a few ulp and costs
less than one `exp` call. Results below `floatmin(Float64)` flush to zero, so
`wrightomega` returns exactly `0.0` for `x < -708.4`.

```
julia> wrightomega(1.0)            # 1 + log(1) == 1
1.0

julia> wrightomega(0.0)            # the omega constant
0.5671432904097839
```

Method and coefficients: Fukushima (2020). Complex arguments use TOMS Algorithm 917.
"""
@inline function wrightomega(x::Float64)
    # Non-finite arguments are screened OUTSIDE the @fastmath kernel. @fastmath
    # implies `ninf`/`nnan`, so a guard placed inside it is optimised away -- which
    # is why the reference implementation returns NaN for +Inf.
    isfinite(x) || return x > 0 ? Inf : (x < 0 ? 0.0 : x)
    return _wrightomega(x)
end

@inline @fastmath function _wrightomega(x::Float64)
    if x < WO_X1                        # x < 0
        if x < WO_X0
            # W0 is analytic at z = 0, so approximate in z = exp(x). Underflow of
            # exp is benign here: z -> 0 gives R0 -> 0, which is the correct limit.
            return _wo_rat(exp(x), WO_P0, WO_Q0)
        else
            return _wo_rat(x, WO_P1, WO_Q1)
        end
    elseif x < WO_X4                    # 0 <= x < 38.008
        if x < WO_X2
            return _wo_rat(x, WO_P2, WO_Q2)
        elseif x < WO_X3
            return _wo_rat(x, WO_P3, WO_Q3)
        else
            return _wo_rat(x, WO_P4, WO_Q4)
        end
    elseif x < WO_X8                    # 38.008 <= x < 1.047e5
        if x < WO_X5
            return _wo_rat(x, WO_P5, WO_Q5)
        elseif x < WO_X6
            return _wo_rat(x, WO_P6, WO_Q6)
        elseif x < WO_X7
            return _wo_rat(x, WO_P7, WO_Q7)
        else
            return _wo_rat(x, WO_P8, WO_Q8)
        end
    elseif x < WO_X12                   # 1.047e5 <= x < 3.728e11
        if x < WO_X9
            return _wo_rat(x, WO_P9, WO_Q9)
        elseif x < WO_X10
            return _wo_rat(x, WO_P10, WO_Q10)
        elseif x < WO_X11
            return _wo_rat(x, WO_P11, WO_Q11)
        else
            return _wo_rat(x, WO_P12, WO_Q12)
        end
    else
        # One Newton correction to the asymptote omega -> x. Fukushima eq. (6).
        # Already below double-precision epsilon for x > 9.484e5, well under X12.
        return x * (1.0 - log(x) / (x + 1.0))
    end
end

# Every other real type round-trips through the Float64 kernel and converts back to its
# own float type (Float32 -> Float32, Integer/Rational -> Float64, ...). Big types get
# Float64 accuracy. `oftype(float(x), ...)` rather than `wrightomega(float(x))` so that
# types with `float(T) === T` (e.g. BigFloat) cannot recurse.
wrightomega(x::T) where {T<:Real} = oftype(float(x), wrightomega(Float64(x)))

wrightomega(::Missing) = missing
