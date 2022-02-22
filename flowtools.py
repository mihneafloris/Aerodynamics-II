def flowisentropic2(GAMMA, VAR, MTYPE):
    #   FLOWISENTROPIC2(GAMMA, VAR, MTYPE) returns the
    #   isentropic flow Mach number, MACH, temperature ratio, T,
    #   pressure ratio, P, density ratio, RHO, and area ratio, A.
    #   FLOWISENTROPIC2 calculates these arrays given the specific heat
    #   ratio, GAMMA, and any one of the isentropic flow variables.  The
    #   isentropic flow variable is selected by the string, MTYPE. The
    #   temperature, pressure, and density ratios are a comparison of the local
    #   static conditions over the stagnation (or total) conditions.  The area
    #   ratio is a comparison of the instantaneous streamtube area for a
    #   quasi-one-dimensional flow over the reference area (throat area or
    #   minimum area) where the Mach number of the flow becomes unity.
    #
    #   Inputs for FLOWISENTROPIC2 are:
    #
    #   GAMMA :  Specific heat ratio.
    #
    #   VAR   :  Numerical values for one of the isentropic flow relations.
    #
    #            MACH :  Mach number. MACH is used with MTYPE variable 'mach'. 
    #
    #            T    :  Temperature ratio.  The temperature ratio is defined 
    #                    as the local static temperature over the stagnation
    #                    temperature. T is used with MTYPE variable 'temp'.
    #
    #            P    :  Pressure ratio.  The pressure ratio is defined as the
    #                    local static pressure over the stagnation pressure. P 
    #                    is used with MTYPE variable 'pres'.
    #
    #            RHO  :  Density ratio.  The density ratio is defined as the 
    #                    local density over the stagnation density.  RHO is 
    #                    used with MTYPE variable 'dens'.
    #
    #            A    :  Area ratio.  A is used with MTYPE variables 'sub' or
    #                    'sup'.
    #
    #   MTYPE :  A string for selecting the isentropic flow variable
    #             represented by VAR.
    #
    #             'mach'  :  Default value.  Indicates that the function is in
    #                        Mach number input mode.
    #   
    #             'temp'  :  Indicates that the function is in temperature
    #                        ratio input mode.
    #
    #             'pres'  :  Indicates that the function is in pressure ratio
    #                        input mode.
    #
    #             'dens'  :  Indicates that the function is in density ratio
    #                        input mode.
    #
    #             'sub'   :  Indicates that the function is in subsonic area
    #                        ratio input mode.  The subsonic area ratio is
    #                        defined as the local subsonic streamtube area over
    #                        the reference streamtube area for sonic
    #                        conditions.
    #
    #             'sup'   :  Indicates that the function is in supersonic area
    #                        ratio input mode.  The supersonic area ratio is
    #                        defined as the local supersonic streamtube area
    #                        over the reference streamtube area for sonic
    #                        conditions.
    #
    #   Outputs calculated for FLOWISENTROPIC2 are a tuple arranged as:
    #   (MACH, T, P, RHO, A)
    #
    #   MACH    :  Mach number.
    #
    #   T       :  Temperature ratio.  The temperature ratio is defined as the 
    #              local static temperature over the stagnation temperature.
    #
    #   P       :  Pressure ratio.  The pressure ratio is defined as the local 
    #              static pressure over the stagnation pressure.
    #
    #   RHO     :  Density ratio.  The density ratio is defined as the local 
    #              density over the stagnation density.
    #
    #   A       :  Area ratio.  The area ratio is defined as the local 
    #              streamtube area over the reference streamtube area for
    #              sonic conditions.
    #
    # Kyle Lynch - k.p.lynch at tudelft.nl

    from math import sqrt

    # Parse arguments and calculate the Mach number.
    if MTYPE == 'mach':
        M = VAR;
    elif MTYPE == 'sub':
        # Subsonic area ratio solution.
        # http://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
        P = 2/(GAMMA+1)
        Q = 1-P
        R = VAR**2
        a = P**(1/Q)
        r = (R-1)/(2*a)
        X = 1/( (1+r) + sqrt(abs(r*(r+2))) )
        for i in range(0,3,1):
            X = P*(X-1)/(1-R*(P+Q*X)**(-P/Q))
        M = sqrt(X)
    elif MTYPE == 'sup':
        # Supersonic area ratio solution.
        # http://www.grc.nasa.gov/WWW/winddocs/utilities/b4wind_guide/mach.html
        P = 2/(GAMMA+1)
        Q = 1-P
        R = VAR**(2*Q/P)
        a = Q**(1/P)
        r = (R-1)/(2*a)
        X = 1/( (1+r) + sqrt(r*(r+2)) )
        for i in range(0,3,1):
            X = Q*(X-1)/(1-R*(Q+P*X)**(-Q/P))
        M = 1/sqrt(X)
    elif MTYPE == 'pres':
        # Pressure ratio solution.
        M = sqrt((2/(GAMMA-1))*(VAR**(-(GAMMA-1)/GAMMA)-1))
    elif MTYPE == 'dens':
        # Density ratio solution.
        M = sqrt((2/(GAMMA-1))*(VAR**(-(GAMMA-1)/1)-1))
    elif MTYPE == 'temp':
        # Temperature ratio solution.
        M = sqrt((2/(GAMMA-1))*(VAR**(-1)-1))
    else:
        print('Unsupported mode.')

    # Organize outputs.
    MACH = M;
    T = 1 / (1 + 0.5*(GAMMA-1)*M**2);
    P = 1 / ((1 + 0.5*(GAMMA-1)*M**2)**(GAMMA/(GAMMA-1)));
    RHO = 1 / ((1 + 0.5*(GAMMA-1)*M**2)**(1/(GAMMA-1)));
    A = sqrt( (1/M**2)*((2/(GAMMA+1))*(1+0.5*(GAMMA-1)*M**2))**((GAMMA+1)/(GAMMA-1)) );

    return (MACH, T, P, RHO, A)

def flownormalshock2(GAMMA, VAR, MTYPE):
    #   FLOWNORMALSHOCK2(GAMMA, VAR, MTYPE, MACH, T, P, RHO, M, P0, P1)
    #   returns the normal shock relations for a given specific heat ratio, 
    #   GAMMA, and any one of the normal shock relations.  The normal shock 
    #   relations are selected by the string, MTYPE.  The outputs are the Mach
    #   number, MACH, temperature ratio, T, the (static) pressure ratio, P, the
    #   density ratio, RHO, downstream Mach number, M, and the total
    #   (stagnation) pressure ratio, P0.  All of these ratios are downstream
    #   value over upstream value.  P1 is the Rayleigh-Pitot ratio, which is
    #   the ratio of upstream static pressure over the downstream stagnation
    #   pressure.  Note that upstream is said to be "before" or "ahead" of the
    #   shock.  Downstream is said to be "after" or "behind" the shock.

    #   Inputs for FLOWNORMALSHOCK2 are:
    #
    #   GAMMA :  Specific heat ratio.
    #
    #   VAR   :  Numerical values for one of the isentropic flow relations.
    #
    #            MACH :  Mach number. MACH is used with MTYPE variable 'mach'. 
    #
    #   MTYPE  :  A string for selecting the isentropic flow variable
    #             represented by VAR.
    #
    #             'mach'  :  Default value.  Indicates that the function is in
    #                        Mach number input mode.
    #
    #   Outputs calculated for FLOWNORMALSHOCK are given as a tuple in the format
    #   (MACH, T, P, RHO, M, P0, P1)
    #
    #   MACH    :  Upstream Mach number.
    #
    #   T       :  Temperature ratio.  The temperature ratio is defined as the
    #              static temperature downstream of the shock over the static 
    #              temperature upstream of the shock.
    #
    #   P       :  Pressure ratios.  The pressure ratio is defined as the 
    #              static pressure downstream of the shock over the static
    #              pressure upstream of the shock.
    #
    #   RHO     :  Density ratio.  The density ratio is defined as the density
    #              of the fluid downstream of the shock over the density 
    #              upstream of the shock.
    #
    #   M       :  Downstream Mach number.
    #
    #   P0      :  Total pressure ratio.  The total pressure ratio is defined
    #              as the total pressure downstream of the shock over the total
    #              pressure upstream of the shock.
    #
    #   P1      :  Rayleigh-Pitot ratio.  The Rayleigh-Pitot ratio is defined 
    #              as the static pressure upstream of the shock over the total 
    #              pressure downstream of the shock.
    #
    # Kyle Lynch - k.p.lynch at tudelft.nl

    from math import sqrt, log, exp

    # Parse arguments and calculate the upstream Mach number.
    if MTYPE == 'mach':
        MACH = VAR;
    else:
        print('Unsupported mode.')

    # Make sure upstream Mach number at least 1.
    if (MACH < 1):
        print('Upstream Mach number less than 1.')

    # Constants.
    cp = 1.01 * 1000;
    R = 287;

    # Organize outputs.
    t1 = (1+(2*GAMMA)/(GAMMA+1)*(MACH**2-1))
    t2 = (2+(GAMMA-1)*MACH**2) / ((GAMMA+1)*MACH**2)
    T = t1*t2

    P = 1 + ((2*GAMMA)/(GAMMA+1))*(MACH**2-1)
    RHO = ((GAMMA+1)*MACH**2)/(2+(GAMMA-1)*MACH**2)
    M = sqrt( (1+((GAMMA-1)/2)*MACH**2) / (GAMMA*MACH**2-(GAMMA-1)/2) )

    ds = cp*log( (1+((2*GAMMA)/(GAMMA+1))*(MACH**2-1)) * (2+(GAMMA-1)*MACH**2)/((GAMMA+1)*MACH**2) ) - R*log(1+(2*GAMMA/(GAMMA+1))*(MACH**2-1))

    P0 = exp(-ds/R)

    p1 = ( (GAMMA+1)**2*MACH**2 / (4*GAMMA*MACH**2-2*(GAMMA-1)) )**(GAMMA/(GAMMA-1))
    p2 = ( 1 - GAMMA + 2*GAMMA*MACH**2 ) / ( GAMMA + 1 )
    P1 = 1/(p1*p2);

    return (MACH, T, P, RHO, M, P0, P1)
