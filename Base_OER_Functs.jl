#=
File contents:
Mircokinetic model and supporting functions for OER with 4-step adsorbate evolving mechanism.
Standard states used: pH = 0, T = 298 K, p = 1 bar.

Required packages: using CSV, DataFrames, DifferentialEquations, NBInclude, PyPlot, Trapz, DiffEqSensitivity, ForwardDiff 
Written in Julia 1.8.4

Author: Sallye R. Gathmann, University of Minnesota
Funding: This work was supported as part of the Center for Programmable Energy Catalysis, an Energy Frontier Research Center funded by the U.S. Department of Energy, Office of Science, Basic Energy Sciences at the University of Minnesota under award #DE-SC0023464.
DOI of published work: doi.org/10.26434/chemrxiv-2024-gs6zn
=#


#==================================================
  STANDARD FUNCTIONS FOR (ALMOST) ALL SIMULATIONS
==================================================#
function get_constants()
    T = 298.              # [K], temperature
    kB = 8.61733034e-5    # [eV/K], Boltzmann constant
    h = 4.1357e-15        # [eV-s], Planck constant
    e = 1.                # [eV/V], elementary charge

    R = 8.314             # [J/mol-K], gas constant
    F = 96485.3321        # [C/mol], Faraday constant

    Na = 6.0221408e23     # [molecule/mol], Avogadro constant
    Nsites = 1e15         # [site/cm^2], oxide catalyst site density

    return T, kB, h, e, R, F, Na, Nsites
end


"""
    calc_rxn_thermo(E_OH,U)

Calculates thermodynamic parameters based on calatalyst descriptor `E_OH` [eV] and applied potential `U` [V].
Note that in this model, `ΔG2 = E_OH + 0.659 eV`


Returns:
- Free energy of all species, `G` [eV]
- Free energy of each elementary step, `ΔGrxn` [eV]
- Reversible potential of each elementary step, `U_rev` [V]
"""
function calc_thermo(ΔGvolc,U)
    # Import constants
    T, kB, h, e, R, F, Na, Nsites = get_constants()
    
    # Reference potential used in DFT calculations
    U_DFT = 0. # [V]
    
    # Calculate 0V formation energies [eV] using scaling relations
    # Transform ΔGvolc into E_OH to use in scaling relations
    E_OH = ΔGvolc - 0.659
    E_O = 2*E_OH + 0.91
    E_OOH = E_OH + 3.20
    
    # Calculate DFT (U = 0V) reaction free energies [eV]
    ΔG_corr = [0.295, 0.044, 0.377]  # [OH*, O*, OOH*]
    ΔGrxn_DFT = Vector{Float64}(undef,4)
    ΔGrxn_DFT[1] = E_OH + ΔG_corr[1]
    ΔGrxn_DFT[2] = (E_O + ΔG_corr[2]) - (E_OH + ΔG_corr[1])
    ΔGrxn_DFT[3] = (E_OOH + ΔG_corr[3]) - (E_O + ΔG_corr[2])
    ΔGrxn_DFT[4] = 4.92 - sum(ΔGrxn_DFT[1:3])
    # Finally, adjust for potential
    ΔGrxn = ΔGrxn_DFT .- e*(U .- U_DFT)
    
    # Calculate the reversible potential of the reaction
    U_eq = ΔGrxn_DFT ./ e
    
    return ΔGrxn, U_eq
end


"""
    calc_ks_matAgEaeq(ΔGrxn,U_eq,U)

Calculates rate constants for the AEM mechanism assuming that Ea_eq is independent of material and step. Using this as a first approximation of the kinetics.

Potential is accounted for using Butler-Volmer with β = 0.5; max(k) = k_eq.
"""
function calc_ks_matAgEaeq(U_eq,U,ΔGrxn,Ea_eq)
    # Initialize arrays
    kf = Vector{Float64}(undef,4)
    # Import constants
    T, kB, h, e, R, F, Na, Nsites = get_constants()
    
    # Standard values for PCET reactions, material independent:
    A = 1.0e9    # [1/s], kinetic prefector 
    β = 0.5      # [-], transfer coefficient
    
    # Equilibrium constant
    K = exp.(-ΔGrxn/(kB*T))
    
    # Calculate elem. rate constant:
    # Intrinsic/reversible rate constant
    k_eq = A*exp(-Ea_eq/(kB*T))
    # Forward rate constant @ applied electrochemical potential U:
         # added 3rd term to account for ΔGrxn > TS calculated using Ea_eq method
    kf = min.(A, k_eq*exp.((β*e*(U.-U_eq))/(kB*T)), A*exp.(-ΔGrxn/(kB*T)))
    # Reverse rate constant [1/s]
    kr = kf ./ K   
    # Assemble into one vector
    k = vcat(kf,kr)
    
    return k, K

end


"""
    calc_rate(u,p)

Calculates the OER elementary step rates taking inputs `u=[species]` and `p=k`.
"""
function calc_rate(u,k)
    # Activities (defined by standard state)
    a_w = 1.  # water
    a_h = 1.  # proton, H+
    a_O2 = 2.69e-4 / 4.03e-2  # O2(g), ref: Energy Environ. Sci., 2020, 13, 4962.   
    # Rate constants
    kf = k[1:4]
    kr = k[5:8]
    # Rates
    r = Vector(undef,4)
    r[1] = kf[1]*u[1]*a_w  -  kr[1]*u[2]*a_h
    r[2] = kf[2]*u[2]      -  kr[2]*u[3]*a_h
    r[3] = kf[3]*u[3]*a_w  -  kr[3]*u[4]*a_h
    r[4] = kf[4]*u[4]      -  kr[4]*u[1]*a_h*a_O2
    return r
end


"""
    calc_current(r)
Calculates overall OER current density (instantaneous, mA/cm^2) from rate data.
"""
function calc_current(r)
    # Import constants
    T, kB, h, e, R, F, Na, Nsites = get_constants()
    # Calculate OER current density, [mA/cm^2]
    i_OER = F*(Nsites/Na)*sum(r)*1e3
    return i_OER
end


"""
    state_equation(du,u,k,t)

Defines the OER species balance; use with static or dynamic catalysts.  
ODEProblem solution `sol.u = [θ_*, θ_OH*, θ_O*, θ_OOH*]`

- RESTYLED as ODE; just using d/dt[θ_*] as site balance. Works as long as ∑θ_j* = 1 as initial condition.
- User-defined parameter `k` defines the catalyst state. k[1:4] = forward rate constants, k[5:8] = reverse rate constants.
"""
function state_equation(du,u,k,t)
    #=
    Mechanism:
        (1) H2O(l) + *  <-->  OH* + H+ + e-
        (2) OH*  <-->  O* + H+ + e-
        (3) H2O(l) + O*  <-->  OOH* + H+ + e-
        (4) OOH* <--> O2(g) + * + H+ + e-
    =#
    # Calculate elem. rates
    r = calc_rate(u,k)  
    # Species balances
    du[1] = r[4] - r[1] # d/dt θ*
    du[2] = r[1] - r[2]  # d/dt θ_OH*
    du[3] = r[2] - r[3]  # d/dt θ_O*
    du[4] = r[3] - r[4]  # d/dt θ_OOH*
    nothing
end




#==================================================
        FUNCTIONS FOR DYNAMIC SIMULATIONS
==================================================#
"""
    build_squarewave_callbacks(tau,nosc,k)

Builds a descrete callback that will switch the catalyst state (integrator.p) at 
integrator.t ∊ tau for nosc oscillations. 
"""
function build_squarewave_callbacks(tau,nosc,k)
    # Timepoints for callback
    switch_to_S2 = [n for n in tau[1]:sum(tau):nosc*sum(tau)]
    switch_to_S1 = [n for n in sum(tau):sum(tau):nosc*sum(tau)]
    state_switch_times = sort!(vcat(switch_to_S1,switch_to_S2))
    
    # Build individual conditions, affects => callbacks
    # Switching to state 2
    function condition_to_S2(u,t,integrator) 
        t ∈ switch_to_S2
    end
    function affect_to_S2!(integrator)
        integrator.p = big.(k[:,2])
    end
    cb_to_S2 = DiscreteCallback(condition_to_S2,affect_to_S2!); 
    # Switching to state 1
    function condition_to_S1(u,t,integrator) 
        t ∈ switch_to_S1
    end
    function affect_to_S1!(integrator)
        integrator.p = big.(k[:,1])
    end
    cb_to_S1 = DiscreteCallback(condition_to_S1,affect_to_S1!)
    
    # Combine into single discrete callback
    cbs_state_switches = CallbackSet(cb_to_S2,cb_to_S1)
    
    return state_switch_times, cbs_state_switches
end


function build_ss_callbacks(tau,nosc,k,state_switch_times)
    # Some parameters
    check_ss_tstops = state_switch_times[nosc_min*2:nosc_step*2:end] # guarentees there's a tstop in integration @ this timepoint
    rates_check_max = Vector{Float64}(undef,length(check_ss_tstops)) .= 0.
    success = false
    cov_check = false
    
    function condition_check_ss(u,t,integrator)
        if t ∈ check_ss_tstops
            # Preallocate vectors
            idx_mid = Vector{Int}(undef,2)
            idx_end = Vector{Int}(undef,2)
            theta_avg_mid = Vector{BigFloat}(undef,4) 
            theta_avg_end = Vector{BigFloat}(undef,4)
            ss_check = Vector{BigFloat}(undef,4)
            
            # Pull info from integration up to this point
            curr_osc = floor(Int,integrator.t/sum(tau))
            if rem(curr_osc,nosc_step) != 0
                curr_osc = ceil(Int,integrator.t/sum(tau))
            end
            num_osc = 5
            t = integrator.sol.t
            u = transpose(hcat((integrator.sol.u)...))

            # Find switch indices & populate
            switch_indices = find_state_switches(state_switch_times[1:curr_osc*2],t)
            idx_mid[1] = switch_indices[convert(Int,2*curr_osc/2)]+1
            idx_mid[2] = switch_indices[convert(Int,(2*curr_osc/2)+num_osc*2)]  
            idx_end[1] = switch_indices[findfirst(state_switch_times .== t[end])-(num_osc*2)+1]+1
            idx_end[2] = switch_indices[findfirst(state_switch_times .== t[end])+1]
    
            ## Determine if coverages have reached steady-state
            # Calculate time-averaged surface coverages at midpt, endpt  
            for n in 1:4
                theta_avg_mid[n] = trapz(t[idx_mid[1]:idx_mid[2]], u[idx_mid[1]:idx_mid[2],n]) / (t[idx_mid[2]] - t[idx_mid[1]])
                theta_avg_end[n] = trapz(t[idx_end[1]:idx_end[2]], u[idx_end[1]:idx_end[2],n]) / (t[idx_end[2]] - t[idx_end[1]])
            end
    
            # Calculate ss_check value = |x_mid - x_end| / mean(x_mid, x_end) for coverages
            ss_check .= abs.(theta_avg_mid - theta_avg_end) ./ (0.5.*(theta_avg_mid + theta_avg_end))         
            if any(ss_check .> tols[3])
                cov_check == false
                return false
            else
                ## After passing coverage check, check if rates also pass
                cov_check = true
                # Preallocate 
                rates_avg_end = Vector{BigFloat}(undef,4)
                rates_check = Vector{BigFloat}(undef,4)
                r = Matrix{BigFloat}(undef,length(t),4)
                state_vector = Vector{Int64}(undef,length(t))         
                
                # Calculate instantneous rates
                r, sv = calc_rates_squarewave(t,u,state_switch_times[1:curr_osc*2],k)   
     
                # Time-averaged rates at midpt, endpt
                for n in 1:4
                    rates_avg_end[n] = trapz(t[idx_end[1]:idx_end[2]], r[idx_end[1]:idx_end[2],n]) / (t[idx_end[2]] - t[idx_end[1]])
                end
                # Return false if any negatives
                if any(rates_avg_end .< 0) == true
                    return false
                else
                    # Check for steady-state => series reaction so all MUST match or there's accumulation somewhere.
                    for n in 1:4
                        rates_check[n] = (round(convert(Float64,rates_avg_end[1]),sigdigits=4) - round(convert(Float64,rates_avg_end[n]),sigdigits=4)) / (sum(round.(convert.(Float64, rates_avg_end), sigdigits=4))/length(rates_avg_end))
                    end
                    rates_check_max[convert(Int,curr_osc/nosc_step)] = findmax(abs.(rates_check))[1]
                    if any(abs.(rates_check) .> tols[3])
                        return false
                    else
                        success = true
                        return true
                    end
                end
            end       
        
        else
            return false
        end
    end
    
    
    function affect_check_ss!(integrator)
        terminate!(integrator)
    end

    
    function condition_check_rates(u,t,integrator)
        if t ∈ check_ss_tstops[2:end] && cov_check == true
            # Pull info from integration up to this point
            curr_osc = floor(Int,integrator.t/sum(tau))
            if rem(curr_osc,nosc_step) != 0
                curr_osc = ceil(Int,integrator.t/sum(tau))
            end
            
            # Check if rate_check_max is improving
            if abs(rates_check_max[convert(Int,curr_osc/nosc_step)] - rates_check_max[convert(Int,((curr_osc - nosc_step)/nosc_step))]) < 1e-5
                return true
            # Check for oscillatory behavior in rates_check_max
            elseif curr_osc >= nosc_step*3 && abs(rates_check_max[convert(Int,curr_osc/nosc_step)] - rates_check_max[convert(Int,((curr_osc - nosc_step*2)/nosc_step))]) < 1e-5
                return true
            # If we've reached the end & still hasn't met sim tols
            elseif curr_osc == nosc .&& rates_check_max[convert(Int,curr_osc/nosc_step)] > tols[3]
                return true
            else
                return false
            end
        
        else
            return false
        end
    end


    function affect_check_rates!(integrator)
        # Get current tolerance
        curr_tols = [integrator.opts.reltol, integrator.opts.abstol]
    
        # Tighten tolerances for next simulations
        integrator.opts.reltol = curr_tols[1]/100
        integrator.opts.abstol = curr_tols[2]/100

        if success == true
            terminate!(integrator)
        else
            if integrator.opts.abstol < 1e-40
                # Kill if tolerances get too tight
                println("! Maximum tolerance reached; integrator terminated.")
                terminate!(integrator)
            else
                # Reinitialize integrator w/ new tolerances
                reinit!(integrator, erase_sol=true)
                rates_check_max .= 0.
            end
        end
    end
    
    # Assemble into one CallbackSet
    cb_check_ss1 = DiscreteCallback(condition_check_ss, affect_check_ss!)
    cb_check_ss2 = DiscreteCallback(condition_check_rates, affect_check_rates!)
    cb_check_ss = CallbackSet(cb_check_ss1, cb_check_ss2)
    
    return cb_check_ss
end


"""
    find_state_switches(state_switch_times,t)

Returns a vector of t = sol.t indices at which the catalyst state switches.
"""
function find_state_switches(state_switch_times,t)
    # Preallocate
    switch_indices = Vector{Int64}(undef,length(state_switch_times)+1) 
    
    # Starts simulation at catalyst state 1
    switch_indices[1] = 1
    
    # Find other switch points
    for n in range(1,length(state_switch_times))
        switch_indices[n+1] = switch_indices[n+1] = searchsortedfirst(t, state_switch_times[n])
    end
    
    return switch_indices
end


"""
    calc_rates_squarewave(t,theta_inst,state_switch_times,k)

Returns instantaneous elementary rates and catalyst state vector.
"""
function calc_rates_squarewave(t,theta_inst,state_switch_times,k)
    # Create a vector of k(t) based on catalyst, waveform params
    # Initialize arrays
    k_inst = Matrix{BigFloat}(undef,length(t),8)
    state_vector = Vector{Int64}(undef,length(t)); state_vector .= 0
    # Find where in t vector cat state switches
    switch_indices = find_state_switches(state_switch_times,t)
    # Populate a matrix of k(t), catalyst state vector
    for n in range(1,length(state_switch_times),step=2)
        # State 1
        if n > 1
            k_inst[switch_indices[n]+1:switch_indices[n+1],:] .= k[:,1]'
            state_vector[switch_indices[n]+1:switch_indices[n+1]] .= 1
        else
            k_inst[switch_indices[n]:switch_indices[n+1],:] .= k[:,1]'
            state_vector[switch_indices[n]:switch_indices[n+1]] .= 1
        end
        # State 2
        k_inst[switch_indices[n+1]+1:switch_indices[n+2],:] .= k[:,2]'
        state_vector[switch_indices[n+1]+1:switch_indices[n+2]] .= 2
    end
    # Simulation ends at state 1 (fill any remaining #undef in k)
    fin_idx = findfirst(x->!(x in (1,2)), state_vector)[end]
    k_inst[fin_idx:end,:] .= k[:,1]'
    state_vector[fin_idx:end] .= 1
    
    ## CALCULATE INSTANTANEOUS RATES!
    # Initialize arrays
    elem_rates_inst = Matrix{BigFloat}(undef,length(t),4)
    # Calculate elementary rate [1/s]
    for i in 1:length(t)
        elem_rates_inst[i,:] = calc_rate(theta_inst[i,:], k_inst[i,:])
    end
    
    return elem_rates_inst, state_vector
end




#==================================================
           RUNS A STATIC OER SIMULATION
==================================================#
""" runs static OER model """
function run_staticOER(u0,tspan,tols,k)
    # Builds model
    #= Modify coverages using a fudge factor so initial solution point isn't quite as still 
        --> prevents rate[1,:] < 0, which will cause DRC etc calcs to fail. =#
    ff = 1e-6
    u0[1] = u0[1] - ff
    u0[2:4] = u0[2:4] .+ ff/3

    # Set up and solve ODE
    mkm = ODEProblem(state_equation, big.(u0), big.(tspan), big.(k))
    sol = solve(mkm, RadauIIA5(), reltol=tols[1], abstol=tols[2], 
        dtmin=1e-50, maxiters=1e6, dense=false)
    t = sol.t
    u = Array(Array(sol)')
    sol = []

    # Calculate instantaneous rates [1/s]
    len = length(t)
    r = Matrix{BigFloat}(undef,len,4)
    for i in 1:len
        r[i,:] = calc_rate(u[i,:],k)
    end

    # Calculate instantaneous current density [mA/cm^2]
    i_OER = Vector{BigFloat}(undef,len)
    for i in 1:len
        i_OER[i] = calc_current(r[i,:])
    end
    
    return t, u, r, i_OER
end




#==================================================
           RUNS A DYNAMIC OER SIMULATION
==================================================#
""" runs single dynamic OER simulation with user-defined inputs. 
    Catalyst state (ΔGvolc) is oscillated using a square waveform.
""" 
function run_dynamicOER_OscCat(Ea_eq,ΔGvolc,ΔΔG,f,D,η,tols,nosc,Figs=false,threaded=nothing)
    
    # Set up remaining simulation parameters
    tspan = [0., nosc/f]      # [s], integration time limits
    u0 = [1., 0., 0., 0.] # initial coverages, [θ_*, θ_OH*, θ_O*, θ_OOH*]

   
    # Calculate thermodynamic & kinetic parameters
    # Electrochemical potential
    U = 1.23 + η # [V]
    
    # Catalyst thermodynamics, state 1 & state 2
    ΔGrxn1, U_eq1 = calc_thermo((ΔGvolc - ΔΔG/2),U);
    ΔGrxn2, U_eq2 = calc_thermo((ΔGvolc + ΔΔG/2),U);
    
    # Calculate kinetics, state 1 & state 2
    k1, K1 = calc_ks_matAgEaeq(U_eq1,U,ΔGrxn1,Ea_eq);
    k2, K2 = calc_ks_matAgEaeq(U_eq2,U,ΔGrxn2,Ea_eq);
    # Assemble into one array for running dynamic simulation
    k_dyn = hcat(k1,k2)
    
    # Calculate tau
    tau = Vector{Float64}(undef,2)
    tau[1] = (1/f)*D
    tau[2] = (1/f)*(1-D)
    

    # Builds model --> use callbacks to implement dynamics
    global nosc_step = 25
    global nosc_min = 25
    # Get state switch times
    state_switch_times, cbs_state_switches = build_squarewave_callbacks(tau,nosc,k_dyn)
    cb_check_ss = build_ss_callbacks(tau,nosc,k_dyn,state_switch_times)
    cbs = CallbackSet(cbs_state_switches,cb_check_ss)
    
    # Build ODE problem to solve
    mkm = ODEProblem(state_equation, big.(u0), big.(tspan), big.(k_dyn[:,1]))
    
    # Solve model
    sol = solve(mkm, RadauIIA5(), callback=cbs, tstops=state_switch_times, 
    reltol=tols[1], abstol=tols[2], dtmin=1e-50, dtmax=sum(tau)/50, maxiters=1e8, dense=false)
    t_dyn = sol.t
    u_dyn = Array(Array(sol)')
    sol = [];

    
    ## Analyze results
    # Calculate switch indices for returned values
    curr_osc = floor(Int,t_dyn[end]/sum(tau))
    if rem(curr_osc,nosc_step) != 0
        curr_osc = ceil(Int,t_dyn[end]/sum(tau))
    end
    if curr_osc == nosc
        println(" ! Reached end of simulation time; results may not be at steady-state.")
    end
    sstend = findfirst(state_switch_times .== t_dyn[end])
    if typeof(sstend) != Nothing
        state_switch_times = state_switch_times[1:sstend]
        switch_indices = find_state_switches(state_switch_times,t_dyn)

        # Calculate rate
        r_dyn, state_vector = calc_rates_squarewave(t_dyn,u_dyn,state_switch_times,k_dyn)

        # Calcualte current density
        i_OER_dyn = Vector{BigFloat}(undef,size(r_dyn,1)) .= 0
        for i in 1:length(i_OER_dyn)
            i_OER_dyn[i] = calc_current(r_dyn[i,:])
        end

        ## Calculate time-averaged quantities
        # Find range of t_dyn to include
        num_osc = min(5,convert(Int,round(curr_osc/3,digits=0)))
        idx = [switch_indices[convert(Int, length(switch_indices) - num_osc*2)]+1, switch_indices[end]]
        # Preallocate arrays
        u_avg = Vector{Float64}(undef,4) .= 1e10
        r_avg = Vector{Float64}(undef,4) .= 1e10
        for n in 1:4
            u_avg[n] = trapz(t_dyn[idx[1]:idx[2]], u_dyn[idx[1]:idx[2],n]) / (t_dyn[idx[2]] - t_dyn[idx[1]])
            r_avg[n] = trapz(t_dyn[idx[1]:idx[2]], r_dyn[idx[1]:idx[2],n]) / (t_dyn[idx[2]] - t_dyn[idx[1]])
        end
        i_OER_avg = trapz(t_dyn[idx[1]:idx[2]], i_OER_dyn[idx[1]:idx[2]]) / (t_dyn[idx[2]] - t_dyn[idx[1]])

        # Verify time-averaged rates match for all steps & simulation didn't just reach max nosc/told
        rates_check = Matrix{Float64}(undef,4,4)
        for n in 1:4
            rates_check[n,:] = (convert.(Float64, round.(r_avg[n], sigdigits=4)) .- convert.(Float64, round.(r_avg[1:4], sigdigits=4))) / (sum(convert.(Float64, round.(r_avg, sigdigits=4))) / length(r_avg))
        end
        if findmax(rates_check)[1] > tols[3]*1.1
            # Rates didn't reach steady-state => export ridic value (note: can't use zeros or it'll break the run_dynamic function ... inf wait lol rip)
            i_OER_avg = 0
        end

        # If running in non-parallelized loops
        if threaded == nothing       
            # Pull indices that correspond to last 3 oscillations --> should be at steady-state
            idx1 = switch_indices[(length(switch_indices)-4)]+1
            idx2 = switch_indices[(length(switch_indices)-0)]

            # Collect results to return
            # Initialize dataframe for storing time-averaged, steady state data
            labels = ["ΔGvolc", "ΔΔG", "f", "DC", "η", "Ea_eq [eV]", "i_OER", "θ_star", "θ_OH", "θ_O", "θ_OOH", "r_1", "r_2", "r_3", "r_4"]
            df_timeAvg = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]], labels)
            push!(df_timeAvg, hcat(ΔGvolc, ΔΔG, f, D, η, Ea_eq, i_OER_avg, u_avg', r_avg'))
            df_SSinst = DataFrame(hcat(t_dyn[idx1:idx2],state_vector[idx1:idx2],u_dyn[idx1:idx2,:],r_dyn[idx1:idx2,:], i_OER_dyn[idx1:idx2]), vcat(["Time"], ["CatState"], labels[8:end], ["i_OER"]))

            # Return results
            if Figs == true
                # Plot and export figures
                SSplot = plot_dynamic((t_dyn[idx1:idx2] .- t_dyn[idx1]),u_dyn[idx1:idx2,:], r_dyn[idx1:idx2,:], i_OER_dyn[idx1:idx2,:], state_vector[idx1:idx2,:]);
                # Clear out variables
                t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
                return df_timeAvg, df_SSinst, SSplot
            else
                # Clear out variables
                t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
                # dataframes only
                return df_timeAvg, df_SSinst
            end

        # If running in threaded (parallelized) loops
        else
            # Clear out variables
            t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
            return u_avg, r_avg, i_OER_avg
        end

    else # couldn't find state switch time of tend -- weird bug I don't feel like hunting down rn
        if threaded == nothing       

            # Return results
            if Figs == true
                # Clear out variables
                t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
                fig, ax = plt.subplots()
                return DataFrame(), DataFrame(), fig
            else
                # Clear out variables
                t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
                # dataframes only
                return DataFrame(), DataFrame()
            end

        # If running in threaded (parallelized) loops
        else
            # Clear out variables
            t_dyn = []; u_dyn = []; r_dyn = []; state_vector = [];
            return [0,0,0,0], [0,0,0,0], 0
        end
    end
end




#====================================================
        FUNCTIONS FOR (static) DRC SIMULATIONS
    adapted from: https://doi.org/10.1002/aic.17653
====================================================#
function calc_rate_DRC(u,kf)
    # Activities (defined by standard state)
    a_w = 1.  # water
    a_h = 1.  # proton, H+
    a_O2 = 2.69e-4 / 4.03e-2  # O2(g), ref: Energy Environ. Sci., 2020, 13, 4962.   
    # Rate constants
    global K
    kr = kf ./ K
    # Calculates rate
    r = Vector(undef,4)
    r[1] = kf[1]*u[1]*a_w  -  kr[1]*u[2]*a_h
    r[2] = kf[2]*u[2]      -  kr[2]*u[3]*a_h
    r[3] = kf[3]*u[3]*a_w  -  kr[3]*u[4]*a_h
    r[4] = kf[4]*u[4]      -  kr[4]*u[1]*a_h*a_O2
    return r
end


function state_equation_DRC(du,u,kf,t)
    # Calculate elem. rates
    r = calc_rate_DRC(u,kf)  
    # Species balances
    du[1] = r[4] - r[1] # d/dt θ*
    du[2] = r[1] - r[2]  # d/dt θ_OH*
    du[3] = r[2] - r[3]  # d/dt θ_O*
    du[4] = r[3] - r[4]  # d/dt θ_OOH*
    nothing
end


function calc_overall_rate(u,kf)
    # Activities (defined by standard state)
    a_w = 1.  # water
    a_h = 1.  # proton, H+
    a_O2 = 2.69e-4 / 4.03e-2  # O2(g), ref: Energy Environ. Sci., 2020, 13, 4962.
    # Rate constants
    global K
    kr = kf ./ K
    # Overall rate of O2 production
    r = kf[4]*u[4] - kr[4]*u[1]*a_h*a_O2
    
    return r
end


function rate_wrapper(p)
    p = exp.(p)
    _prob = remake(prob, p=p)
    sol =solve(_prob, RadauIIA5(), sensealg=ForwardDiffSensitivity(), reltol=1e-16, abstol=1e-18)  # [nt, nu]
    p_matrix = reshape(p, 1, size(p, 1))  # [1, np]
    # Parse solution => export time datapoint for plotting later
    theta_sol = Array(Array(sol)');
    # Size of solution => not using saveat conditions
    nt = size(theta_sol,1)
    p_repeat = repeat(p_matrix, nt, 1)  # [nt, np]
    r = Array{Real, 1}(undef, nt)
    for i in 1:nt
        r[i] = calc_overall_rate(theta_sol[i, :], p_repeat[i, :])
        if any(r[i] .<= (0))
            r[i] = 1e-20
        end
    end
    return log.(r)
end


function finite_diff(prob,kf,ε)   
    # Start with base case:
    sol_ss = Array(solve(prob, RadauIIA5(), reltol=1e-16, abstol=1e-18))'[end,:];
    r_ss = calc_overall_rate(sol_ss,kf)
    
    # preallocate results vector
    fd_drc = Vector{BigFloat}(undef,4);

    for step in 1:4
        # Nudge kf value
        kf_new = kf[1:4]
        kf_new[step] = kf[step]*(1+ε)
        # solve ODE with nudged kf value
        _prob = remake(prob, p=big.(kf_new))
        sol_ss_new = Array(solve(_prob, RadauIIA5(), reltol=1e-16, abstol=1e-18))'[end,:];
        r_ss_new = calc_overall_rate(sol_ss_new, kf_new)
        # Calculate FD DRC
        if any(r_ss_new .<= (0))
            fd_drc[step] = 1e-20
        else
            fd_drc[step] = (log(r_ss_new) - log(r_ss)) / (log(kf_new[step]) - log(kf[step]))
        end
    end

    return fd_drc
end




#==================================================
              FUNCTIONS FOR PLOTTING
==================================================#
function plot_reaction_coord(ΔGrxn,U,Ea_eq)    
    # Calcualte free energy of each species [eV]
    #  Note that E(H+/e-) is wrapped into E_DFT calcs via CHE 
    G = Matrix{Float64}(undef,5,size(ΔGrxn,2))
    G[1,:] .= 0.               # H2O + *
    for n in 1:size(ΔGrxn,2)
        G[2,n] = G[1,n] + ΔGrxn[1,n]  # OH*
        G[3,n] = G[2,n] + ΔGrxn[2,n]  # O*
        G[4,n] = G[3,n] + ΔGrxn[3,n]  # OOH*
        G[5,n] = G[4,n] + ΔGrxn[4,n]  # O2(g) + *
    end

    # Calculate transition state energies
    TS = Matrix{Float64}(undef,4,size(ΔGrxn,2))
    for n in 1:size(ΔGrxn,2)
        for step in 1:4
            β = 0.5
            e = 1
            if length(U) > 1
                U_eq = ΔGrxn[step,n] .+ e*(U[n] - 0)
                TS[step,n] = max(max(0,Ea_eq - β*e*(U[n] - U_eq)) + G[step,n], G[step+1,n])
            else
                U_eq = ΔGrxn[step,n] .+ e*(U - 0)
                TS[step,n] = max(max(0,Ea_eq - β*e*(U - U_eq)) + G[step,n], G[step+1,n])
            end
        end
    end
               

    # Create figure
    fig, ax = plt.subplots(1,1, figsize = (5.75,4.2), dpi = 100)
    
    # Formatting
    ax.set_ylabel("Gibbs Free Energy  /  eV", fontsize=12)
    ax.set_ylim([findmin(G)[1] - 0.5, findmax(TS)[1] + 0.5])
    ax.set_xticks([2.5,22.5,42.5,62.5,82.5]);
    ax.set_xticklabels([])
    ax.grid(which="major", linestyle=":")
    ax.tick_params(axis="y", which="major", labelsize=11, width=1.1, direction="in", length=4, left="on", right="on", pad=5)
    colors = ["navy"; "darkorange"];
    for axis in ["top","bottom","left","right"]
        ax.spines[axis].set_linewidth(1.)
    end
 
    # Textboxes
    if length(U) == 1
#        ax.text(0.04,0.95,"\$\\eta\$ = $(round(convert(Float64, U - 1.23), sigdigits=3)) V", transform=ax.transAxes, fontsize=11, verticalalignment="top", horizontalalignment="left")
    else
        ctr = ΔGrxn[2,1] + U[1];
#        ax.text(0.04,0.95,"\$ΔG_{2, ctr}^{0V}\$ = $(round(convert(Float64, ctr), sigdigits=3)) eV", transform=ax.transAxes, fontsize=11, verticalalignment="top", horizontalalignment="left")
    end
    bot = -0.08; rot = 0;
    ax.text(-0.01,bot, "H\$_2\$O + *", transform=ax.transAxes, fontsize=12, rotation=rot)
    ax.text(0.26,bot, "HO*", transform=ax.transAxes, fontsize=12, rotation=rot)
    ax.text(0.47,bot, "O*", transform=ax.transAxes, fontsize=12, rotation=rot)
    ax.text(0.66,bot, "HOO*", transform=ax.transAxes, fontsize=12, rotation=rot)
    ax.text(0.88,bot, "O\$_{2(g)}\$ + *", transform=ax.transAxes, fontsize=12, rotation=rot)
    
    
    lw = [2., 1.75, 0.75]
        # Initialize so legend is correct
        ax.plot([-1,6],  [0,0], color=colors[1], linewidth=lw[2], label="State 1")
        if size(ΔGrxn,2) > 1
            ax.plot([-1,6],  [0,0], color=colors[2], linewidth=lw[2], label="State 2")
            ax.legend(loc="lower left", bbox_to_anchor=(0.01,0.01), ncol=1, fontsize=11, frameon=false)
        end
    # Plot remainder  
    for n in range(1,size(ΔGrxn,2),step=1)  
        # Plot stuff
        x = [0,19,39,59,79]
        for m in 2:5
            # Intermediate energies
            ax.plot([x[m],x[m]+7],  [G[m,n],G[m,n]], color=colors[n], linewidth=lw[1])
        end
        
        # Changed transition states so they aren't plotted if there is no fwd activation barrier
        x1 = [10,30,50,70];
        for m in 1:4
            ax.plot([x1[m],x1[m]+5],[TS[m,n],TS[m,n]], color=colors[n], linewidth=lw[2])
            ax.plot([x1[m]-4,x1[m]], [G[m,n],TS[m,n]],  color=colors[n], linewidth=lw[3], linestyle="dashed")
            ax.plot([x1[m]+5,x1[m]+9],[TS[m,n],G[m+1,n]], color=colors[n], linewidth=lw[3], linestyle="dashed")
            ax.text(x1[m]+1, TS[m,n]+0.05, "\$\\ddag_{$(m)}\$", color=colors[n])
        end

    end
    
    plt.tight_layout()
    return fig     
end


function plot_static(t,theta_inst,rates_inst,cd_inst)
    
    clrs = ["indigo","mediumvioletred", "salmon","maroon"]
    st = ["solid", "dashed", "solid", "dashdot"]
    plt.rcParams["test.usetex"] = "True"

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10,4.2))

    for n in 1:4
        ax1.plot(t,theta_inst[:,n], color=clrs[n], linestyle=st[n], linewidth=1.5)
        ax2.plot(t,rates_inst[:,n], color=clrs[n], linestyle=st[n], linewidth=1.5)
    end

    ax1.set_ylabel("Surface Coverage, \$\\theta_{j^*}\$  /  -", fontsize=12)
    ax1.legend(ncol=4, ["\$\\theta_{*}\$", "\$\\theta_{OH^*}\$",
        "\$\\theta_{O^*}\$", "\$\\theta_{OOH^*}\$"], 
        fontsize=11, loc = "lower center", bbox_to_anchor=(0.5,1.),
        edgecolor="none")
    ax2.set_ylabel("Elementary Rates, \$r_i\$  /  s\$^{-1}\$", fontsize=12)
    ax2.legend(ncol=4, ["\$r_1\$", "\$r_2\$", "\$r_3\$", "\$r_4\$"],
        fontsize=11, loc = "lower center", bbox_to_anchor=(0.5,1.),
        edgecolor="none")
    ax2.set_yscale("log")

    for axes in [ax1,ax2]
        axes.set_xscale("log")
        axes.set_xlabel("Time on Stream  /  s", fontsize=12)
        axes.grid(which="major", linestyle=":")
        for axis in ["top","bottom","right","left"]
            axes.spines[axis].set_linewidth(1.)
            axes.tick_params(axis="both", which="major", labelsize=11, 
                width=1.1, direction="in", length=4, left="on", right="on", pad=5)
        end 
    end

    # Twinx for current density
    ax2a = ax2.twinx()
    ax2a.plot(t,cd_inst,color="b",linestyle="-",linewidth=2,alpha=0.15) # easier to see : style
    ax2a.plot(t,cd_inst,color="b", linestyle=":",linewidth=1.75)
    ax2a.set_yscale("log")
    ax2a.tick_params(axis="y", which="major", color="b", labelsize=11,
        width=1, direction="in", length=4, pad=5)
    ax2a.tick_params(axis="y", which="minor", color="b", labelsize=11,
        width=1, direction="in", length=2, pad=5)
    ax2a.spines["right"].set_color("b")
    ax2a.tick_params(colors="b")
    ax2a.set_ylabel("Current Density, \$i\$\$_{OER}\$  /  mA cm\$^{-2}\$", color="b", fontsize=12)

    plt.tight_layout()
    
    return fig
end


function plot_DRC_results(ode_sol,ad_drc)
    # Line colors, styles
    clrs = ["indigo","mediumvioletred", "salmon","maroon"]
    st = ["solid", "dashed", "solid", "dashdot"]
    
    # Create figure
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(10,4.2))
    
    # Plot normal ODE solution
    for n in 1:4
        ax1.plot(ode_sol.t, Array(Array(ode_sol)')[:,n],
            color=clrs[n], linestyle=st[n], linewidth=1.5)
    end
    # Format
    ax1.set_xlabel("Time  /  s", fontsize=12)
    ax1.set_xscale("log")
    ax1.set_xticks([1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1e1])
    ax1.set_xlim([1e-12,1e2])
    ax1.grid(which="major", linestyle=":")
    ax1.set_ylabel("Coverages, \$\\theta_{j^*}\$  /  -", fontsize=12)
    ax1.legend(ncol=4, ["\$\\theta_{*}\$", "\$\\theta_{OH^*}\$",
        "\$\\theta_{O^*}\$", "\$\\theta_{OOH^*}\$"], 
        fontsize=11, loc = "lower center", bbox_to_anchor=(0.5,1.),
        edgecolor="none")
   
    # Plot AD DRC results --> bar plot
    ax2.bar([1,2,3,4], height=ad_drc[end,:], color="darkcyan", edgecolor="darkslategray")
    # Format
    ax2.set_ylabel("Degree of Rate Control, \$X_{RC,i}\$  /  -", fontsize=12)
    ax2.set_xscale("linear")
    ax2.grid(which="major", linestyle=":")
    ax2.set_xlabel("Step ID  /  -", fontsize=12)
    ax2.text(0.5, 1.05, "∑\$X_{RC,i}\$ = $(round(convert(Float64,sum(ad_drc[end,:])), sigdigits=2))",
        horizontalalignment="left")
    ax2.set_ylim([0,1])
    ax2.set_xticks([1,2,3,4])
    
    # Last formatting
    for axes in [ax1,ax2]
        for axis in ["top","bottom","right","left"]
            axes.spines[axis].set_linewidth(1.)
            axes.tick_params(axis="both", which="major", labelsize=11, 
                width=1.1, direction="in", length=4, left="on", right="on", pad=5)
        end 
    end
    plt.tight_layout()
    
    return fig
end


function plot_dynamic(t,u,r,i_OER,state_vector)
    
    clrs = ["indigo","mediumvioletred", "salmon","maroon"]
    st = ["solid", "dashed", "solid", "dashdot"]
    plt.rcParams["test.usetex"] = "True"

    fig, (ax3, ax1, ax2) = plt.subplots(1,3, figsize=(15,4.2))

    for n in 1:4
        ax1.plot(t,u[:,n], color=clrs[n], linestyle=st[n], linewidth=1.5)
        ax2.plot(t,r[:,n], color=clrs[n], linestyle=st[n], linewidth=1.5)
    end
    ax3.plot(t,state_vector, color="grey", linestyle="-", linewidth=1.5)

    # Formatting
    ax1.set_ylabel("Surface Coverage, \$\\theta_{j^*}\$  /  -", fontsize=12)
    ax1.legend(ncol=4, ["\$\\theta_{*}\$", "\$\\theta_{OH^*}\$",
        "\$\\theta_{O^*}\$", "\$\\theta_{OOH^*}\$"], 
        fontsize=11, loc = "lower center", bbox_to_anchor=(0.5,1.),
        edgecolor="none")
    ax2.set_ylabel("Elementary Rates, \$r_i\$  /  s\$^{-1}\$", fontsize=12)
    ax2.legend(ncol=4, ["\$r_1\$", "\$r_2\$", "\$r_3\$", "\$r_4\$"],
        fontsize=11, loc = "lower center", bbox_to_anchor=(0.5,1.),
        edgecolor="none")
    ax2.set_yscale("log")
    ax3.set_ylabel("Catalyst State, S  /  -", fontsize=12)
    ax3.set_yticks([])

    for axes in [ax1,ax2,ax3]
        axes.set_xlabel("Time  /  -", fontsize=12)
        axes.grid(which="major", linestyle=":")
        for axis in ["top","bottom","right","left"]
            axes.spines[axis].set_linewidth(1.)
            axes.tick_params(axis="both", which="major", labelsize=11, 
                width=1.1, direction="in", length=4, left="on", right="on", pad=5)
        end 
    end
    
    # Twinx for current density
    ax2a = ax2.twinx()
    ax2a.plot(t,i_OER,color="b",linestyle="-",linewidth=2,alpha=0.15) # easier to see : style
    ax2a.plot(t,i_OER,color="b", linestyle=":",linewidth=1.75)
    ax2a.tick_params(axis="y", which="major", color="b", labelsize=11,
        width=1, direction="in", length=4, pad=5)
    ax2a.tick_params(axis="y", which="minor", color="b", labelsize=11,
        width=1, direction="in", length=2, pad=5)
    ax2a.spines["right"].set_color("b")
    ax2a.tick_params(colors="b")
    ax2a.set_ylabel("Current Density, \$i\$\$_{OER}\$  /  mA cm\$^{-2}\$", color="b", fontsize=12)
    ax2a.set_yscale("log")

    plt.tight_layout()


    return fig
end