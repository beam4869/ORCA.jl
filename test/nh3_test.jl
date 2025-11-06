# test/nh3_test.jl
using Test
using ORCA
using JuMP
using Statistics
using XLSX, DataFrames

@testset "NH3 toy test" begin
    # --- Load small test data (relative to this file) ---
    xlsx_path = joinpath(@__DIR__, "Jan_2022.xlsx")
    @test isfile(xlsx_path)  # ensure the file exists in test/
    df1 = DataFrame(XLSX.readtable(xlsx_path, "test_for_pareto_frontier"))
    ceb_total = Float64.(df1[!, "Jan 21st 11pm price"])
    gridemm_total = Float64.(df1[!, "Jan 21st 11pm emission"])
    

    # --- Build your model (trimmed for brevity) ---
    tl = 48; il = 3
    t = collect(1:tl)
    wind = zeros(tl); wpow = zeros(tl)
    avgwind = mean(wind)  # just to confirm Statistics is available

    rho   = [60.0, 0.8, 2.12]
    mmin  = [9.99, 0.0, 750.0]
    mmax  = [33.3, 2250.0, 1100.0]
    nemax = 15
    msmin = [4035.0, 22601.0, 0.0]
    msmax = [807197.0, 4.52e6, 0.0]
    mnh3tar = 1000.0
    rmin, rmax, rampt, damax = 100.0, 40.0, 4, 1.0
    ceb = copy(ceb_total[1:tl])
    gridemm = copy(gridemm_total[1:tl])
    println("ceb: ", ceb)
    println("gridemm: ", gridemm)
    cnh3, con, cel, cpsa = 1.46, 50.0, 0.42, 0.10
    cbh2, h2perh = 2.30, 100.0
    mi0 = copy(msmin); ne0 = 0; mnh30 = mnh3tar
    da0 = [0, 0, 0]
    h2fossil=9.3  # carbon emission factor from hydrogen derived from steam reforming of methane. units: kg CO2/Kg H2
    el_start_risk=5
    el_op_risk=1

    nh3w = Model()
    set_silent(nh3w)  # Suppress solver output (optional)

    @variable(nh3w, da[-2:tl], Bin)
    @variable(nh3w, de[1:tl] >= 0, Int)
    @variable(nh3w, ms[1:il, 0:tl] >= 0)
    @variable(nh3w, mp[1:il, 0:tl] >= 0)
    @variable(nh3w, mnh3dev[1:tl])
    @variable(nh3w, ne[0:tl], Int)
    @variable(nh3w, peb[1:tl] >= 0)
    @variable(nh3w, pei[1:il, 1:tl] >= 0)
    @variable(nh3w, pes[1:tl] >= 0)
    @variable(nh3w, u[1:il] >= 0)
    @variable(nh3w, bh2[1:tl] >= 0)

    #variables lUB
    @constraint(nh3w, demin[i in 1:tl], -de[i]<=0)
    @constraint(nh3w, msmini[j in 1:il,i in 1:tl], -ms[j,i]<=0)
    @constraint(nh3w, mpmin[j in 1:il,i in 1:tl], -mp[j,i]<=0)
    @constraint(nh3w, pebmin[i in 1:tl], -peb[i]<=0)
    @constraint(nh3w, peimin[j in 1:il,i in 1:tl], -pei[j,i]<=0)
    @constraint(nh3w, pesmin[i in 1:tl], -pes[i]<=0)
    @constraint(nh3w, umin[i in 1:il], -u[i]<=0)
    @constraint(nh3w, bh2min[i in 1:tl], -bh2[i]<=0) # confusing...


    #steady constraints. section 2.2
    @constraint(nh3w, energybal[i in 1:tl], wpow[i]+peb[i]-pes[i]-pei[1,i]-pei[2,i]-pei[3,i]==0) # energy balance
    @constraint(nh3w, energyunt[i in 1:il, j in 1:tl], pei[i,j] == rho[i]*mp[i,j]) #the power used by different chemical units is propotional to the mass of the materials
    @constraint(nh3w, deviation[j in 1:tl], mnh3dev[j] == mnh3tar-mp[3,j]) #economic penalty for underproduction

    #dynamic constraints, section 2.3
    @constraint(nh3w, h2stor[j in 1:tl], ms[1,j]==ms[1,j-1]+bh2[j]+mp[1,j]-(3/17)*mp[3,j]) #ms[i,0] the iteration process of H2
    @constraint(nh3w, n2stor[j in 1:tl], ms[2,j]==ms[2,j-1]+mp[2,j]-(14/17)*mp[3,j]) # the iteration process of N2
    @constraint(nh3w, totalstor[i in 1:il], u[i]==ms[i,0]-ms[i,48]) #ms[i,0] # the total depleted H2 and N2 gas during the production period.

    #steady constraints. section 2.2
    @constraint(nh3w, minflow[i in 2:il, j in 1:tl], mmin[i]-mp[i,j]<=0)
    @constraint(nh3w, maxflow[i in 2:il, j in 1:tl], -mmax[i]+mp[i,j]<=0)
    @constraint(nh3w, minelec[j in 1:tl], -ne[j]<=0)
    @constraint(nh3w, maxelec[j in 1:tl], ne[j]-nemax<=0)
    @constraint(nh3w, minh2pr[j in 1:tl], mmin[1]*ne[j]- mp[1,j]<=0)
    @constraint(nh3w, maxh2pr[j in 1:tl], mp[1,j]-mmax[1]*ne[j]<=0)
    @constraint(nh3w, minstor[i in 1:il, j in 1:tl], msmin[i]-ms[i,j]<=0)
    @constraint(nh3w, maxstor[i in 1:il, j in 1:tl], -msmax[i]+ms[i,j]<=0)
    @constraint(nh3w, maxbh2[j in 1:tl], bh2[j]-h2perh<=0) # should be equvlient as before, nu still in a different form...

    #dynamic constraints, section 2.3
    @constraint(nh3w, elecon[j in 1:tl], -de[j]+ne[j]-ne[j-1]<=0) #ne[0]
    @constraint(nh3w, rampmin[j in 1:tl], -da[j]*rmin - mp[3,j]+mp[3,j-1]<=0) #mp[i,0]
    @constraint(nh3w, rampmax[j in 1:tl],  mp[3,j]-mp[3,j-1] - da[j]*rmax<=0)
    @constraint(nh3w, ramplim[j in 1:tl], sum(da[i] for i in (j+1-rampt):j)-damax<=0) #da[-2,-1,0]

    #dynamic constraint feed-in data constraints
    @constraint(nh3w, ms_t0[i in 1:il], ms[i,0] == mi0[i]) #make the initial into the variable.
    @constraint(nh3w, ne_t0, ne[0]==ne0) #make the initial into the variable.
    @constraint(nh3w, mnh3_t0, mp[3,0]==mnh30) #make the initial into the variable.
    @constraint(nh3w, da_t[j in -2:0], da[j]==da0[j+3]) 


    @expression(nh3w, costs, sum(sum(ceb[j]*rho[i]*mp[i,j] for i in 1:il)
        +cel*mp[1,j]+cbh2*bh2[j]+cpsa*mp[2,j] for j in 1:tl))
    @expression(nh3w, emms, 
        sum(sum(rho[j]*mp[j,i]*gridemm[i] for j in 1:il) for i in 1:tl)
        +sum(h2fossil*bh2[i] for i in 1:tl))
    @expression(nh3w, watuse, sum(mp[1,j]*9.0 for j in 1:tl)+sum(bh2[j]*3 for j in 1:tl)) # water usage; 
    @expression(nh3w, safetyelectro, sum(el_op_risk*(mp[1,i])/mmax[1] for i in 1:tl))

    # --- Call ORCA entry point ---
    res = ORCA.main(nh3w, [costs, emms, watuse, safetyelectro])
    println("ORCA.main returned: ", res.groups)

    @test !isnothing(res)
end


