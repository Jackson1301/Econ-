module KMCournotSp2021Assignment

import Random
import Distributions; const Dist = Distributions
import StatsBase; const Sb = StatsBase
import DataFrames; const Df = DataFrames
import PrettyTables; const Pt = PrettyTables

"""
You don't need to change this function.

Market demand function: A - b * p.

"""
function paramvalues(; numfirms = 2,
                       epochlens = Dict(i => 30 for i = 1:numfirms),
                       initqs = Dict(i => 40 for i = 1:numfirms),
                       deltas = Dict(i => 3.0 for i = 1:numfirms),
                       epsilons = Dict(i => 0.7 for i = 1:numfirms),
                       costs = Dict(i => 0.0 for i = 1:numfirms),
                       A = 400,
                       b = 2,
                       T = 6600)


    Dict(
        :numfirms => numfirms,
        :epochlens => epochlens,
        :initqs => initqs,
        :deltas => deltas,
        :epsilons => epsilons,
        :costs => costs,
        :A => A, :b => b, :T => T
    )
end

"""
You don't need to change this function.

Calculate the Cournot-Nash outcome

Return a dictionary consisting of:
    :firmqtys   - Each firms output in the Cournot Nash equilibrium.
    :Q          - Aggregate output in the Cournot Nash equilibrium.
    :P          - Market price in the Cournot Nash equilibrium.

"""
function cournot_nash(pv)
    A = pv[:A]
    b = pv[:b]
    costs = pv[:costs]
    aggcosts = costs |> values |> sum
    nfirms = length(costs)
    qtynum(i) = A + aggcosts - (nfirms + 1) * costs[i]
    qtyden = (nfirms + 1) * b
    firmqtys = Dict(i => qtynum(i)/qtyden for i = 1:nfirms)

    res = Dict{Symbol, Any}(:firmqtys => firmqtys)
    res[:Q] = firmqtys |> values |> sum

    res[:P] = A - b * res[:Q]
    return res
end

"""
Calculate firm's profit given:
    firmidx - one of {1, 2, ..., numfirms}.
    firmqty - qty produced by the firm.
    aggqty  - aggregate quantity in the market.
    pv      - paramvalues

Return value is of type Float64.
"""
function calc_firm_profit(firmidx, firmqty, aggqty, pv)
    A = pv[:A]
    b = pv[:b]
    Q = aggqty
    ## calculate price using the inverse demand
    price = A - b*Q
    ## get firms cost of production form pv using firmidx
    cost = pv[:costs][firmidx]
    ## calculate and return profit.
    profit = (price - cost)*firmqty
    return Float64(profit)
end

# pv = paramvalues()
# c = calc_firm_profit(2, 30, 70, pv)


"""
Draw the bidQuantity given current quantity (currq) and delta.
bidQuantity is drawn from uniform distribution U[currq - delta, currq + delta].

Return value is of type Float64
"""
function bidquantity(currq, delta)
    ## Define the uniform distribution with the correct parameters.
    unif_dist = Dist.Uniform(currq-delta, currq+delta)
    ## return a random value from this uniform distribution
    # Random.seed!(1)
    return Float64(rand(unif_dist))
end

# bidquantity(30, 3.0)

"""
Calculate profits of all firms given the quantities dictionary and paramvalues.
quantity dictionary has firm index (1,2,3...) as keys and the
quantity (Float64) as values. (firmidx (Int64) => quantity (Float64))

Return value: profits (Dictionary with firm index as keys and the profits (Float64) as values.)
              The keys should go from 1,2,...,numfirms
"""
function calc_profits(quantities, pv)
    ## calculate aggregate quantity using the quantities dictionary.
    aggqty = Sb.sum(values(quantities))

    ## use the aggregate quantity and quantities dict to calculate and
    ## return profits of all firms (see the return value described above)
	# profits = Dict()
    # for i in keys(quantities)
    #     profit = calc_firm_profit(i, quantities[i], aggqty, pv)
    #     push!(profits, i => profit)
    # end
    profits = Dict(i => Float64(calc_firm_profit(i, quantities[i], aggqty, pv)) for i in keys(quantities))
    return profits
end

# pv = Dict(
#     :T => 6600,
# 	:b => 2,
# 	:numfirms => 2,
# 	:initqs => Dict(2 => 40, 1 => 40),
# 	:A => 400,
# 	:deltas => Dict(2 => 3.0, 1 => 3.0),
# 	:epochlens => Dict(2 => 30, 1 => 30),
# 	:epsilons => Dict(2 => 0.7, 1 => 0.7),
# 	:costs => Dict(2 => 5.0, 1 => 10.0)
# )
# calc_profits(Dict(2 => 45.0, 1 => 85.0), pv)

"""
Calculate the amount by which the current quantity should be changed
(increased or decreased) by firm `firmidx`, given the `firmreturnsdict` and
`pv` (paramvalues).
`firmreturnsdict` is a dictionary in the following format
    :up     => Array of returns when firm produced a quantity the same or higher than its current currentQuantity
    :down   => Array of returns when firm produced a quantity lower than its current currentQuantity

Return value is of type Float64.
"""
function calc_currqchange(firmidx, firmreturnsdict, pv)
    ## calculate average profit from the firmreturnsdict[:up] array.
	if isempty(firmreturnsdict[:up])
		avg_prof_up = 0
	else
		avg_prof_up = Sb.mean(firmreturnsdict[:up])
	end
	## calculate average profit from the firmreturnsdict[:down] array
	if isempty(firmreturnsdict[:down])
		avg_prof_down = 0
	else
		avg_prof_down = Sb.mean(firmreturnsdict[:down])
	end
    # avg_prof_up = Sb.mean(firmreturnsdict[:up]).
    # avg_prof_down = Sb.mean(firmreturnsdict[:down])
    ## if same or the former is bigger return the epsilon for this firm from pv
    ## if the latter is bigger return -1 * epsilon for this firm from pv
    if avg_prof_up >= avg_prof_down
        return Float64(pv[:epsilons][firmidx])
    else
        return Float64(-1*pv[:epsilons][firmidx])
    end
end

# pv = Dict(
#     :T => 6600,
# 	:b => 2,
# 	:numfirms => 2,
# 	:initqs => Dict(2 => 40, 1 => 40),
# 	:A => 400,
# 	:deltas => Dict(2 => 3.0, 1 => 3.0),
# 	:epochlens => Dict(2 => 30, 1 => 30),
# 	:epsilons => Dict(2 => 0.6, 1 => 0.8),
# 	:costs => Dict(2 => 0.0, 1 => 0.0)
# )
# firmreturnsdict = Dict(:up => [32.43029950291791, 68.85052016067101, 86.97162407893, 24.98321328972841, 34.82407658654518, 6.602576509850033, 31.828621996325325, 38.6709458875224, 0.3996605889154514, 57.944561918710846], :down => [2.535956655709395, 28.373583343973543, 40.22865954750362, 14.79136923237368, 68.04999906060488, 19.223985296945614, 56.45116985159268])
# calc_currqchange(2, firmreturnsdict, pv)

"""
Generate new currqtys dictionary and returnddict dictionary for all firms.
How? In the following manner:
Loop over each firmidx. For each firm check if it is a new epoch (epoch lengths
can vary across firms).

    If it is not a new epoch then keep the currqty the same for the firmidx and keep
    the returnsdict of the firm the same.

    If it is a new epoch, calculate the new currqty using the current currqty and
    calc_currqchange function. Also, create new firmreturndict with empty :up and :down
    arrays.

To check whether it is a new epoch for a firm, calculate (`t` % firm's epoch length).
It is a new epoch if this values is 1 and t != 1.

"""

function change_currqtys_reset_returnsdict(t, currqtys, returnsdict, pv)
    newcurrqtys = Dict()
    newreturnsdict = Dict()
    epochlens = pv[:epochlens]
    for firmidx in keys(currqtys)
        ## Is `t` the start of a new epoch? If yes,
        if (t % epochlens[firmidx]) == 1 && t != 1
            ## calculate new currqty
            currqty = Float64(currqtys[firmidx]) + calc_currqchange(firmidx, returnsdict[firmidx], pv)
            push!(newcurrqtys, firmidx => Float64(currqty))
            ## reset firm's retursdict
            push!(newreturnsdict, firmidx => Dict(:up => [], :down => []))
			# replace!(newreturnsdict, firmidx => Dict(:up => [], :down => []))
        else
        ## if not,
            ## add current currqty to newcurrqtys as it is
            push!(newcurrqtys, firmidx => Float64(currqtys[firmidx]))
            ## add current firm's returnsdict to newreturnsdict as it is.
            push!(newreturnsdict, firmidx => returnsdict[firmidx])
		end
    end
    return newcurrqtys, newreturnsdict
end

# change_currqtys_reset_returnsdict(11, Dict(2 => 45.0, 1 => 53.0),
# 	Dict{Any, Any}(
# 	2 => Dict(:up => [20.947237319807076, 25.137920979222493, 2.037486871266725, 28.77015122756894, 85.9512136087661, 7.695088688120899, 64.03962459899388, 87.35441302706855], :down => [27.85824200287785, 75.13126327861701, 64.4883353942093, 7.782644396003469, 84.81854810000327, 8.56351682044918, 55.32055454580578]),
# 	1 => Dict(:up => [23.603334566204694, 34.651701419196044, 31.27069683360675, 0.790928339056074, 48.86128300795012, 21.096820215853597, 95.1916339835734, 99.99046588986135], :down => [25.166218303197184, 98.66663668987997, 55.57510873245724, 43.71079746096251, 42.471785049513144, 77.3223048457377, 28.11902322857298])),
# 	Dict(
#     	:T => 6600,
# 		:b => 2,
# 		:numfirms => 2,
# 		:initqs => Dict(2 => 40, 1 => 40),
# 		:A => 400,
# 		:deltas => Dict(2 => 3.0, 1 => 3.0),
# 		:epochlens => Dict(2 => 30, 1 => 10),
# 		:epsilons => Dict(2 => 0.7, 1 => 0.7),
# 		:costs => Dict(2 => 0.0, 1 => 0.0)
# ))

"""
Runs one episode.
Returns
    - new currqtys dict
    - newreturnsdict
    - bidquantity

    For the first two return values, uses the function defined above.

    currqtys is a Dict of type firmidx => firm's current quantity
    returnsdict is Dict of type firmidx => (Dict of :up|:down => returns array)
    quantities is a Dict: firmidx => Array of quantities produced in the epoch
"""
function run_episode(t, currqtys, returnsdict, pv)
    newcurrqtys, newreturnsdict = change_currqtys_reset_returnsdict(t, currqtys, returnsdict, pv)

    bidqs = Dict()
    for firmidx in keys(newcurrqtys)
        ## populate bidqs with each firms bidQuantity for this episode
		push!(bidqs, firmidx => bidquantity(newcurrqtys[firmidx],pv[:deltas][firmidx]))
    end

    ## calculate profits of all firms in this episode.
	profits = calc_profits(currqtys, pv)

    for (firmidx, newcurrqty) in newcurrqtys
        ## update the newreturnsdict dictionary appropriately with the
        ## calculated profits.

		# if bidqs[firmidx] >= currqtys[firmidx]

		if bidqs[firmidx] >= currqtys[firmidx]
			push!(newreturnsdict[firmidx][:up], bidqs[firmidx])
		else
			push!(newreturnsdict[firmidx][:down], bidqs[firmidx])
		end
    end
    return (newcurrqtys, newreturnsdict, bidqs)
end

# Random.seed!(50)
# run_episode(11, Dict(2 => 45.0, 1 => 53.0),
# 	Dict{Any, Any}(
# 	2 => Dict(:up => [20.947237319807076, 25.137920979222493, 2.037486871266725, 28.77015122756894, 85.9512136087661, 7.695088688120899, 64.03962459899388, 87.35441302706855], :down => [27.85824200287785, 75.13126327861701, 64.4883353942093, 7.782644396003469, 84.81854810000327, 8.56351682044918, 55.32055454580578]),
# 	1 => Dict(:up => [23.603334566204694, 34.651701419196044, 31.27069683360675, 0.790928339056074, 48.86128300795012, 21.096820215853597, 95.1916339835734, 99.99046588986135], :down => [25.166218303197184, 98.66663668987997, 55.57510873245724, 43.71079746096251, 42.471785049513144, 77.3223048457377, 28.11902322857298])),
# 	Dict(
#     	:T => 6600,
# 		:b => 2,
# 		:numfirms => 2,
# 		:initqs => Dict(2 => 40, 1 => 40),
# 		:A => 400,
# 		:deltas => Dict(2 => 3.0, 1 => 3.0),
# 		:epochlens => Dict(2 => 30, 1 => 10),
# 		:epsilons => Dict(2 => 0.7, 1 => 0.7),
# 		:costs => Dict(2 => 0.0, 1 => 0.0)
# ))

"""
Returns the quantities array for the market simulation (all epochs and all
episodes in one market).
"""
function simulate_mkt(pv)
    quantities = Dict(i => [] for i = 1:pv[:numfirms])
    ## initialize the currqtys dictionary
	currqtys = Dict(i => pv[:initqs][i] for i = 1:pv[:numfirms])
    ## initialize returnsdict dictionary for each firm
	returnsdict = Dict(i => Dict(:up => [], :down => []) for i = 1:pv[:numfirms])
    for t = 1:pv[:T]
        ## obtain new currqtys and new returnsdict and bidqs for episode t
		currqtys, returnsdict, bidqs = run_episode(t, currqtys, returnsdict, pv)
        ## update quantities dictionary using the bidqs data.
		for i in keys(quantities)
			push!(quantities[i], bidqs[i])
		end
    end
    return quantities
end

pv = Dict(
    :T => 1500,
	:b => 2,
	:numfirms => 2,
	:initqs => Dict(2 => 40, 1 => 40),
	:A => 400,
	:deltas => Dict(2 => 3.0, 1 => 3.0),
	:epochlens => Dict(2 => 30, 1 => 10),
	:epsilons => Dict(2 => 0.7, 1 => 0.7),
	:costs => Dict(2 => 5.0, 1 => 0.0))
Random.seed!(50)
s = simulate_mkt(pv)
Sb.mean(s[1])
Sb.mean(s[2])

"""
Return value looks like (given N firms)
    Dict(
        :q1     => Average of last obs_to_use observation from firm 1's quantities array,
        :q2     => Average of last obs_to_use observation from firm 2's quantities array,
        :q3     => Average of last obs_to_use observation from firm 3's quantities array,
        .
        .
        .
        :qN     => Average of last obs_to_use observation from firm N's quantities array,
        :aggQ   => Sum of values of the keys :q1, :q2, :q3 etc. described above. This is the
                   average aggregate quantity.
    )
"""
function calc_avg_and_aggavg_qty(quantities, obs_to_use)

    ## Create a dictionary that contains keys such as
    ## q1, q2, q3 etc. where the number denotes the firmidx.
    ## The value should be the average of firm's output in the *last* obs_to_use observations.
	qtys = Dict()
	for firmidx in 1:length(quantities)
		obs = []
		for i in (length(quantities[firmidx]) - obs_to_use + 1):length(quantities[firmidx])
			# push!(obs, quantities[firmidx][i])
			# push!(obs, Sb.sample(quantities[firmidx]))
			push!(obs, rand(quantities[firmidx]))
		end
		# print(obs)
		c = string("q",firmidx)
		push!(qtys, Symbol(c) => Sb.mean(obs))
	end
    ## find the average aggregate quantity, :aggQ, using the aggregate of firm qty averages.
	agg_avg = Sb.sum(values(qtys))
    ## add this to dictionary above with key :aggQ
	push!(qtys, :aggQ => agg_avg)
    ## return the dictionary
	return qtys
end

pv = Dict(
    :T => 2000,
	:b => 2,
	:numfirms => 2,
	:initqs => Dict(2 => 40, 1 => 40),
	:A => 400,
	:deltas => Dict(2 => 3.0, 1 => 3.0),
	:epochlens => Dict(2 => 30, 1 => 30),
	:epsilons => Dict(2 => 0.7, 1 => 0.7),
	:costs => Dict(2 => 0.0, 1 => 0.0))
Random.seed!(10)
s = simulate_mkt(pv)
c = calc_avg_and_aggavg_qty(s,500)

"""
Should not have to change any code here...
"""
function run_reps(numreps, pv)
    ## This is the empty dataframe to store data
    df = Df.DataFrame( rseed = Int64[], aggQ = Float64[] )
    numfirms = pv[:numfirms]
    for firmidx in 1:numfirms
        df[:, Symbol("q$(firmidx)")] = Float64[]
    end

    ## number of observation used for calculating average
    obs_to_use = 1000
    for rseed in 1:numreps
        # println("Running rep $(rseed)")
        Random.seed!(rseed)
        quantities = simulate_mkt(pv)
        d = Dict(:rseed => rseed)
        d = merge(d, calc_avg_and_aggavg_qty(quantities, obs_to_use))
        push!(df, d)
    end
    return df
end

"""
Should not have to change any code here...
This puts the output in a format that can be passed to pretty_table.
"""
function output(nasheqm, km_mkts, numfirms)
    vars = vcat([:aggQ], [Symbol("q$i") for i = 1:numfirms])
    nashoutput = vcat([nasheqm[:Q]], [nasheqm[:firmqtys][i] for i = 1:numfirms])
    kmoutput = [Sb.mean(km_mkts[:, col]) for col in vars]
    ([:Var, :Cournot, :KM], hcat(vars, nashoutput, kmoutput))
end


function main()
    pv = paramvalues(numfirms = 4)
    pv[:costs][1] = 10.0
    println("paramvalues: $(pv)")
    nasheq = cournot_nash(pv)
    km_mkts_df = run_reps(10, pv)
    headers, data = output(nasheq, km_mkts, pv[:numfirms])
    pt = Pt.pretty_table(data, headers)
    println(pt)
    return nasheq, km_mkts_df
end

end ## end module
