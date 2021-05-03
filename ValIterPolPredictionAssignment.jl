## For the problem given in Example 4.1 of Sutton and Barto, Reinforcement Learning
## implements value iteration algorithm for estimating value of a policy, and

module ValIterPolPredictionAssignment

import DataFrames; const Df = DataFrames
import CSV

"""
** You don't need to change this **
Returns a dictionary containing parameter values
"""
function paramvalues()
    allstates = Set([(r, c) for r=1:4 for c=1:4])
    terminal_states = Set([(1,1), (4,4)])
    nonterminal_states = setdiff(allstates, terminal_states)
    actions = [:up, :down, :left, :right]

    Dict(
        :theta => 0.001,
        :gamma => 1,
        :allstates => allstates,
        :terminal_states => terminal_states,
        :nonterminal_states => nonterminal_states,
        :initv => Dict(s => 0.0 for s in allstates),

        ## since the possible actions are the same in every state we can define
        ## set of actions independent of the states.
        :actions => actions,
        :poss_rewards => [-1],

        ## Parameter :initpol is only required for policy iteration for finding
        ## optimal policies. We are finding all optimal actions for each state
        :initpol => Dict(s => Dict(a => (a == :up) ? 1.0 : 0.0 for a in actions)
                              for s in nonterminal_states)
    )
end


"""
** You don't have to change this function **

Translate action (:up, :down, :left, :right) to [rowchange, colchange]
This function does not worry about feasibility. It just translates it
into [rowchange, colchange] value that can be worked with on the grid.
"""
function action_to_movement(action :: Symbol)
    Dict(
        :up => [-1, 0],
        :down => [1, 0],
        :left => [0, -1],
        :right => [0, 1]
    )[action]
end


"""
Calculate [row, col] of the next state given current state `st`,
action `a` and the set of feasible states `allstates`.

State `st` is an array of two integers [row, col]
Action `a` is either :up, :down, :left, or :right.

Return value: Array of two integers [row, col] denoting the next state.
              If the action results in an infeasible state, the agent remains in the
              same state `st`

            Some examples for the specific Gridworld problem we are studying.
            Example 1: st = [1,2], a = :up should return [1,2]
            Example 2: st = [1,2], a = :down should return [2,2]
            Example 3: st = [2,4], a = :right should return [2,4]
"""
function nextst(st, a, allstates)
    ## use action_to_movement function to convert a to [rowchange, colchange]
    ## use st, [rowchange, colchange] and allstates to determine next state.
end


"""
This is the p(s', r'|s, a) function in the algorithm.
It gives the probability of that next state is sprime and reward is rprime
given current state s and action a.

In the Gridworld problem reward is
-1 for all transitions. Hence, for any rprime value other than -1
this function should return 0.

You also have to check that nextst(s, a, allstates) = sprime. If
it is, and rprime = -1, then probability is 1, or it is zero.

Return value: Float64.

Some Examples:
    prob([1, 2], -1, [1, 3], :left, allstates) = 1
    prob([1, 2], -1, [1, 3], :down, allstates) = 0
    prob([1, 2], 0, [1, 3], :left, allstates) = 0
    prob([2, 3], -1, [1, 3], :down, allstates) = 1
"""
function prob(sprime, rprime, s, a, allstates)
    ## your code here.
end


"""
Use a dictionary to denote the policy function.

This function should return a dictionary with possible nonterminal
states (`nontermstates`) as keys. For each key (which is a nonterminal state),
the value should be another dictionary with actions (from passed in argument possactions)
as keys, and for each action key, its probability as value.

Since in this exercise we are evaluating an equiprobable random policy,
the probability of each action will be 1/(number of actions) in possactions.

Partial example (for the Gridworld problem):
    equiprob_rndpolicy(nontermstates, possactions) = Dict(
        [1, 2] => Dict(:up => 1/4, :down => 1/4, :left => 1/4, :right = 1/4),
        [1, 3] => Dict(:up => 1/4, :down => 1/4, :left => 1/4, :right = 1/4)
        [1, 4] => Dict(:up => 1/4, :down => 1/4, :left => 1/4, :right = 1/4)
        .
        .
        .
        [4, 3] => Dict(:up => 1/4, :down => 1/4, :left => 1/4, :right = 1/4)
    )
"""
function equiprob_rndpolicy(nontermstates, possactions)
    ## your code here.
end


"""
Method 1: Calculate expected value for a given current state `currst` and `action`,
given the current estimate of the value function `valest`.

Note that this is the expected value for a specific *action* in a specific
state (the term in the square brackets in the Value Iteration Algorithm for
Policy Evaluation in your notes).

See the paramvalues for what allstates and gamma look like.
        allrewards = [-1]

Below, we will calculate expected value for a given *policy* in a specific state.
"""
function expvalue(currst, action, allstates, allrewards, valest, gamma)
    ## Your code here
end



"""
Method 2: Calculate expected value for a given current state `currst` and `policy`,
given the current estimate of the value function `valest`.

To find what the `policy` parameter looks like see the output of equiprob_rndpolicy function.

See the paramvalues for what allstates and gamma look like.
        allrewards = [-1],
        allactions is the same as actions in paramvalues

This calculates the right hand side of the equation in the Value Iteration
Algorithm for Policy Evaluation in your notes. You will find it convenient to use
Method 1 you described above to calculate this expected value.
"""
function expvalue(currst, allactions, policy,
                  allstates, allrewards, valest, gamma)
    ## Write your code here...
end


"""
Calculate max-norm for the value function. Used to determine convergence.
See step 2 of the Value Iteration Algorithm for Policy Evaluation

Inputs: vest0 and vest1 are two value estimates.
A value estimate is a Dictionary with all states as keys and the value estimate (a Float 64)
as values. See `initv` in paramvalues to get an idea about the structure.
"""
function calc_maxnorm(vest0, vest1)
    ## Your code here.
end


"""
This uses the above functions to implement the value iteration algorithm.
Inputs are params (see the output of paramvalues function)
and policy (see the output of equiprob_rndpolicy function)

The return is a value estimate after the algorithm has converged
(see the `:initv` in paramvalues for an idea of what
the return type of this function.)
"""
function valiter(params, policy)
    ## Make the appropriate initializations and code the steps for the
    ## algorithm.
end


function main()
    params = paramvalues()
    ##-------------------------------------------------------------------------
    ## Value of random policy
    gamma = params[:gamma]
    policy = equiprob_rndpolicy(params[:nonterminal_states], params[:actions])
    res = valiter(params, policy)
    # ## when type is specified as Union{Float64, Missing} instead of only Float64
    # ## the empty array shows `missing` for unassigned values. I find this
    # ## better than seeing some arbitrary number when the grid is printed out.
    grid = Array{Union{Float64, Missing}}(undef, 4, 4)
    for st in keys(res)
        grid[st[1], st[2]] =  round(res[st], digits = 1)
    end
    grid
end

end ## end module
