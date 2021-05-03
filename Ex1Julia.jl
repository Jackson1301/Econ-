##-------------------------------------------------------------------------
## Computational Economics Spring 2021
## Assignment 1: Basics of Julia
##
## Instructions:
## 1. Don't change the module name or the filename.
## 2. Read the docstring before each function for what you need to do.
## 3. Don't change function name or signature.
##-------------------------------------------------------------------------
module Ex1Julia

import Random

##-------------------------------------------------------------------------
"""
    freqdist(iter)

Write a function `freqdist` that takes an
iterable and returns a dictionary that shows the frequency count of each
value in the iterable.

For example:

    freqdist([5, 1, 5, 3, 5, 2, 6, 2, 2, 5]) =
        Dict(1 => 1, 2 => 3, 3 => 1, 5 => 4, 6 => 1)

    freqdist("hello world") = Dict('w' => 1, 'h' => 1, 'd' => 1, 'l' => 3,
                    'e' => 1, 'r' => 1, 'o' => 2, ' ' => 1)
"""
function freqdist(iter)
    valuelist = Dict()
    for c in iter
        if haskey(valuelist, c)
            valuelist[c] += 1
        else
            valuelist[c] = 1
        end
    end
    return(valuelist)
end
##-------------------------------------------------------------------------
"""
    has_duplicates(arr)

Write a function called `has_duplicates` that takes an array and returns true if
there is any element that appears more than once. It should not modify the
original array.

Example:

    has_duplicates([10, 15, 5, 8]) ## should return false.
    has_duplicates([8, 9, 4, 9, 10]) ## should return true.
"""
function has_duplicates(arr)
    valuelist = freqdist(arr)
    det = (sum(values(valuelist))-sum(length(valuelist)))
    if det >= 1
        ans = true
    else
        ans = false
    end
    return ans
end
##-------------------------------------------------------------------------
"""
    prob_same_bday(numpeople)

Write a function `prob_same_bday` that takes the
number of people in a group, `numpeople`,
as an argument and uses *simulation* to calculate the probability that any two
people have the same birthday.

Assume that there are only 365 days in the year (don't worry about the leap years),
and anyone in the group has a birthday on one of these 365 days.

To use *simulation*, do the following
(I will use `numpeople` = 25 in this example, but your code should work
for other values of `numpeople`. Your functions should not depend on value of 25,
but should instead use `numpeople`).:
0. Create a variable `samebday` and initialize it to zero.
1. Generate a random sample of 25 numbers (because `numpeople` = 25)
   from 1 through 365. Each number represents a person's birthday.
   Hint: Look up the `rand` function to draw a random sample of integers
   from a given range.
2. If any two numbers in the sample of 25 numbers is the same, increase
   `samebday` by 1.
3. Repeat steps 1 and 2 for 5000 times.
4. The probability is given by `samebday`/5000: your function should return this value.

Example:

    prob_same_bday(5) = 0.0278
    prob_same_bday(25) = 0.5596

"""
import Random
function prob_same_bday(numpeople)
    Random.seed!(1)
    samebday = 0
    for rep = 1:5000
        bdays = []
        for i in 1:numpeople
            push!(bdays, rand(1:365))
        end
        ## I think this is a more accurate way
        # counter = freqdist(bdays)
        # repeats = (sum(values(counter))-sum(length(counter)))
        # samebday += repeats
        if has_duplicates(bdays)
            samebday += 1
        end
    end
    return samebday/5000
end

prob_same_bday(25)
##-------------------------------------------------------------------------
"""
    sort_tuples(arr_tuples, idx)

Write a function sort_tuples(arr_tuples) that takes an array of tuples
and sorts them in the ascending order, by the value in position `idx` in the tuple.

Examples:

    sort_tuples([(2,5), (3,4), (3,6), (1,6)], 1)
    ## returns [(1, 6), (2, 5), (3, 4), (3, 6)]

    sort_tuples([(2,5), (3,4), (3,6), (1,6)], 2)
    ## returns [(3, 4), (2, 5), (3, 6), (1, 6)]

    sort_tuples([(35, 60, 'c'), (31, 96, 'h'), (2, 25, 'd'), (17, 75, 'a')], 3)
    ## returns [(17, 75, 'a'), (35, 60, 'c'), (2, 25, 'd'), (31, 96, 'h')]
"""
function sort_tuples(arr_tuples, idx)
    return sort(arr_tuples, by=arr_tuples->by=arr_tuples[idx])
end
##-------------------------------------------------------------------------
"""
    simpsons(f, a, b, n)

Calculate integral of function `f` (should be a function of one variable)
between `a` and `b` (a < b) using the Simpson's rule.

Parameter `n` is an *even* integer - higher the `n` the better is the approximation
of the integral.

(This exercise is taken from SICP, although the goal there and here is
different.)

Given `f`, `a`, `b`, `n`, the function should return the value:

    h/3 * (y0 + 4*y1 + 2*y2 + 4*y3 + 2*y4 + ... + 4*y_{n-1} + y_{n})

    where, h = (b - a)/n and yk = f(a + k*h)

Example:

    simpsons(x -> x^2, 0, 2, 100) ## 2.666666666666667
    simpsons(x -> x^2, 0, 2, 1000) ## 2.666666666666665
"""
function simpsons(f, a, b, n)
    simp = []
    h = (b - a)/n
    for i in 0:n
        if i == (0 | n)
            push!(simp, f(a + i*h))
        elseif isodd(i)
            push!(simp, 4*f(a + i*h))
        else
            push!(simp, 2*f(a + i*h))
        end
    end
    return(h*sum(simp)/3)
end
##-------------------------------------------------------------------------
"""

    sqroot(num, initguess, tol)

Write a function sqroot(num, initguess, tol)
that takes a positive number `num`, initial guess `initguess`
and returns its square root using the Newton's method given tolerance `tol`.


See the notes on topic 'What is Computational Economics' which describes the
Newton's algorithm for finding zero of a function.
You have to use that algorithm to solve for the value of `x`
for which `f(x) = x^{2} - num = 0`.

"""
function sqroot(num, initguess, tol)
    p = (initguess + num/initguess)/2
    if (abs(p^2 - num) > tol)
        sqroot(num, p, tol)
    else
        return p
    end
end

end ## end module
