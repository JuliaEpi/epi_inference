# SIR inference examples

This is a collection of methods for performing inference using contuinuous-time Markov chain versions of the SIR model.
What sets these apart from other methods is that we assume low- or no-noise observations of the system. Typically this is the number of infection or recovery events over a day.

I am in the process of updating lots of old matlab code to Julia to make it faster and easier to share. 

## Examples

See [pmHM inference examples](SIR_examples.md) for examples of these in use

I've implemented three different approaches for calculating or estimating the likelihood. 
The first example below is based on integrating the forward (master) equation. 
This uses a novel method that doesn't require specifiying the stochastic transistion matrix.
The other two use particle filters. 

1. [Forward equation approach](SIR_noQ.jl)
2. [Alive particle filter](SIR_alive.jl)
3. [Importance sampling particle filter](SIR_IS.jl)
