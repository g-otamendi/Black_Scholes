# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math
from scipy.stats import norm

def black_scholes_option(st, k, r, t, sigma):
    
    d1 = ((math.log( st / k ) ) + ( r + ( ( sigma * sigma ) / 2 ) * t ) ) / ( sigma * math.sqrt( t ) )
    
    d2 = d1 - ( sigma*math.sqrt( t ) )
    
    option_call = st * norm.cdf( d1 ) - k * math.exp( - r * t ) * norm.cdf( d2 )
    
    return (option_call)

def greeks_from_blackscholes(st, k, r, t, sigma):
    
    d1 = ((math.log( st / k ) ) + ( r + ( ( sigma * sigma ) / 2 ) * t ) ) / ( sigma * math.sqrt( t ) )
    
    d2 = d1 - ( sigma*math.sqrt( t ) )
    
    delta = norm.cdf( d1 )
    
    gamma = norm.pdf( d1 )/(st*sigma*math.sqrt(t))
    
    theta = -(st*norm.pdf( d1 )*sigma)/(2*math.sqrt(t))-r*k*math.exp(-r*t)*norm.cdf(d2)
    
    vega = st*math.sqrt(t)*norm.pdf(d1)
    
    rho = k*t*math.exp(-r*t)*norm.cdf(d2)
    
    return (delta, gamma, theta, vega, rho)

st = 100

k = 100

r = 0.05

t = 1

sigma = 0.2

print(greeks_from_blackscholes(st, k, r, t, sigma))

