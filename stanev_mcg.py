""" Pretty-fier - MCG's code assignment 
by Roger Stanev  11/12/2020 v.0
General Instructions:
This is a project intended to give you the opportunity to show off your design, problem solving and testing skills.
Your solution should be as complete and unit tested as possible, we want to see a first simple pass at how you think and how you work.  
The problem is intended to take a couple of hours at most to complete. 

Code Exercise - Write a number prettifier:
Write tested code in Python (no notebooks please) that accepts a numeric type and returns a truncated, "prettified" string version.
The prettified version should include one number after the decimal when the truncated number is not an integer.
It should prettify numbers greater than 6 digits and support millions, billions and trillions. """
# Import libraries
import sys
import math
from decimal import Decimal

# Prettyfier function
def prettyfy(n, precision=1):
    """ Turn numbers into 'pretty' human readable values e.g. 2.5M, 1.1B """

    # K: thousand(s), M: million(s), B: billion(s), T: trillion(s)
    dbase = ['', 'K', 'M', 'B', 'T']
    idx = max(0, min(len(dbase)-1, int(math.floor(0 if n == 0 else math.log10(abs(n)) / 3))))
    result = '{:.{precision}f}'.format(n / 10**(3 * idx), precision=precision)
    
    # Decimal floating point class ref.: https://docs.python.org/3.7/library/decimal.html#
    dec = Decimal(result)

    # Remove the exponent rounding to the nearest integer
    if dec == dec.to_integral():
        # Rounding number to a fixed exponent
        result = dec.quantize(Decimal(1))
    else:
        # Reducing to its simplest form
        result = dec.normalize()

    return '{0}{dx}'.format(result, dx=dbase[idx])

# Main program entry
""" Program takes a user numeric type and 
    returns a truncated 'pretty-fied' string version of the user input 
"""
while True:
    try:
        val = float(input('Please enter a number:'))
        print(prettyfy(val))
    except (ValueError, KeyError) as err:
        print('Invalid input:', err)
        sys.exit(1)