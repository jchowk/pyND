"""
 Convert a velocity in one velocity definition to an equivalent
 velocity in another definition.

 @param vel {in}{required}{type=double} The velocity, in m/s, to convert.

 @param toveldef {in}{required}{type=string} The desired velocity
 definition.  Must be one of 'TRUE', 'OPTICAL' and 'RADIO'.

 @param fromveldef {in}{required}{type=string} The input velocity
 definition.  Must be one of 'TRUE', 'OPTICAL' and 'RADIO'

 @returns the converted velocity in m/s.

 @version $Id$
"""

from __future__ import print_function

from numpy import *

def veltovel(vel, toveldef, fromveldef):
    n_params = 3
    def _ret():  return None

    # COMPILE_OPT IDL2

    if (fromveldef == toveldef):
        return vel # nothing to do

    result = vel / light_speed
    # convert to true
    _expr = fromveldef
    if _expr == RADIO:
        result = (2.e0 * result - result * result) / (2.e0 - 2.e0 * result + result * result)

    elif _expr == OPTICAL:
        result = (2.e0 * result + result * result) / (2.e0 + 2.e0 * result + result * result)
    else:
        pass # nothing to do


    # convert from true

    _expr = toveldef
    if _expr == RADIO:
        result = 1.e0 - sqrt((1.e0 - result) / (1.e0 + result))

    elif _expr == OPTICAL:
        result = sqrt((1.e0 + result) / (1.e0 - result)) - 1.e0

    else:
        pass # nothing to do

    return result * light_speed
