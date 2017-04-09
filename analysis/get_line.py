from numpy import *
from linetools.lists.linelist import LineList


def get_line(input_line):
    """
    output = get_line(input_line)

    :param input_line: linetools-style input, either symbol ('CIV 1548') or wavelength (1548.1)  
    :return: linetools line information. 
    """
    # Load in the linetools line lists.
    line_list = LineList('ISM', verbose=False, closest=True)
    user_line = line_list[input_line]

    # Check that we've actually got a line:
    bad = False
    if not user_line:
        ## BAIL!!
        print('No line information.')
        bad = True

    return user_line