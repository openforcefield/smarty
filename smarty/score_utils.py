import numpy

def load_trajectory( trajFile):
    """Load data from a specified smarty trajectory .csv file and return a summary.

    Note that any SMARTS patterns which do not match anything are ignored in the resulting summary.

    Parameters
    ----------
    
        trajFile (str) : filename to read from

    Returns
    -------
        timeseries (dict) : status by iteration number
            Dictionary, keyed by iteration, storing the state at each iteration
            Subsequent keys are by reference atom types, i.e. timeseries[1]['HO']
            gives data at step 1 on what (if anything) matches 'HO'. Subsequent
            keys are 'smarts', 'matches', 'molecules', 'atomsmatched', 'index' (serial # 
            of match), `ParNum` (parent number), `ParentParNum` (parent of parent)
            `denominator` (number of possible matches of this type), `fraction`
            (fraction of this type matched).


    TO DO: Shift to pandas data frame 
    """

    file = open(trajFile, 'r')
    header = file.readline()
    tmp = header.split(',')
    
    # If the number if headers is not as expected, this is a different version and we can't parse
    if len(tmp) != 10:
        raise Exception("Number of headers in trajectory not as expected; can't parse.")

    # Read rest of file
    text = file.readlines()

    # Initialize storage
    timeseries = {}

    # How many iterations are we looking at?
    lastline_elements = text[-1].split(',')
    max_its = int( lastline_elements[0])
    # Number of lines
    max_lines = len(text)

    
    # Process file    
    linenr = 0
    while linenr < max_lines:
        line_elements = text[linenr].split(',')
        iteration = int(line_elements[0])

        # Pull elements from line and store
        if not timeseries.has_key( iteration): timeseries[iteration] = {} 
        reftype = line_elements[5]

        if not reftype=="'NONE'":
            timeseries[iteration][reftype]={}
            timeseries[iteration][reftype]['smarts'] = line_elements[2]
            timeseries[iteration][reftype]['index'] = int(line_elements[1])
            timeseries[iteration][reftype]['ParNum'] = int(line_elements[3])
            timeseries[iteration][reftype]['ParentParNum'] = int(line_elements[4])
            timeseries[iteration][reftype]['matches'] = int(line_elements[6])
            timeseries[iteration][reftype]['molecules'] = int(line_elements[7])
            timeseries[iteration][reftype]['atomsmatched'] = int(line_elements[8])
            timeseries[iteration][reftype]['denominator'] = int(line_elements[9])
            try:
                timeseries[iteration][reftype]['fraction'] = timeseries[iteration][reftype]['atomsmatched']/float(timeseries[iteration][reftype]['denominator'])
            except ZeroDivisionError:
                print("At iteration %s, found %s matched atoms and a denominator of %s for reftype %s..." % (iteration, timeseries[iteration][reftype]['atomsmatched'], timeseries[iteration][reftype]['denominator'], reftype))
                raise 

        linenr+=1


    return timeseries


def scores_vs_time(timeseries):
    """Process a timeseries as read by load_trajectory and return the fraction of each reference atom type found at each time.


    Parameters
    ----------
    trajectory : dict
        Trajectory information as output by load_trajectory

    Returns
    -------
    time_fractions : dict
        Dictionary of NumPy arrays, keyed by reference type. 
        The full score across all types is under `all`.

    """

    # How many iterations are present?
    max_its = numpy.max(timeseries.keys())

    # Retrieve keys of all reference types
    reftypes = set()
    for it in timeseries.keys():
        for reftype in timeseries[it].keys():
            if reftype not in reftypes:
                 reftypes.add(reftype)

    # Allocate storage
    time_fractions = {}
    time_fractions['all'] = numpy.zeros( max_its, float)
    for reftype in reftypes:
        time_fractions[reftype] = numpy.zeros( max_its, float)

    
    # Update with data
    for it in range(max_its):
        # Update reference types occuring at this iteration
        denom = 0
        numer = 0
        for reftype in reftypes:
            if timeseries[it].has_key(reftype):
                try:
                    time_fractions[reftype][it] = timeseries[it][reftype]['fraction']
                except KeyError:
                    print("Can't find key set %s, %s, %s for timeseries." % (it, reftype, 'fraction'))
                    print("Available keys:", timeseries[it][reftype].keys())
                #print("%s, %s - %s" % (it, reftype, timeseries[it][reftype]['denominator']))
                denom += timeseries[it][reftype]['matches']
                numer += timeseries[it][reftype]['atomsmatched']
                
            # Any reference type which does not appear at this time point has zero matches so we just leave the value at zero
        
        # Handle 'all' case last
        time_fractions['all'][it] = numer/float(denom)

    return time_fractions
