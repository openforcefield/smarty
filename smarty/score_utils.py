import numpy
import pandas as pd
import matplotlib
matplotlib.use('pdf')
import pylab as pl

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
            Subsequent keys are by reference types, (i.e. timeseries[1]['HO'])
            and an entry for total if included in the trajectory file at timeseries[1]['total']
            gives data at step 1 on what (if anything) matches 'HO'. Subsequent
            keys are 'smarts', 'matches', 'molecules', 'fractionmatched', 'index' (serial #
            of match), `ParNum` (parameter number/label), `ParentParNum` (parameter number/label of parent)
            `denominator` (number of possible matches of this type), `fraction`
            (fraction of this type matched).

    """
    data = pd.read_csv(trajFile, quotechar="'")
    data_dict = data.to_dict()
    # If the number if headers is not as expected, this is a different version and we can't parse
    if len(data.columns) != 10:
        raise Exception("Number of headers in trajectory not as expected; can't parse.")

    # Initialize storage
    timeseries = {}

    # Number of lines
    max_lines = data.index[-1]

    # How many iterations are we looking at?
    max_its = data.Iteration[max_lines]

    keys = list(data.columns)
    keys.remove('RefType')
    keys.remove('Iteration')

    numerator = data.columns[-2].lower()
    denominator = data.columns[-1].lower()
    # Process file
    for linenr in data.index:
        iteration = data.Iteration[linenr]

        # Pull elements from line and store
        if not iteration in timeseries: timeseries[iteration] = {}
        reftype = data.RefType[linenr]

        if not reftype=="'NONE'":
            timeseries[iteration][reftype]={}
            for k in keys:
                if k in ['ParNum', 'ParentParNum']:
                    timeseries[iteration][reftype][k] = data_dict[k][linenr]
                else:
                    timeseries[iteration][reftype][k.lower()] = data_dict[k][linenr]
            try:
                den = float(timeseries[iteration][reftype][denominator])
                timeseries[iteration][reftype]['fraction'] = timeseries[iteration][reftype][numerator]/den
            except ZeroDivisionError:
                print("At iteration %s, found %s matched atoms and a denominator of %s for reftype %s..." % (iteration, timeseries[iteration][reftype][numerator], timeseries[iteration][reftype][denominator], reftype))
                raise

    return timeseries

def scores_vs_time(timeseries, numerator = 'fractionmatched'
        ):
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
            'all' is from the total list if available or calculated from other references
    """

    # How many iterations are present?
    max_its = numpy.max([k for k in timeseries])

    # Retrieve keys of all reference types
    reftypes = set()
    for it in timeseries:
        for reftype in timeseries[it]:
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
            if reftype in timeseries[it]:
                try:
                    time_fractions[reftype][it] = timeseries[it][reftype]['fraction']
                except KeyError:
                    print("Can't find key set %s, %s, %s for timeseries." % (it, reftype, 'fraction'))
                    print("Available keys:", timeseries[it][reftype])
                denom += timeseries[it][reftype]['denominator']
                numer += timeseries[it][reftype][numerator]

            # Any reference type which does not appear at this time point has zero matches so we just leave the value at zero

        # Handle 'all' case last
        if time_fractions['all'][it] == 0:
            time_fractions['all'][it] = numer/float(denom)

    return time_fractions

def create_plot_file(trajFile, plot_filename, plot_others=False, verbose = False):
    """
    Creates plot to demonstrate performance of smarty or smirky

    trajFile - csv file generated by smarty, smarty_elemental, or smirky
    plot_filename - pdf to save plot file to
    plot_others - if True plots data for all reftypes separately, optional
    """

    data = pd.read_csv(trajFile, quotechar="'")
    numerator = data.columns[-2].lower()

    timeseries = load_trajectory(trajFile)
    time_fractions = scores_vs_time(timeseries, numerator)

    max_score = max(time_fractions['all']) *100.0
    if verbose: print("Maximum score was %.1f %%" % max_score)
    # plot overall score
    pl.plot( time_fractions['all'], 'k-', linewidth = 2.0)

    if plot_others:
        reftypes = [k for k in time_fractions]
        reftypes.remove('all')

        # Plot scors for individual types
        for reftype in reftypes:
            pl.plot(time_fractions[reftype])

        pl.legend(['all']+reftypes, loc='lower right')

    pl.xlabel('Iterations')
    pl.ylabel('Fraction of reference type found')
    pl.ylim(-0.1, 1.1)

    pl.savefig(plot_filename)

