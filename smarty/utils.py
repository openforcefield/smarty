#!/usr/bin/env python

"""
Utility subroutines for SMARTY atom type sampling

"""
#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import sys
import string

from optparse import OptionParser # For parsing of command line arguments

import os
import math
import copy
import re
import numpy
import random

import openeye.oechem
import openeye.oeomega
import openeye.oequacpac

from openeye.oechem import *
from openeye.oeomega import *
from openeye.oequacpac import *

import time
from simtk import unit

#=============================================================================================
# UTILITY ROUTINES
#=============================================================================================

def get_data_filename(relative_path):
    """Get the full path to one of the reference files in testsystems.

    In the source distribution, these files are in ``smarty/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the repex folder).

    """

    from pkg_resources import resource_filename
    fn = resource_filename('smarty', os.path.join('data', relative_path))

    if not os.path.exists(fn):
        raise ValueError("Sorry! %s does not exist. If you just added it, you'll have to re-install" % fn)

    return fn


def parse_odds_file(filename, verbose = False):
    """
    parses files that have the form
    decorator       odds
    if only one column odds will be assumed equally probable

    Parameters
    -----------
    filename: string or file object
    may be an absolute file path, a path relative to the current working directory, a path relative to this module's data subdirectory (for built in decorator files), or an opten file-like object with a readlines() method.

    Returns
    --------
    choices: 2-tuple of the form ( [decorators], [odds] )
    """
    if verbose:
        if isinstance(filename, file):
            print("Attempting to parse file '%s'" % filename.name)
        else:
            print("Attempting to parse file '%s'" % filename)

    # if no file return None
    if filename is None:
        return None

    # if input is a file object
    try:
        input_lines = filename.readlines()
        if verbose: print("Attempting to parse file '%s'" % filename.name)
    except AttributeError:
        if verbose: print("Attempting to parse file '%s'" % filename)
        try:
            ifs = open(filename, 'r')
            input_lines = ifs.readlines()
        except IOError:
            ifs = get_data_filename(filename)
            ifs = open(ifs, 'r')
            input_lines = ifs.readlines()
        except Exception as e:
            raise Exception("%s\nProvided file (%s) could not be parsed" % (str(e), filename))
    except Exception as e:
        msg = str(e) + '\n'
        msg += "Could not read data from file %s" % filename
        raise Exception(msg)

    # close file
    ifs.close()

    decorators = []
    odds = []
    noOdds = False
    for l in input_lines:
        # skip empty lines
        if len(l) == 0:
            continue
        # check for and remove comments
        comment = l.find('%')
        if comment == -1: # no comment
            entry = l.split()
        elif comment > 0: # remove trailing comment
            entry = l[:comment].split()
        else: # whole line is a comment skip
            continue

        # add decorator
        if entry[0] == "''" or entry[0] == '""':
            decorators.append('')
        else:
            decorators.append(entry[0])

        if len(entry) == 2:
            odds.append(float(entry[1]))
        elif len(entry) == 1:
            noOdds = True
        else:
            raise Exception("Error entry (%s) in decorator file '%s' is invalid" % (l, filename))

    if (odds.count(0) == len(odds)) or noOdds:
        odds = None
        #TODO: handle case where 1 line is missing odds entry

    return (decorators, odds)

