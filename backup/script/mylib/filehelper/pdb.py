#! /usr/bin/env python

def beg_model (i):
    return ("MODEL"                ).ljust(6)  +\
           (""                     ).ljust(4)  +\
           ("%d" % i               ).rjust(4)

def end_model():
    return "ENDMDL"

def ter ( serial, resName, chainID, resSeq, iCode ):
    return ("TER"                  ).ljust(6)  +\
           ("%d" % serial          ).rjust(5)  +\
           (""                     ).ljust(6)  +\
           ("%s" % resName         ).rjust(3)  +\
           (""                     ).ljust(1)  +\
           ("%s"%chainID           )           +\
           ("%d" % resSeq          ).rjust(4)  +\
           ("%s" % iCode           )

def atom ( serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge ):
    return ( "ATOM"                ).ljust(6)  +\
           ( "%d" % serial         ).rjust(5)  +\
           ( " "                   )           +\
           ( "%s" % name           ).rjust(4)  +\
           ( "%s" % altLoc         ).ljust(1)  +\
           ( "%s" % resName        ).rjust(3)  +\
           ( " "                   )           +\
           ( "%s" % chainID        )           +\
           ( "%d" % resSeq         ).rjust(4)  +\
           ( "%s" % iCode          ).rjust(1)  +\
           ( ""                    ).ljust(3)  +\
           ( "%8.3f" % float(x)    )           +\
           ( "%8.3f" % float(y)    )           +\
           ( "%8.3f" % float(z)    )           +\
           ( "%6.2f" % occupancy   )           +\
           ( "%6.2f" % tempFactor  )           +\
           ( ""                    ).ljust(10) +\
           ( "%s" %  element       ).rjust(2)  +\
           ( "%s" %  charge        ).rjust(2)

def connect ( me,n1='',n2 = '',n3 = '', n4 = '' ):
    return ( "CONECT"              )           +\
           ( "%s" % str(me)        ).rjust(5)   +\
           ( "%s" % str(n1)        ).rjust(5)   +\
           ( "%s" % str(n2)        ).rjust(5)   +\
           ( "%s" % str(n3)        ).rjust(5)   +\
           ( "%s" % str(n4)        ).rjust(5)
