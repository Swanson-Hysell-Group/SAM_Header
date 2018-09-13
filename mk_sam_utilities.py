import sys
import numpy.linalg
import textwrap
import argparse
import time
from datetime import datetime as dt
from contextlib import ContextDecorator


################################################################################
#                               PmagPy functions                               #
################################################################################

def igrf(input_list):
    """
    prints out Declination, Inclination, Intensity data
    from an input list with format:
        [Date, Altitude, Latitude, Longitude]
    Date must be in format
        XXXX.XXXX
    with years and decimals of a year (A.D.)
    """
    x, y, z, f = doigrf(input_list[3] % 360.,
                        input_list[2], input_list[1], input_list[0])
    Dir = cart2dir((x, y, z))
    return Dir


def doigrf(long, lat, alt, date, **kwargs):
    """
    called with doigrf(long,lat,alt,date,**kwargs) calculates the interpolated
    (<2010) or extrapolated (>2010) main field and secular variation
    coefficients and passes these to the Malin and Barraclough routine to
    calculate the IGRF field. dgrf coefficients for 1945 to 2005, igrf for pre
    1945 and post 2010 from http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html

    for dates prior to between 1900 and 1600, this program uses coefficients
    from the GUFM1 model of Jackson et al. 2000 prior to that, it uses either
    arch3k or one of the cals models

    input:
       long  = east longitude in degrees (0 to 360 or -180 to 180)
       lat   = latitude in degrees (-90 to 90)
       alt   = height above mean sea level in km (itype = 1 assumed)
       date  = Required date in years and decimals of a year (A.D.)
    output:
       x     = north component of the magnetic force in nT
       y     = east component of the magnetic force in nT
       z     = downward component of the magnetic force in nT
       f     = total magnetic force in nT

    To check the results you can run the interactive program at the NGDC
    http://www.ngdc.noaa.gov/geomagmodels/IGRFWMM.jsp
    """
    import coefficients as cf
    gh, sv = [], []
    colat = 90.-lat
    # convert to colatitude for MB routine
    if long > 0:
        long = long+360.
    # ensure all positive east longitudes
    itype = 1
    models, igrf12coeffs = cf.get_igrf12()
    if 'mod' in kwargs.keys():
        if kwargs['mod'] == 'arch3k':
            psvmodels, psvcoeffs = cf.get_arch3k()  # use ARCH3k coefficients
        elif kwargs['mod'] == 'cals3k':
            # default: use CALS3K_4b coefficients between -1000,1900
            psvmodels, psvcoeffs = cf.get_cals3k()
        else:
            psvmodels, psvcoeffs = cf.get_cals10k()  # use prior to -1000
    # use geodetic coordinates
    if 'models' in kwargs:
        if 'mod' in kwargs.keys():
            return psvmodels, psvcoeffs
        else:
            return models, igrf12coeffs
    if date < -8000:
        print('too old')
        sys.exit()
    if date < -1000:
        model = date-date % 50
        gh = psvcoeffs[psvcoeffs.index(model)]
        sv = (psvcoeffs[psvmodels.index(model+50)]-gh)/50.
        x, y, z, f = magsyn(gh, sv, model, date, itype, alt, colat, long)
    elif date < 1900:
        model = date-date % 50
        gh = psvcoeffs[psvmodels.index(model)]
        if model+50 < 1900:
            sv = (psvcoeffs[psvmodels.index(model+50)]-gh)/50.
        else:
            field2 = igrf12coeffs[models.index(1900)][0:120]
            sv = (field2-gh)/float(1900-model)
        x, y, z, f = magsyn(gh, sv, model, date, itype, alt, colat, long)
    else:
        model = date-date % 5
        if date < 2015:
            gh = igrf12coeffs[models.index(model)]
            sv = (igrf12coeffs[models.index(model+5)]-gh)/5.
            x, y, z, f = magsyn(gh, sv, model, date, itype, alt, colat, long)
        else:
            gh = igrf12coeffs[models.index(2015)]
            sv = igrf12coeffs[models.index(2015.20)]
            x, y, z, f = magsyn(gh, sv, model, date, itype, alt, colat, long)
    if 'coeffs' in kwargs.keys():
        return gh
    else:
        return x, y, z, f


# def unpack(gh):
#     """
#     unpacks gh list into l m g h type list
#     """
#     data = []
#     k, l = 0, 1
#     while k+1 < len(gh):
#         for m in range(l+1):
#             if m == 0:
#                 data.append([l, m, gh[k], 0])
#                 k += 1
#             else:
#                 data.append([l, m, gh[k], gh[k+1]])
#                 k += 2
#     return data


def magsyn(gh, sv, b, date, itype, alt, colat, elong):
    """
    Computes x, y, z, and f for a given date and position, from the
    spherical harmonic coeifficients of the International Geomagnetic
    Reference Field (IGRF).
    From Malin and Barraclough (1981), Computers and Geosciences, V.7, 401-405.

    Input:
          date  = Required date in years and decimals of a year (A.D.)
          itype = 1, if geodetic coordinates are used, 2 if geocentric
          alt   = height above mean sea level in km (if itype = 1)
          alt   = radial distance from the center of the earth (itype = 2)
          colat = colatitude in degrees (0 to 180)
          elong = east longitude in degrees (0 to 360)
                  gh        = main field values for date (calc. in igrf subroutine)
                  sv        = secular variation coefficients (calc. in igrf subroutine)
                  begin = date of dgrf (or igrf) field prior to required date

    Output:
          x     - north component of the magnetic force in nT
          y     - east component of the magnetic force in nT
          z     - downward component of the magnetic force in nT
          f     - total magnetic force in nT

          NB: the coordinate system for x,y, and z is the same as that specified
          by itype.

    Modified 4/9/97 to use DGRFs from 1945 to 1990 IGRF
    Modified 10/13/06 to use  1995 DGRF, 2005 IGRF and sv coefficient
    for extrapolation beyond 2005. Coefficients from Barton et al. PEPI, 97: 23-26
    (1996), via web site for NOAA, World Data Center A. Modified to use
    degree and
    order 10 as per notes in Malin and Barraclough (1981).
    coefficients for DGRF 1995 and IGRF 2005 are from
    http://nssdcftp.gsfc.nasa.gov/models/geomagnetic/igrf/fortran_code/
    igrf subroutine calculates
    the proper main field and secular variation coefficients (interpolated between
    dgrf values or extrapolated from 1995 sv values as appropriate).
    """
    p = numpy.zeros((66), 'f')
    q = numpy.zeros((66), 'f')
    cl = numpy.zeros((10), 'f')
    sl = numpy.zeros((10), 'f')
    begin = b
    t = date - begin
    r = alt
    one = colat*0.0174532925
    ct = numpy.cos(one)
    st = numpy.sin(one)
    one = elong*0.0174532925
    cl[0] = numpy.cos(one)
    sl[0] = numpy.sin(one)
    x, y, z = 0.0, 0.0, 0.0
    cd, sd = 1.0, 0.0
    l, ll, m, n = 1, 0, 1, 0
    if itype != 2:
        # if required, convert from geodectic to geocentric
        a2 = 40680925.0
        b2 = 40408585.0
        one = a2 * st * st
        two = b2 * ct * ct
        three = one + two
        rho = numpy.sqrt(three)
        r = numpy.sqrt(alt*(alt+2.0*rho) + (a2*one+b2*two)/three)
        cd = (alt + rho) / r
        sd = (a2 - b2) / rho * ct * st / r
        one = ct
        ct = ct*cd - st*sd
        st = st*cd + one*sd
    ratio = 6371.2 / r
    rr = ratio * ratio

    # compute Schmidt quasi-normal coefficients p and x(=q)
    p[0] = 1.0
    p[2] = st
    q[0] = 0.0
    q[2] = ct
    for k in range(1, 66):
        if n < m:  # else go to 2
            m = 0
            n = n + 1
            rr = rr * ratio
            fn = n
            gn = n - 1
        # 2
        fm = m
        if k != 2:  # else go to 4
            if m == n:   # else go to 3
                one = numpy.sqrt(1.0 - 0.5/fm)
                j = k - n - 1
                p[k] = one * st * p[j]
                q[k] = one * (st*q[j] + ct*p[j])
                cl[m-1] = cl[m-2]*cl[0] - sl[m-2]*sl[0]
                sl[m-1] = sl[m-2]*cl[0] + cl[m-2]*sl[0]
            else:  # 3
                gm = m * m
                one = numpy.sqrt(fn*fn - gm)
                two = numpy.sqrt(gn*gn - gm) / one
                three = (fn + gn) / one
                i = k - n
                j = i - n + 1
                p[k] = three*ct*p[i] - two*p[j]
                q[k] = three*(ct*q[i] - st*p[i]) - two*q[j]

        # synthesize x, y, and z in geocentric coordinates.
        # 4
        one = (gh[l-1] + sv[ll+l-1]*t)*rr
        if m != 0:  # else go to 7
            two = (gh[l] + sv[ll+l]*t)*rr
            three = one*cl[m-1] + two*sl[m-1]
            x = x + three*q[k]
            z = z - (fn + 1.0)*three*p[k]
            if st != 0.0:  # else go to 5
                y = y + (one*sl[m-1] - two*cl[m-1])*fm*p[k]/st
            else:
                # 5
                y = y + (one*sl[m-1] - two*cl[m-1])*q[k]*ct
            l = l + 2
        else:  # 7
            x = x + one*q[k]
            z = z - (fn + 1.0)*one*p[k]
            l = l + 1
        m = m + 1
    # convert to coordinate system specified by itype
    one = x
    x = x*cd + z*sd
    z = z*cd - one*sd
    f = numpy.sqrt(x*x + y*y + z*z)
    return x, y, z, f


def cart2dir(cart):
    """
    converts a direction to cartesian coordinates
    """
    cart = numpy.array(cart)
    rad = numpy.pi/180.  # constant to convert degrees to radians
    if len(cart.shape) > 1:
        Xs, Ys, Zs = cart[:, 0], cart[:, 1], cart[:, 2]
    else:  # single vector
        Xs, Ys, Zs = cart[0], cart[1], cart[2]
    Rs = numpy.sqrt(Xs**2+Ys**2+Zs**2)  # calculate resultant vector length
    # calculate declination taking care of correct quadrants (arctan2) and making modulo 360.
    Decs = (numpy.arctan2(Ys, Xs)/rad) % 360.
    try:
        # calculate inclination (converting to degrees) #
        Incs = numpy.arcsin(Zs/Rs)/rad
    except:
        print('trouble in cart2dir')  # most likely division by zero somewhere
        return numpy.zeros(3)
    # return the directions list
    return numpy.array([Decs, Incs, Rs]).transpose()


def sundec(sundata):
    """
    returns the declination for a given set of suncompass data

    INPUT:
      sundata={'date':'yyyy:mm:dd:hr:min','delta_u':DU,'lat':LAT,'lon':LON,'shadow_angle':SHADAZ}
      where:
         DU is the hours to subtract from local time to get Greenwich Mean Time
         LAT,LON are the site latitude,longitude (negative for south and west respectively)
         SHADAZ is the shadow angle of the desired direction with respect to the sun.

    OUTPUT:
      the declination of the desired direction wrt true north.
    """
    rad = numpy.pi/180.
    iday = 0
    timedate = sundata["date"]
    timedate = timedate.split(":")
    year = int(timedate[0])
    mon = int(timedate[1])
    day = int(timedate[2])
    hours = float(timedate[3])
    min = float(timedate[4])
    du = int(sundata["delta_u"])
    hrs = hours-du
    if hrs > 24:
        day += 1
        hrs = hrs-24
    if hrs < 0:
        day = day-1
        hrs = hrs+24
    julian_day = julian(mon, day, year)
    utd = (hrs+min/60.)/24.
    greenwich_hour_angle, delta = gha(julian_day, utd)
    H = greenwich_hour_angle+float(sundata["lon"])
    if H > 360:
        H = H-360
    lat = float(sundata["lat"])
    if H > 90 and H < 270:
        lat = -lat
    # now do spherical trig to get azimuth to sun
    lat = (lat)*rad
    delta = (delta)*rad
    H = H*rad
    ctheta = numpy.sin(lat)*numpy.sin(delta)+numpy.cos(lat) * \
        numpy.cos(delta)*numpy.cos(H)
    theta = numpy.arccos(ctheta)
    beta = numpy.cos(delta)*numpy.sin(H)/numpy.sin(theta)
    # check which beta
    beta = numpy.arcsin(beta)/rad
    if delta < lat:
        beta = 180-beta
    sunaz = 180-beta
    suncor = (sunaz+float(sundata["shadow_angle"])) % 360.  # mod 360
    return suncor


def gha(julian_day, f):
    """
    returns greenwich hour angle
    """
    rad = numpy.pi/180.
    d = julian_day-2451545.0+f
    L = 280.460 + 0.9856474*d
    g = 357.528 + 0.9856003*d
    L = L % 360.
    g = g % 360.
    # ecliptic longitude
    lamb = L+1.915*numpy.sin(g*rad)+.02*numpy.sin(2*g*rad)
    # obliquity of ecliptic
    epsilon = 23.439 - 0.0000004*d
    # right ascension (in same quadrant as lambda)
    t = (numpy.tan((epsilon*rad)/2))**2
    r = 1/rad
    rl = lamb*rad
    alpha = lamb-r*t*numpy.sin(2*rl)+(r/2)*t*t*numpy.sin(4*rl)
    # alpha=mod(alpha,360.0)
    # declination
    delta = numpy.sin(epsilon*rad)*numpy.sin(lamb*rad)
    delta = numpy.arcsin(delta)/rad
    # equation of time
    eqt = (L-alpha)
    #
    utm = f*24*60
    H = utm/4+eqt+180
    H = H % 360.0
    return H, delta


def julian(mon, day, year):
    """
    returns julian day
    """
    ig = 15+31*(10+12*1582)
    if year == 0:
        print("Julian no can do")
        return
    if year < 0:
        year = year+1
    if mon > 2:
        julian_year = year
        julian_month = mon+1
    else:
        julian_year = year-1
        julian_month = mon+13
    j1 = int(365.25*julian_year)
    j2 = int(30.6001*julian_month)
    j3 = day+1720995
    julian_day = j1+j2+j3
    if day+31*(mon+12*year) >= ig:
        jadj = int(0.01*julian_year)
        julian_day = julian_day+2-jadj+int(0.25*jadj)
    return julian_day


def to_year_fraction(date):
    """authored by ninjagecko on stackoverflow:
    http://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    """
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction




################################################################################
#                            Command Line Interface                            #
################################################################################

################
#  Docstrings  #
################

def get_prog_desc():
    # define docstrings to print with -h or --help options
    prog_desc = textwrap.dedent("""\
            Using a formatted CSV, creates and writes a .sam header file and a
            set of sample files. The program can be run on multiple files
            provided either as an explicit sequence of names or as a glob
            pattern (use --help to view examples).
            """)
    return prog_desc


def get_prog_epilog():
    prog_epilog = textwrap.dedent("""\

    Examples
    --------
    Create header and sample files from a single .csv for site `P1` containing
    samples `1a` through `6a`:

        $ mk_sam_file.py P1.csv

        Output ('.' = current directory):
            ./P1.csv(new) ./P1.sam  ./P1.inp  ./P1-1a  ./P1-2a
            ./P1-3a       ./P1-4a   ./P1-5a   ./P1-6a

    Same as above, but write files to another directory:

        $ mk_sam_file.py P1.csv --dirout P1_files

        Output:
            ./P1_files/P1.csv(new)  ./P1_files/P1.sam  ./P1_files/P1.inp
            ./P1_files/P1-1         ./P1_files/P1-2    [...]

    Run script on multiple .csv files:

        $ mk_sam_file.py P1.csv P2.csv P3.csv
        --OR--
        $ mk_sam_file.py P[1-3].csv

        Output:
            ./P1.csv(new)  ./P1.sam  ./P1.inp  ./P1-1a  ./P1-2a  [...]
            ./P2.csv(new)  ./P2.sam  ./P2.inp  ./P2-1a  ./P2-2a  [...]
            ./P3.csv(new)  ./P3.sam  ./P3.inp  ./P3-1a  ./P3-2a  [...]

    Other options --all and --auto-dirs:

        $ mk_sam_file.py --all --auto-dirs

        Output:
            ./P1/P1.csv(new)  ./P1/P1.sam  ./P1/P1.inp  ./P1/P1-1a  [...]
            ./P2/P2.csv(new)  ./P2/P2.sam  ./P2/P2.inp  ./P2/P2-1a  [...]
            ...
            ./Z15/Z15.csv(new)  ./Z15/Z15.sam  ./Z15/Z15.inp  ./Z15/Z15.1a  [...]
            ./CF10/CF10.csv(new)  ./CF10/CF10.sam  ./CF10/CF10.inp  [...]
            etc.

    """)
    return prog_epilog


##################################
#  Command line argument parser  #
##################################

def get_parser(with_examples=False):
    # get docstrings
    prog_desc = get_prog_desc()
    if with_examples:
        prog_epilog = get_prog_epilog()
    else:
        prog_epilog = None
    # set up argument parser for command line use
    parser = argparse.ArgumentParser(prog="mk_sam_file.py", add_help=False,
                                     description=prog_desc,
                                     epilog=prog_epilog,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    config_main = parser.add_argument_group(title="Positional arguments")
    config_main.add_argument('csv_file', nargs='*',
                             help=""".csv file(s). Required unless the --all
                             option is given""")
    config_opts = parser.add_argument_group(title="Optional arguments")
    config_opts.add_argument('-h', '--help', action='help',
                             help="""Show this help message and exit. Use long
                             form (--help) to include examples.""")
    config_opts.add_argument('-a', '--all', action='store_true',
                             help="""Create .sam header files for all CSV files
                             in the current directory.""")
    config_opts.add_argument('-d', '--dirout', dest='output_directory', metavar='',
                             help="""Output directory. If path does not exist,
                             it will be created. To automatically configure
                             output directories (based on csv name), use
                             --auto-dirs.""")
    config_opts.add_argument('-ad', '--auto-dirs', action='store_true',
                             help="""Write contents to a directory with same
                             name as the csv file.""")
    # CLI options
    cli_opts = parser.add_argument_group(title="Additional options")
    cli_opts.add_argument('-y', '--yes', action='store_true',
                          help="""Do not pause on error notifications.""")
    cli_opts.add_argument('-q', '--quiet', dest='quiet', action='store_true',
                          help="""Show minimal output and write the full output
                          to log file (helpful with a large number of files).
                          Default display method if number of files > 10, unless
                          --verbose option is used.""")
    cli_opts.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                          help="""Display long output (default if number of
                          files <= 10).""")
    cli_opts.add_argument('--no-color', dest='color', action='store_false',
                          help="""Disable color output (enabled by default).""")
    return parser


#######################################################
#  CLI display (controls for color, verbosity, etc.)  #
#######################################################

class tclr:
    """define colors for output to terminal"""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARN = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARN = ''
        self.FAIL = ''
        self.ENDC = ''
        self.BOLD = ''


class Logger(object):
    """Logger
    If quiet==True, replace long output with progress bar. A summary of warnings
    and all exceptions will still be shown. Otherwise, Logger simply relays all
    print statements back to stdout.
    """
    def __init__(self, quiet=False):
        self.quiet = quiet
        self.terminal = sys.stdout
        self.log = open("mk_sam.log", "w+")
        self.log.write('\n{:-^70}\n\n'.format('  Started at {}  '.format(time.asctime())))

    def write(self, message):
        """main write method for redirection of sys.stdout"""
        if not self.quiet:
            self.terminal.write(message)
        else:
            # if applicable, capture and show the final warning message and
            # count for large declination difference
            if ('WARNING: ' in str(message) and
                    'WARNING: declinations' not in str(message)):
                # self.log.write("\033[2K\r")
                message = message.replace('\033[93m', '')
                message = message.replace('\033[0m', '')
                self.log.write(message)
                self.log.write("\n")

    def write_progress(self, n, num, colors = None):
        i = int((n/num)*100)
        width = int((i + 1) / 2)
        bar = "[" + "#" * width + " " * (50 - width) + "]"
        count = f"( {n} / {num} )"
        if n == num and colors is not None:
            count = colors.OKGREEN + f"( {n} / {num} )" + colors.ENDC
        self.terminal.write("\033[1000D" + bar + " "*4 + count)
        self.flush()

    def override_quiet(self):
        """override quiet option for important messages (e.g. errors)"""
        # clear line and carriage return
        self.terminal.write("\033[2K\r")
        self.quiet = False

    def silence(self):
        """restore quiet after override_quiet()"""
        # start on new line
        self.terminal.write("\n")
        self.quiet = True

    def flush(self):
        self.terminal.flush()


class loggercontext(ContextDecorator):
    """loggercontext
    A context manager for Logger (simplifies the implementation in the code to a
    `with` statement). This also ensures that sys.stdout is returned to its
    original state.
    """
    def __init__(self, quiet=False):
        self.quiet = quiet

    def __enter__(self):
        sys.stdout = Logger(quiet=self.quiet)
        return sys.stdout

    def __exit__(self, *exc):
        sys.stdout.log.close()
        sys.stdout = sys.__stdout__
        return False
