import sys
import numpy.linalg
import time
from datetime import datetime as dt
from contextlib import ContextDecorator
import traceback


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


def unpack(gh):
    """
    unpacks gh list into l m g h type list
    """
    data = []
    k, l = 0, 1
    while k+1 < len(gh):
        for m in range(l+1):
            if m == 0:
                data.append([l, m, gh[k], 0])
                k += 1
            else:
                data.append([l, m, gh[k], gh[k+1]])
                k += 2
    return data


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


#def measurements_methods(meas_data, noave):
#    """
#    get list of unique specs
#    """
#    version_num = get_version()
#    sids = get_specs(meas_data)
#    # list  of measurement records for this specimen
#    #
#    # step through spec by spec
#    #
#    SpecTmps, SpecOuts = [], []
#    for spec in sids:
#        TRM, IRM3D, ATRM, CR = 0, 0, 0, 0
#        expcodes = ""
#        # first collect all data for this specimen and do lab treatments
#        # list  of measurement records for this specimen
#        SpecRecs = get_dictitem(meas_data, 'er_specimen_name', spec, 'T')
#        for rec in SpecRecs:
#            if 'measurement_flag' not in rec.keys():
#                rec['measurement_flag'] = 'g'
#            tmpmeths = rec['magic_method_codes'].split(":")
#            meths = []
#            if "LP-TRM" in tmpmeths:
#                TRM = 1  # catch these suckers here!
#            if "LP-IRM-3D" in tmpmeths:
#                IRM3D = 1  # catch these suckers here!
#            elif "LP-AN-TRM" in tmpmeths:
#                ATRM = 1  # catch these suckers here!
#            elif "LP-CR-TRM" in tmpmeths:
#                CR = 1  # catch these suckers here!
#            #
#            # otherwise write over existing method codes
#            #
#            # find NRM data (LT-NO)
#            #
#            elif float(rec["measurement_temp"]) >= 273. and float(rec["measurement_temp"]) < 323.:
#                # between 0 and 50C is room T measurement
#                if ("measurement_dc_field" not in rec.keys() or
#                        float(rec["measurement_dc_field"]) == 0 or
#                        rec["measurement_dc_field"] == "") and \
#                    ("measurement_ac_field" not in rec.keys() or
#                        float(rec["measurement_ac_field"]) == 0 or
#                        rec["measurement_ac_field"] == ""):
#                    # measurement done in zero field!
#                    if "treatment_temp" not in rec.keys() or \
#                        rec["treatment_temp"].strip() == "" or \
#                        (float(rec["treatment_temp"]) >= 273. and
#                         float(rec["treatment_temp"]) < 298.):
#                        # between 0 and 50C is room T treatment
#                        if "treatment_ac_field" not in rec.keys() or \
#                                rec["treatment_ac_field"] == "" or \
#                                float(rec["treatment_ac_field"]) == 0:
#                            # no AF
#                            if "treatment_dc_field" not in rec.keys() or \
#                                    rec["treatment_dc_field"] == "" or \
#                                    float(rec["treatment_dc_field"]) == 0:  # no IRM!
#                                if "LT-NO" not in meths:
#                                    meths.append("LT-NO")
#                            elif "LT-IRM" not in meths:
#                                meths.append("LT-IRM")  # it's an IRM
#                        #
#                        # find AF/infield/zerofield
#                        #
#                        elif "treatment_dc_field" not in rec.keys() or \
#                                rec["treatment_dc_field"] == "" or \
#                                float(rec["treatment_dc_field"]) == 0:  # no ARM
#                            if "LT-AF-Z" not in meths:
#                                meths.append("LT-AF-Z")
#                        else:  # yes ARM
#                            if "LT-AF-I" not in meths:
#                                meths.append("LT-AF-I")
#                    #
#                    # find Thermal/infield/zerofield
#                    #
#                    elif float(rec["treatment_temp"]) >= 323:  # treatment done at  high T
#                        if TRM == 1:
#                            if "LT-T-I" not in meths:
#                                # TRM - even if zero applied field!
#                                meths.append("LT-T-I")
#                        elif "treatment_dc_field" not in rec.keys() or \
#                                rec["treatment_dc_field"] == "" or \
#                                float(rec["treatment_dc_field"]) == 0.:  # no TRM
#                            if "LT-T-Z" not in meths:
#                                # don't overwrite if part of a TRM experiment!
#                                meths.append("LT-T-Z")
#                        else:  # yes TRM
#                            if "LT-T-I" not in meths:
#                                meths.append("LT-T-I")
#                    #
#                    # find low-T infield,zero field
#                    #
#                    else:  # treatment done at low T
#                        if "treatment_dc_field" not in rec.keys() or \
#                                rec["treatment_dc_field"] == "" or \
#                                float(rec["treatment_dc_field"]) == 0:  # no field
#                            if "LT-LT-Z" not in meths:
#                                meths.append("LT-LT-Z")
#                        else:  # yes field
#                            if "LT-LT-I" not in meths:
#                                meths.append("LT-LT-I")
#                if "measurement_chi_volume" in rec.keys() or "measurement_chi_mass" in rec.keys():
#                    if "LP-X" not in meths:
#                        meths.append("LP-X")
#                elif "measurement_lab_dc_field" in rec.keys() and \
#                        rec["measurement_lab_dc_field"] != 0:
#                    # measurement in presence of dc field and not susceptibility; hysteresis!
#                    if "LP-HYS" not in meths:
#                        hysq = raw_input("Is this a hysteresis experiment? [1]/0")
#                        if hysq == "" or hysq == "1":
#                            meths.append("LP-HYS")
#                        else:
#                            metha = raw_input("Enter the lab protocol code "
#                                              "that best describes this experiment ")
#                            meths.append(metha)
#                methcode = ""
#                for meth in meths:
#                    methcode = methcode+meth.strip()+":"
#                rec["magic_method_codes"] = methcode[:-1]  # assign them back

#            # done with first pass, collect and assign provisional method codes
#            if "measurement_description" not in rec.keys():
#                rec["measurement_description"] = ""
#            rec["er_citation_names"] = "This study"
#            SpecTmps.append(rec)

#    # ready for second pass through, step through specimens,
#    # check whether ptrm, ptrm tail checks, or AARM, etc.
#    for spec in sids:
#        MD, pTRM, IZ, ZI = 0, 0, 0, 0  # these are flags for the lab protocol codes
#        expcodes = ""
#        NewSpecs, SpecMeths = [], []
#        experiment_name, measnum = "", 1
#        if IRM3D == 1:
#            experiment_name = "LP-IRM-3D"
#        if ATRM == 1:
#            experiment_name = "LP-AN-TRM"
#        if CR == 1:
#            experiment_name = "LP-CR"
#        NewSpecs = get_dictitem(SpecTmps, 'er_specimen_name', spec, 'T')
#        # first look for replicate measurements
#        Ninit = len(NewSpecs)
#        if noave != 1:
#            # averages replicate measurements, returns treatment keys that are being used
#            vdata, treatkeys = vspec_magic(NewSpecs)
#            if len(vdata) != len(NewSpecs):
#                print(spec, 'started with ', Ninit, ' ending with ', len(vdata))
#                NewSpecs = vdata
#                print("Averaged replicate measurements")

#        # now look through this specimen's records - try to figure out what experiment it is
#        if len(NewSpecs) > 1:  # more than one meas for this spec - part of an unknown experiment
#            SpecMeths = get_list(NewSpecs, 'magic_method_codes').split(":")
#            # TRM steps, could be TRM acquisition, Shaw or a Thellier experiment or TDS experiment
#            if "LT-T-I" in SpecMeths and experiment_name == "":
#                # collect all the infield steps and look for changes in dc field vector
#                Steps, TI = [], 1
#                for rec in NewSpecs:
#                    methods = get_list(
#                        NewSpecs, 'magic_method_codes').split(":")
#                    if "LT-T-I" in methods:
#                        Steps.append(rec)  # get all infield steps together
#                rec_bak = Steps[0]
#                if "treatment_dc_field_phi" in rec_bak.keys() and \
#                        "treatment_dc_field_theta" in rec_bak.keys():
#                    # at least there is field orientation info
#                    if rec_bak["treatment_dc_field_phi"] != "" and \
#                            rec_bak["treatment_dc_field_theta"] != "":
#                        phi0 = rec_bak["treatment_dc_field_phi"]
#                        theta0 = rec_bak["treatment_dc_field_theta"]
#                        for k in range(1, len(Steps)):
#                            rec = Steps[k]
#                            phi = rec["treatment_dc_field_phi"]
#                            theta = rec["treatment_dc_field_theta"]
#                            if phi != phi0 or theta != theta0:
#                                # if direction changes, is some sort of anisotropy experiment
#                                ANIS = 1
#                if "LT-AF-I" in SpecMeths and "LT-AF-Z" in SpecMeths:  # must be Shaw :(
#                    experiment_name = "LP-PI-TRM:LP-PI-ALT-AFARM"
#                elif TRM == 1:
#                    experiment_name = "LP-TRM"
#            else:
#                TI = 0  # no infield steps at all
#            if "LT-T-Z" in SpecMeths and experiment_name == "":  # thermal demag steps
#                if TI == 0:
#                    experiment_name = "LP-DIR-T"  # just ordinary thermal demag
#                # heart pounding - could be some kind of TRM
#                # normalized paleointensity or LP-TRM-TD experiment
#                elif TRM != 1:
#                    Temps = []
#                    # check through the infield steps -
#                    # if all at same temperature, then must be a demag of a total TRM with checks
#                    for step in Steps:
#                        if step['treatment_temp'] not in Temps:
#                            Temps.append(step['treatment_temp'])
#                    if len(Temps) > 1:
#                        experiment_name = "LP-PI-TRM"  # paleointensity normalized by TRM
#                    else:
#                        # thermal demag of a lab TRM (could be part of a LP-PI-TDS experiment)
#                        experiment_name = "LP-TRM-TD"
#                TZ = 1
#            else:
#                TZ = 0  # no zero field steps at all
#            if "LT-AF-I" in SpecMeths:  # ARM steps
#                Steps = []
#                for rec in NewSpecs:
#                    tmp = rec["magic_method_codes"].split(":")
#                    methods = []
#                    for meth in tmp:
#                        methods.append(meth.strip())
#                    if "LT-AF-I" in methods:
#                        Steps.append(rec)  # get all infield steps together
#                rec_bak = Steps[0]
#                if "treatment_dc_field_phi" in rec_bak.keys() and \
#                        "treatment_dc_field_theta" in rec_bak.keys():
#                    # at least there is field orientation info
#                    if rec_bak["treatment_dc_field_phi"] != "" and \
#                            rec_bak["treatment_dc_field_theta"] != "":
#                        phi0 = rec_bak["treatment_dc_field_phi"]
#                        theta0 = rec_bak["treatment_dc_field_theta"]
#                        ANIS = 0
#                        for k in range(1, len(Steps)):
#                            rec = Steps[k]
#                            phi = rec["treatment_dc_field_phi"]
#                            theta = rec["treatment_dc_field_theta"]
#                            if phi != phi0 or theta != theta0:
#                                # if direction changes, is some sort of anisotropy experiment
#                                ANIS = 1
#                        if ANIS == 1:
#                            experiment_name = "LP-AN-ARM"
#                if experiment_name == "":  # not anisotropy of ARM - acquisition?
#                        field0 = rec_bak["treatment_dc_field"]
#                        ARM = 0
#                        for k in range(1, len(Steps)):
#                            rec = Steps[k]
#                            field = rec["treatment_dc_field"]
#                            if field != field0:
#                                ARM = 1
#                        if ARM == 1:
#                            experiment_name = "LP-ARM"
#                AFI = 1
#            else:
#                AFI = 0  # no ARM steps at all
#            if "LT-AF-Z" in SpecMeths and experiment_name == "":  # AF demag steps
#                if AFI == 0:
#                    experiment_name = "LP-DIR-AF"  # just ordinary AF demag
#                else:  # heart pounding - a pseudothellier?
#                    experiment_name = "LP-PI-ARM"
#                AFZ = 1
#            else:
#                AFZ = 0  # no AF demag at all
#            if "LT-IRM" in SpecMeths:  # IRM
#                Steps = []
#                for rec in NewSpecs:
#                    tmp = rec["magic_method_codes"].split(":")
#                    methods = []
#                    for meth in tmp:
#                        methods.append(meth.strip())
#                    if "LT-IRM" in methods:
#                        Steps.append(rec)  # get all infield steps together
#                rec_bak = Steps[0]
#                if "treatment_dc_field_phi" in rec_bak.keys() and \
#                        "treatment_dc_field_theta" in rec_bak.keys():
#                    # at least there is field orientation info
#                    if rec_bak["treatment_dc_field_phi"] != "" and \
#                            rec_bak["treatment_dc_field_theta"] != "":
#                        phi0 = rec_bak["treatment_dc_field_phi"]
#                        theta0 = rec_bak["treatment_dc_field_theta"]
#                        ANIS = 0
#                        for k in range(1, len(Steps)):
#                            rec = Steps[k]
#                            phi = rec["treatment_dc_field_phi"]
#                            theta = rec["treatment_dc_field_theta"]
#                            if phi != phi0 or theta != theta0:
#                                # if direction changes, is some sort of anisotropy experiment
#                                ANIS = 1
#                        if ANIS == 1:
#                            experiment_name = "LP-AN-IRM"
#                if experiment_name == "":  # not anisotropy of IRM - acquisition?
#                    field0 = rec_bak["treatment_dc_field"]
#                    IRM = 0
#                    for k in range(1, len(Steps)):
#                        rec = Steps[k]
#                        field = rec["treatment_dc_field"]
#                        if field != field0:
#                            IRM = 1
#                    if IRM == 1:
#                        experiment_name = "LP-IRM"
#                IRM = 1
#            else:
#                IRM = 0  # no IRM at all
#            if "LP-X" in SpecMeths:  # susceptibility run
#                Steps = get_dictitem(
#                    NewSpecs, 'magic_method_codes', 'LT-X', 'has')
#                if len(Steps) > 0:
#                    rec_bak = Steps[0]
#                    if "treatment_dc_field_phi" in rec_bak.keys() and \
#                            "treatment_dc_field_theta" in rec_bak.keys():
#                        # at least there is field orientation info
#                        if rec_bak["treatment_dc_field_phi"] != "" and \
#                                rec_bak["treatment_dc_field_theta"] != "":
#                            phi0 = rec_bak["treatment_dc_field_phi"]
#                            theta0 = rec_bak["treatment_dc_field_theta"]
#                            ANIS = 0
#                            for k in range(1, len(Steps)):
#                                rec = Steps[k]
#                                phi = rec["treatment_dc_field_phi"]
#                                theta = rec["treatment_dc_field_theta"]
#                                if phi != phi0 or theta != theta0:
#                                    # if direction changes, is some sort of anisotropy experiment
#                                    ANIS = 1
#                            if ANIS == 1:
#                                experiment_name = "LP-AN-MS"
#            else:
#                CHI = 0  # no susceptibility at all

#            # now need to deal with special thellier experiment problems - first
#            # clear up pTRM checks and  tail checks
#            if experiment_name == "LP-PI-TRM":  # is some sort of thellier experiment
#                rec_bak = NewSpecs[0]
#                tmp = rec_bak["magic_method_codes"].split(":")
#                methbak = []
#                for meth in tmp:
#                    methbak.append(meth.strip())  # previous steps method codes
#                for k in range(1, len(NewSpecs)):
#                    rec = NewSpecs[k]
#                    tmp = rec["magic_method_codes"].split(":")
#                    meths = []
#                    for meth in tmp:
#                        meths.append(meth.strip())  # get this guys method codes

#                    # check if this is a pTRM check
#                    if (float(rec["treatment_temp"]) <
#                            float(rec_bak["treatment_temp"])):  # went backward
#                        if "LT-T-I" in meths and \
#                                "LT-T-Z" in methbak:  # must be a pTRM check after first z

#                            # replace LT-T-I method code with LT-PTRM-I
#                            methcodes = ""
#                            for meth in meths:
#                                if meth != "LT-T-I":
#                                    methcode = methcode+meth.strip()+":"
#                            methcodes = methcodes+"LT-PTRM-I"
#                            meths = methcodes.split(":")
#                            pTRM = 1
#                        # must be pTRM check after first I
#                        elif "LT-T-Z" in meths and "LT-T-I" in methbak:
#                            # replace LT-T-Z method code with LT-PTRM-Z
#                            methcodes = ""
#                            for meth in meths:
#                                if meth != "LT-T-Z":
#                                    methcode = methcode+meth+":"
#                            methcodes = methcodes+"LT-PTRM-Z"
#                            meths = methcodes.split(":")
#                            pTRM = 1
#                    methcodes = ""
#                    for meth in meths:
#                        methcodes = methcodes+meth.strip()+":"
#                    # attach new method code
#                    rec["magic_method_codes"] = methcodes[:-1]
#                    rec_bak = rec  # next previous record
#                    tmp = rec_bak["magic_method_codes"].split(":")
#                    methbak = []
#                    for meth in tmp:
#                        # previous steps method codes
#                        methbak.append(meth.strip())

#                # done with assigning pTRM checks.  data should be "fixed" in
#                # NewSpecs now let's find out which steps are infield zerofield
#                # (IZ) and which are zerofield infield (ZI)
#                rec_bak = NewSpecs[0]
#                tmp = rec_bak["magic_method_codes"].split(":")
#                methbak = []
#                for meth in tmp:
#                    methbak.append(meth.strip())  # previous steps method codes
#                if "LT-NO" not in methbak:  # first measurement is not NRM
#                    if "LT-T-I" in methbak:
#                        IZorZI = "LP-PI-TRM-IZ"  # first pair is IZ
#                    if "LT-T-Z" in methbak:
#                        IZorZI = "LP-PI-TRM-ZI"  # first pair is ZI
#                    if IZorZI not in methbak:
#                        methbak.append(IZorZI)
#                    methcode = ""
#                    for meth in methbak:
#                        methcode = methcode+meth+":"
#                    # fix first heating step when no NRM
#                    NewSpecs[0]["magic_method_codes"] = methcode[:-1]
#                else:
#                    IZorZI = ""  # first measurement is NRM and not one of a pair
#                for k in range(1, len(NewSpecs)):  # hunt through measurements again
#                    rec = NewSpecs[k]
#                    tmp = rec["magic_method_codes"].split(":")
#                    meths = []
#                    for meth in tmp:
#                        meths.append(meth.strip())  # get this guys method codes

#                    # check if this start a new temperature step of a infield/zerofield pair
#                    if (float(rec["treatment_temp"]) >
#                            float(rec_bak["treatment_temp"])) and \
#                            "LT-PTRM-I" not in methbak:  # new pair?
#                        if "LT-T-I" in meths:  # infield of this pair
#                                IZorZI = "LP-PI-TRM-IZ"
#                                IZ = 1  # at least one IZ pair
#                        elif "LT-T-Z" in meths:  # zerofield
#                                IZorZI = "LP-PI-TRM-ZI"
#                                ZI = 1  # at least one ZI pair
#                    # new pair after out of sequence PTRM check?
#                    elif (float(rec["treatment_temp"]) >
#                          float(rec_bak["treatment_temp"])) and \
#                            "LT-PTRM-I" in methbak and IZorZI != "LP-PI-TRM-ZI":
#                        if "LT-T-I" in meths:  # infield of this pair
#                                IZorZI = "LP-PI-TRM-IZ"
#                                IZ = 1  # at least one IZ pair
#                        elif "LT-T-Z" in meths:  # zerofield
#                                IZorZI = "LP-PI-TRM-ZI"
#                                ZI = 1  # at least one ZI pair
#                    # stayed same temp
#                    if float(rec["treatment_temp"]) == float(rec_bak["treatment_temp"]):
#                        if "LT-T-Z" in meths and \
#                            "LT-T-I" in methbak and \
#                                IZorZI == "LP-PI-TRM-ZI":  # must be a tail check
#                            # replace LT-T-Z method code with LT-PTRM-MD
#                            methcodes = ""
#                            for meth in meths:
#                                if meth != "LT-T-Z":
#                                    methcode = methcode+meth+":"
#                            methcodes = methcodes+"LT-PTRM-MD"
#                            meths = methcodes.split(":")
#                            MD = 1
#                    # fix method codes
#                    if "LT-PTRM-I" not in meths and \
#                        "LT-PTRM-MD" not in meths and \
#                            IZorZI not in meths:
#                        meths.append(IZorZI)
#                    newmeths = []
#                    for meth in meths:
#                        if meth not in newmeths:
#                            newmeths.append(meth)  # try to get uniq set
#                    methcode = ""
#                    for meth in newmeths:
#                        methcode = methcode+meth+":"
#                    rec["magic_method_codes"] = methcode[:-1]
#                    rec_bak = rec  # moving on to next record, making current one the backup
#                    # get last specimen's method codes in a list
#                    methbak = rec_bak["magic_method_codes"].split(":")

#                # done with this specimen's records, now  check if any pTRM checks or MD checks
#                if pTRM == 1:
#                    experiment_name = experiment_name+":LP-PI-ALT-PTRM"
#                if MD == 1:
#                    experiment_name = experiment_name+":LP-PI-BT-MD"
#                if IZ == 1 and ZI == 1:
#                    experiment_name = experiment_name+":LP-PI-BT-IZZI"
#                if IZ == 1 and ZI == 0:
#                    experiment_name = experiment_name+":LP-PI-IZ"  # Aitken method
#                if IZ == 0 and ZI == 1:
#                    experiment_name = experiment_name+":LP-PI-ZI"  # Coe method
#                IZ, ZI, pTRM, MD = 0, 0, 0, 0  # reset these for next specimen

#                # fix the experiment name for all recs for this specimen and
#                # save in SpecOuts
#                for rec in NewSpecs:
#                    # assign an experiment name to all specimen measurements from this specimen
#                    if experiment_name != "":
#                        rec["magic_method_codes"] = rec["magic_method_codes"] + \
#                            ":"+experiment_name
#                    rec["magic_experiment_name"] = spec+":"+experiment_name
#                    rec['measurement_number'] = '%i' % (
#                        measnum)  # assign measurement numbers
#                    measnum += 1
#                    SpecOuts.append(rec)
#            elif experiment_name == "LP-PI-TRM:LP-PI-ALT-AFARM":  # is a Shaw experiment!
#                ARM, TRM = 0, 0
#                # fix the experiment name for all recs for this specimen and
#                # save in SpecOuts
#                for rec in NewSpecs:
#                    # assign an experiment name to all specimen measurements
#                    # from this specimen
#                    # make the second ARM in Shaw experiments LT-AF-I-2, stick
#                    # in the AF of ARM and TRM codes
#                    meths = rec["magic_method_codes"].split(":")
#                    if ARM == 1:
#                        if "LT-AF-I" in meths:
#                            del meths[meths.index("LT-AF-I")]
#                            meths.append("LT-AF-I-2")
#                            ARM = 2
#                        if "LT-AF-Z" in meths and TRM == 0:
#                            meths.append("LP-ARM-AFD")
#                    if TRM == 1 and ARM == 1:
#                        if "LT-AF-Z" in meths:
#                            meths.append("LP-TRM-AFD")
#                    if ARM == 2:
#                        if "LT-AF-Z" in meths:
#                            meths.append("LP-ARM2-AFD")
#                    newcode = ""
#                    for meth in meths:
#                        newcode = newcode+meth+":"
#                    rec["magic_method_codes"] = newcode[:-1]
#                    if "LT-AF-I" in meths:
#                        ARM = 1
#                    if "LT-T-I" in meths:
#                        TRM = 1
#                    rec["magic_method_codes"] = rec["magic_method_codes"] + \
#                        ":"+experiment_name
#                    rec["magic_experiment_name"] = spec+":"+experiment_name
#                    rec['measurement_number'] = '%i' % (
#                        measnum)  # assign measurement numbers
#                    measnum += 1
#                    SpecOuts.append(rec)
#            else:  # not a Thellier-Thellier  or a Shaw experiemnt
#                for rec in NewSpecs:
#                    if experiment_name == "":
#                        rec["magic_method_codes"] = "LT-NO"
#                        rec["magic_experiment_name"] = spec+":LT-NO"
#                        rec['measurement_number'] = '%i' % (
#                            measnum)  # assign measurement numbers
#                        measnum += 1
#                    else:
#                        if experiment_name not in rec['magic_method_codes']:
#                            rec["magic_method_codes"] = rec["magic_method_codes"] + \
#                                ":"+experiment_name
#                            rec["magic_method_codes"] = rec["magic_method_codes"].strip(
#                                ':')
#                        rec['measurement_number'] = '%i' % (
#                            measnum)  # assign measurement numbers
#                        measnum += 1
#                        rec["magic_experiment_name"] = spec+":"+experiment_name
#                    rec["magic_software_packages"] = version_num
#                    SpecOuts.append(rec)
#        else:
#            NewSpecs[0]["magic_experiment_name"] = spec+":" + \
#                NewSpecs[0]['magic_method_codes'].split(':')[0]
#            NewSpecs[0]["magic_software_packages"] = version_num
#            # just copy over the single record as is
#            SpecOuts.append(NewSpecs[0])
#    return SpecOuts


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


# from https://stackoverflow.com/a/3173331
def update_progress(progress):
    print('\r[{0}] {1}%'.format('#'*(progress/10), progress))


class logclr:
    """define colors for output to terminal"""
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'


class Logger(object):
    """Logger
    Redirects all standard output and any unhandled exceptions to file
    'setup.log'. It will also return output to the terminal by default, unless
    the --quiet option was specified during setup.
    """

    def __init__(self, show=None, quiet=False, clr_output=True):
        self.quiet = quiet
        self.clr_output = clr_output
        self.terminal = sys.stdout
        self.quiet_count = 0
        self.log = open("setup.log", "w+")
        self.log.write('\n{:-^80}\n\n'.format('  Started setup at {}  '.format(asctime())))
        self.msg_to_term = True  # send output at start, even if --quiet

    def write(self, message):
        """main write method for redirection of sys.stdout"""
        # first classify the message content to determine whether it should be
        # written to the console in addition to being logged
        if self.msg_to_term or "-W-" in str(message) or "-E-" in str(message):
            # should always evaluate True if quiet==False; if quiet==True,
            # should only evaluate True when printing initial messages during
            # file copy and during warnings/errors
            self.msg_to_term = self.msg_type(message)
        self.log.write(message)

    def msg_type(self, message):
        """classify the message and color output to console"""
        msg_lvl = 0
        # if self.clr_output:
        if '---' in str(message):
            self.terminal.write(logclr.BOLD)
            msg_lvl = 1
        elif '-E-' in str(message):
            self.terminal.write(logclr.FAIL)
            msg_lvl = 2
        elif '-W-' in str(message):
            self.terminal.write(logclr.WARNING)
            msg_lvl = 3
        elif '-I-' in str(message):
            self.terminal.write(logclr.OKGREEN)
            # filter out debug_inp messages if quiet
            if any([x in str(message) for x in ('Running', 'Writing')]):
                msg_lvl = 5
            else:
                msg_lvl = 4
        elif str(message).startswith(' | '):
            self.terminal.write(logclr.OKBLUE)
            msg_lvl = 6
        if self.quiet:
            if msg_lvl < 5:
                self.terminal.write(message)
            else:
                return False
        else:
            self.terminal.write(message)
        self.terminal.write(logclr.ENDC)
        return True

    def flush(self):
        self.terminal.flush()



class loggercontext(ContextDecorator):
    def __init__(self, quiet=False):
        self.quiet = quiet

    def __enter__(self, show=None):
        start_logger(show=self.show, quiet=self.quiet)
        return self

    def __exit__(self, *exc):
        stop_logger(*exc)
        return False


def start_logger(show=None, quiet=False):
    sys.stdout = Logger(show=show, quiet=quiet)


def stop_logger(*exc):
    if not any(exc):
        sys.stdout.log.write(
            '{:-^80}\n'.format('  Setup finished successfully at {}  '.format(asctime())))
        finished = True
    elif exc[0] is SystemExit:
        finished = False
        pass
    else:
        sys.stdout.write('-E- Setup failed! Aborting...\n')
        sys.stdout.log.write('-E- The following error occurred:\n' +
                             ''.join(traceback.format_exception(*exc)))
        finished = False
    sys.stdout.log.close()
    sys.stdout = sys.__stdout__
    if finished:
        print("Setup finished! Full record written to setup.log")
