#!/nsls2/projects/bldev/general64/bin/python-bldev-3 

#############################################################################
# SRWLIB Example#6: Calculating spectral flux of undulator radiation 
# by finite-emittance electron beam collected through a finite aperture
# and power density distribution of this radiation (integrated over all photon energies)
# v 0.02
#############################################################################

from __future__ import print_function  # Python 2.7 compatibility
from srwlib import *


# **********************Auxiliary functions

def AuxSaveTrajData(traj, filePath):
    f = open(filePath, 'w')
    resStr = '#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [rad]'
    if (hasattr(traj, 'arBx')):
        resStr += ', Bx [T]'
    if (hasattr(traj, 'arBy')):
        resStr += ', By [T]'
    if (hasattr(traj, 'arBz')):
        resStr += ', Bz [T]'
    f.write(resStr + '\n')
    ctStep = 0
    if traj.np > 0:
        ctStep = (traj.ctEnd - traj.ctStart) / (traj.np - 1)
    ct = traj.ctStart
    for i in range(traj.np):
        resStr = str(ct) + '\t' + repr(traj.arX[i]) + '\t' + repr(traj.arXp[i]) + '\t' + repr(
            traj.arY[i]) + '\t' + repr(traj.arYp[i]) + '\t' + repr(traj.arZ[i]) + '\t' + repr(traj.arZp[i])
        if (hasattr(traj, 'arBx')):
            resStr += '\t' + repr(traj.arBx[i])
        if (hasattr(traj, 'arBy')):
            resStr += '\t' + repr(traj.arBy[i])
        if (hasattr(traj, 'arBz')):
            resStr += '\t' + repr(traj.arBz[i])
        f.write(resStr + '\n')
        ct += ctStep
    f.close()


# Read data comumns from ASCII file:
def AuxReadInDataColumns(filePath, nCol, strSep):
    f = open(filePath, 'r')
    resCols = []
    for iCol in range(nCol):
        resCols.append([])

    curLine = f.readline()
    while len(curLine) > 0:
        curLineParts = curLine.split(strSep)
        for iCol in range(nCol):
            if (iCol < len(curLineParts)):
                resCols[iCol].append(float(curLineParts[iCol]))
            curLine = f.readline()
    f.close()
    return resCols  # attn: returns lists, not arrays!


def AuxTransmAddSurfHeightProfile(optSlopeErr, heightProfData, dim, ang):
    argHeightProfData = heightProfData[0]
    valHeightProfData = heightProfData[1]
    sinAng = sin(ang)
    npData = len(heightProfData[0])

    # xStep = optSlopeErr.rx/(optSlopeErr.nx - 1)
    # yStep = optSlopeErr.ry/(optSlopeErr.ny - 1)
    # y = optSlopeErr.y - 0.5*optSlopeErr.ry

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart) / (auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart) / (auxMesh.ny - 1)

    y = auxMesh.yStart
    hApprox = 0
    ipStart = 0
    # for iy in range(optSlopeErr.ny):
    for iy in range(auxMesh.ny):
        if ('y' in dim):
            hApprox = 0
            y1 = argHeightProfData[ipStart] * sinAng
            for i in range(ipStart + 1, npData):
                y2 = argHeightProfData[i] * sinAng
                if ((y1 <= y) and (y < y2)):
                    hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / (
                        (argHeightProfData[i] - argHeightProfData[i - 1]) * sinAng)) * (y - y1) + valHeightProfData[
                                  i - 1]
                    # sys.stdout.write(ipStart, i, iy, y1, y, y2, argHeightProfData[i-1], argHeightProfData[i], valHeightProfData[i-1], valHeightProfData[i], hApprox)
                    ipStart = i - 1
                    break
                y1 = y2

        # x = optSlopeErr.x - 0.5*optSlopeErr.rx
        x = auxMesh.xStart

        # for ix in range(optSlopeErr.nx):
        for ix in range(auxMesh.nx):
            if ('x' in dim):
                if (ix == 0): ipStart = 0
                hApprox = 0
                x1 = argHeightProfData[ipStart] * sinAng
                for i in range(ipStart + 1, npData):
                    x2 = argHeightProfData[i] * sinAng
                    if ((x1 <= x) and (x < x2)):
                        hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / (
                            (argHeightProfData[i] - argHeightProfData[i - 1]) * sinAng)) * (x - x1) + valHeightProfData[
                                      i - 1]
                        ipStart = i - 1
                        break
                    x1 = x2
            # ofst = 2*ix + (2*optSlopeErr.nx)*iy
            ofst = 2 * ix + (2 * auxMesh.nx) * iy

            optSlopeErr.arTr[ofst] = 1.  # Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0.  # Optical Path Difference
            if (hApprox != 0):
                optSlopeErr.arTr[ofst + 1] = -2 * sinAng * hApprox  # Optical Path Difference (to check sign!)
                # sys.stdout.write(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep


def drift_ebeam(ebeam, dist):
    ebeam.partStatMom1.z += dist
    ebeam.arStatMom2[0] += ebeam.arStatMom2[1] * dist * 2 + ebeam.arStatMom2[2] * dist * dist
    ebeam.arStatMom2[1] += ebeam.arStatMom2[2] * dist
    ebeam.arStatMom2[3] += ebeam.arStatMom2[4] * dist * 2 + ebeam.arStatMom2[5] * dist * dist
    ebeam.arStatMom2[4] += ebeam.arStatMom2[5] * dist
    return ebeam


def parm_func(GUI_data):
    try:
        beamline = GUI_data['beamline']
    except:
        beamline = 'ES2'  # 'ES1' 'ES2'
        print('We use default value for beamline: %s' % beamline)

    BMmode = 'Norm'  # 'Norm'  'LowDiv'
    bump = True  # False

    strBump = '_bump' if bump else ''

    print('modeling SMI beamline ' + beamline + ' with bump = ' + repr(bump))

    print(
        'Calculating spectral flux of undulator radiation by finite-emittance electron beam collected through a finite aperture')

    # **********************Input Parameters:
    strExDataFolderName = 'smi21crlb'  # example data sub-folder name
    strIntSourSE_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                            'sour_se' + beamline + strBump + '.dat')  # file name for output UR intensity Single Electron
    strIntPropSE_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                            'prop_se' + beamline + strBump + '.dat')  # file name for output intensity propagated Single Electron
    strIntPropME_OutFileName = os.path.join(os.getcwd(), strExDataFolderName,
                                            'prop_me' + beamline + strBump + '.dat')  # file name for output intensity propagated Multi Electron
    # strTrajOutFileName = "res_trj.dat"

    strProfileData_default = os.path.join(os.getcwd(), strExDataFolderName, 'mirror2.dat')
    try:
        if GUI_data['mirror_file']:
            strProfileData = GUI_data['mirror_file']
        else:
            strProfileData = strProfileData_default
    except:
        strProfileData = strProfileData_default
        print('Specified file %s not found. Use the default value %s.' % (
            GUI_data['mirror_file'], strProfileData_default))

    try:
        APE = GUI_data['APE']  # 29.5  # longitudinal position [m] at which UR has to be calculated
    except:
        APE = GUI_data['defaults']['APE']

    try:
        DCM = GUI_data['DCM']
    except:
        DCM = GUI_data['defaults']['DCM']  # 2.44

    try:
        HFM = GUI_data['HFM']
    except:
        HFM = GUI_data['defaults']['HFM']  # 3.023696

    try:
        VFM = GUI_data['VFM']
    except:
        VFM = GUI_data['defaults']['VFM']  # 3.42

    try:
        VM = GUI_data['VM']
    except:
        VM = GUI_data['defaults']['VM']  # 0.7

    try:
        SSA = GUI_data['SSA']
    except:
        SSA = GUI_data['defaults']['SSA']  # 8.0

    try:
        CRL_d = GUI_data['CRL']
    except:
        CRL_d = GUI_data['defaults']['CRL']  # 10.4052

    try:
        Sample = GUI_data['Sample']
    except:
        Sample = GUI_data['defaults']['Sample']  # 3.9

    strProfileDataHFM = os.path.join(os.getcwd(), strExDataFolderName,
                                     'HFM08rms.dat')  # HFM08rms.dat - 464 points 0.5 m
    strProfileDataVFM = os.path.join(os.getcwd(), strExDataFolderName,
                                     'VFM03rms.dat')  # VFM03rms.dat - 464 points 0.5 m
    strProfileDataVM = os.path.join(os.getcwd(), strExDataFolderName, 'VM03rms.dat')  # VM03rmsL.dat - 127 points, 0.5 m
    strProfileDataHKB = os.path.join(os.getcwd(), strExDataFolderName,
                                     'FMsine015.dat')  # FMsine015.dat - 127 points 0.5166 m
    strProfileDataVKB = os.path.join(os.getcwd(), strExDataFolderName, 'VM03rms.dat')  # VM03rms.dat - 500 points 0.5 m

    # *********** Undulator
    numPer = 121.5  # Number of ID Periods (without counting for terminations

    undPer = 0.023  # Period Length [m]
    By = 0.577  # 0.955 #Peak Vertical field [T]

    phBy = 0  # Initial Phase of the Vertical field component

    sBy = -1  # Symmetry of the Vertical field component vs Longitudinal position
    xcID = 0  # Transverse Coordinates of Undulator Center [m]
    ycID = 0
    zcID = 0.6  # 0 #Longitudinal Coordinate of Undulator Center [m]

    sys.stdout.write('   Setup Magnetic Field for Undulator ... ');
    sys.stdout.flush()
    und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1)], undPer, numPer)  # Planar Undulator
    magFldCnt = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]),
                            array('d', [zcID]))  # Container of all Field Elements
    sys.stdout.write('done\n')


    # ***********Electron Beam
    sys.stdout.write('   Setup Electron Beam for Undulator ... ');
    sys.stdout.flush()
    elecBeam = SRWLPartBeam()
    elecBeam.Iavg = 0.5  # Average Current [A]
    elecBeam.partStatMom1.x = 0.  # 200e-06 #0. #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
    elecBeam.partStatMom1.y = 0.  # 30e-06 #0. #-0.00025
    elecBeam.partStatMom1.z = -0.9  # 0. #-0.5*undPer*(numPer + 4) #Initial Longitudinal Coordinate (set before the ID)
    elecBeam.partStatMom1.xp = 0.  # 10e-06 #0. #Initial Relative Transverse Velocities
    elecBeam.partStatMom1.yp = 0.
    elecBeam.partStatMom1.gamma = 3. / 0.51099890221e-03  # Relative Energy
    # 2nd order statistical moments
    elecBeam.arStatMom2[0] = (137.113e-06) ** 2  # <(x-x0)^2>
    elecBeam.arStatMom2[1] = -0.0388489e-09  # <(x-x0)*(x'-x'0)>
    elecBeam.arStatMom2[2] = (6.57004e-06) ** 2  # <(x'-x'0)^2>

    elecBeam.arStatMom2[3] = (5.39499e-06) ** 2  # (15.4091e-06)**2 #<(y-y0)^2>
    elecBeam.arStatMom2[4] = -0.00211765e-09  # <(y-y0)*(y'-y'0)>
    elecBeam.arStatMom2[5] = (1.53393e-06) ** 2  # <(y'-y'0)^2>
    elecBeam.arStatMom2[10] = (0.89e-03) ** 2  # <(E-E0)^2>/E0^2

    # elecBeam = drift_ebeam(elecBeam,-0.5)
    sys.stdout.write('done\n')

    meth = 1  # SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    relPrec = 0.008  # 0.001 #relative precision
    zStartInteg = 0  # longitudinal position to start integration (effective if < zEndInteg)
    zEndInteg = 0  # longitudinal position to finish integration (effective if > zStartInteg)

    npTraj = 50000  # Number of points for trajectory calculation

    useTermin = 1  # Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
    sampFactNxNyForProp = 0.3  # 1.0 #0.2 #0.1 #sampling factor for adjusting nx, ny (effective if > 0) ###TODO fix that
    arPrecParSpec = [meth, relPrec, zStartInteg, zEndInteg, npTraj, useTermin, sampFactNxNyForProp]

    # ***********UR Wavefront Parameters (mesh)
    sys.stdout.write('   Setup Stokes mesh ... ');
    sys.stdout.flush()
    wfr = SRWLWfr()  # for spectral flux vs photon energy
    wfr.allocate(1, 81, 81)  # numbers of points vs photon energy, horizontal and vertical positions

    wfr.mesh.zStart = APE
    wfr.mesh.eStart = 2101.  # initial photon energy [eV]
    wfr.mesh.eFin = 2101.  # final photon energy [eV]wfr.mesh.xStart = -0.0006 #initial horizontal position [m]
    wfr.mesh.xStart = -0.0006  # initial horizontal position [m]
    wfr.mesh.xFin = 0.0006  # final horizontal position [m]
    wfr.mesh.yStart = -0.0006  # initial vertical position [m]
    wfr.mesh.yFin = 0.0006  # final vertical position [m]
    wfr.partBeam = elecBeam
    meshInitPartCoh = deepcopy(wfr.mesh)
    sys.stdout.write('done\n')

    # ***********Optical Elements
    delta = 7.80723421E-05  # Be @ 2.1KeV
    attenLen = 84.3e-06  # [m] #2.1KeV
    diamCRL = 1.e-03  # CRL diameter
    rMinCRL = 200e-06  # CRL radius at the tip of parabola [m]
    nCRL = 1  # number of lenses
    wallThickCRL = 35e-06  # CRL wall thickness [m]

    # Generating a perfect 2D parabolic CRL:
    #   
    CRL = srwl_opt_setup_CRL(3, delta, attenLen, 1, diamCRL, diamCRL, rMinCRL, nCRL, wallThickCRL, 0, 0)
    #
    # Beamline OEs
    # APE='aperture',
    # MOAT='first mirror of Monocromator error shape',
    # VFML='Vertical Focusing Mirror (Spherical) Lens',
    # VFMT='Vertical Focusing Mirror (Spherical) error shape'

    D_APE_MOA = SRWLOptD(DCM)

    # ================introducing Si(111) heat load==============
    heightProf = AuxReadInDataColumns(strProfileData, 2, '\t')
    # MOAT        = SRWLOptT(100, 2001, 2.0e-02, 3.0e-02*sin(1.223866));
    MOAT = SRWLOptT(100, 1000, 2.0e-02, 3.797e-3 * sin(1.223866))
    AuxTransmAddSurfHeightProfile(MOAT, heightProf, 'y', 1.223866)  # incident angle is 70.122373 deg => 1.223866 rad
    opPathDifMOAT = MOAT.get_data(3, 3)
    srwl_uti_save_intens_ascii(opPathDifMOAT, MOAT.mesh,
                               os.path.join(os.getcwd(), strExDataFolderName, 'res_er_mono.dat'), 0,
                               ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                               _arUnits=['', 'm', 'm', 'm'])
    # ================finished Si (111) heat load================

    D_MOA_HFM = SRWLOptD(HFM)
    D_APE_HFM = SRWLOptD(DCM + HFM)  # used for "no bump" option
    if BMmode == 'Norm':
        HFML = SRWLOptL(_Fx=1. / (1. / (APE + DCM + HFM) + 1. / (VFM + VM + SSA + Sample)))  # to focus at ES1
    if BMmode == 'LowDiv':
        #    HFML        = SRWLOptL(_Fx=1./(1./(APE+DCM+HFM)+1./(100))) # to focus at ES2 with a low divergence
        HFML = SRWLOptL(_Fx=34.88)  # to focus at ES2 with a low divergence

    # ================introducing HFM slope error==============
    heightProfHFM = AuxReadInDataColumns(strProfileDataHFM, 2, '\t')
    HFMT = SRWLOptT(464, 200, 0.5 * sin(3.1415927e-03),
                    6.0e-03)  # sinusoidal 0.1, 1.8e-08 both 'h' 'v', angle 3.1415927e-03 rad to correct for horizontal.
    AuxTransmAddSurfHeightProfile(HFMT, heightProfHFM, 'x', 3.1415927e-03)
    opPathDifHFMT = HFMT.get_data(3, 3)
    srwl_uti_save_intens_ascii(opPathDifHFMT, HFMT.mesh,
                               os.path.join(os.getcwd(), strExDataFolderName, 'res_er_HFM.dat'), 0,
                               ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                               _arUnits=['', 'm', 'm', 'm'])
    # ================finished HFM slope error=================
    # HFMT = SRWLOptT()
    D_HFM_VFM = SRWLOptD(VFM)
    if BMmode == 'Norm':
        VFML = SRWLOptL(_Fy=1. / (1. / (APE + DCM + HFM + VFM - 0.6) + 1. / (
            VM + SSA + Sample + 0.9)))  # focus at ES1; if using Bump, VFM must be 3.9+0.9 m (to compensate bump which moves focus 0.8m upstream)
    if BMmode == 'LowDiv':
        #    VFML        = SRWLOptL(_Fy=1./(1./(APE+DCM+HFM+VFM-0.6)+1./(100))) #focus at ES2 with a low divergence
        VFML = SRWLOptL(_Fy=47.5)  # focus at ES2 with a low divergence

    # ================introducing VFM slope error==============
    heightProfVFM = AuxReadInDataColumns(strProfileDataVFM, 2, '\t')
    VFMT = SRWLOptT(200, 464, 6.0e-03, 0.5 * sin(
        6.1086524e-03))  # sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
    AuxTransmAddSurfHeightProfile(VFMT, heightProfVFM, 'y', 6.1086524e-03)
    opPathDifVFMT = VFMT.get_data(3, 3)
    srwl_uti_save_intens_ascii(opPathDifVFMT, VFMT.mesh,
                               os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VFM.dat'), 0,
                               ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                               _arUnits=['', 'm', 'm', 'm'])
    # ================finished VFM slope error=================
    #
    # ================introducing VM slope error==============
    heightProfVM = AuxReadInDataColumns(strProfileDataVM, 2, '\t')
    VMT = SRWLOptT(200, 500, 6.0e-03, 0.5 * sin(
        6.1086524e-03))  # sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
    AuxTransmAddSurfHeightProfile(VMT, heightProfVM, 'y', 6.1086524e-03)
    opPathDifVMT = VMT.get_data(3, 3)
    srwl_uti_save_intens_ascii(opPathDifVMT, VMT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VM.dat'),
                               0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'],
                               _arUnits=['', 'm', 'm', 'm'])
    # ================finished VFM slope error=================
    # VFMT = SRWLOptT()
    # D_VFM_SSA   = SRWLOptD(VM + SSA + Sample) #TODO Where is the SSA?
    D_VFM_VM = SRWLOptD(VM)
    D_VM_SSA = SRWLOptD(SSA)
    D_VFM_SSA = SRWLOptD(VM + SSA)
    if beamline == 'ES1' and BMmode == 'Norm':
        SSA = SRWLOptA('r', 'a', 0.4e-03, 0.4e-03)  # 0.4, 0.4 for NOT low divergence mode;

    if beamline == 'ES2' and BMmode == 'Norm':
        SSA = SRWLOptA('r', 'a', 0.9e-03, 0.9e-03)  # 0.4x0.15 for KB

    if beamline == 'ES2' and BMmode == 'LowDiv':
        SSA = SRWLOptA('r', 'a', 0.9e-03, 0.9e-03)  # 0.4, 0.4 for low divergence mode;

    D_SSA_ES1 = SRWLOptD(Sample)  # TODO Where is the SSA?
    D_SSA_CRL = SRWLOptD(CRL_d)
    D_CRL_ES2 = SRWLOptD(1.5948)
    # D_ES1_HKB   = SRWLOptD(6.7)
    # D_SSA_HKB   = SRWLOptD(Sample + 6.7)

    ApCRL = SRWLOptA('c', 'a', 1.0e-3)
    # CRL = SRWLOptL(_Fx=1./(1./(6.0)+1./(2.1)), _Fy=1./(1./(6.0)+1./(2.1)))

    # angHKB = 3.14e-03 #[rad]
    # HKB = SRWLOptMirEl(_p=6.7, _q=1.40, _ang_graz=angHKB, _r_sag=1.e+40, _size_tang=0.3, _nvx=cos(angHKB), _nvy=0, _nvz=-sin(angHKB), _tvx=-sin(angHKB), _tvy=0, _x=0, _y=0, _treat_in_out=1) #HKB Ellipsoidal Mirror
    # ================introducing HKB slope error==============
    # heightProfHKB  = AuxReadInDataColumns(strProfileDataHKB, 2, '\t')
    # HKBT        = SRWLOptT(127, 200, 0.5166*sin(3.1415927e-03), 6.0e-03) #sinusoidal 0.1, 1.8e-08 both 'h' 'v', angle 3.1415927e-03 rad to correct for horizontal.
    # AuxTransmAddSurfHeightProfile(HKBT, heightProfHKB, 'x', 3.1415927e-03)
    # opPathDifHKBT = HKBT.get_data(3, 3)
    # srwl_uti_save_intens_ascii(opPathDifHKBT, HKBT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_HKB.dat'), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
    # ================finished HKB slope error=================
    # angVKB = 3.14e-03 #grazing angle at VKB center [rad]
    # VKB = SRWLOptMirEl(_p=7.1, _q=1., _ang_graz=angVKB, _r_sag=1.e+40, _size_tang=0.35, _nvx=0, _nvy=cos(angVKB), _nvz=-sin(angVKB), _tvx=0, _tvy=-sin(angVKB), _x=0, _y=0, _treat_in_out=1) #VKB Ellipsoidal Mirror
    # ================introducing VKB slope error==============
    # heightProfVKB = AuxReadInDataColumns(strProfileDataVKB, 2, '\t')
    # VKBT        = SRWLOptT(200, 500, 6.0e-03, 0.5*sin(3.1415927e-03))#sinusoidal equal to HFM. the origina spec is 0.1, 6.75e-09 both 'h' 'v', angle 6.1086524e-03 rad to correct for vertical.
    # AuxTransmAddSurfHeightProfile(VKBT, heightProfVKB, 'y', 3.1415927e-03)
    # opPathDifVKBT = VKBT.get_data(3, 3)
    # srwl_uti_save_intens_ascii(opPathDifVKBT, VKBT.mesh, os.path.join(os.getcwd(), strExDataFolderName, 'res_er_VKB.dat'), 0, ['', 'Horizontal Position', 'Vertical Position', 'Opt. Path Diff.'], _arUnits=['', 'm', 'm', 'm'])
    # ================finished VKB slope error=================

    # D_HKB_VKB = SRWLOptD(0.4) #Distance between centers of Vertically and Horizontally focusing K-B mirrors
    # D_VKB_ES2 = SRWLOptD(1.0)
    # D_ES1_ES2 = SRWLOptD(8.1)

    #               [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
    if beamline == 'ES1':
        # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 2.0, 6.0, 2.0, 6.0,  0,  0,  0 ]   #settings for single electron emission
        ppD_APE_MOA = [0, 0, 1.0, 2, 0, 8.0, 8.0, 2.0, 6.0, 0, 0, 0]  # settings for single electron emission

        # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 8.0, 6.0, 2.0, 10.0,  0,  0,  0 ]  #settings for multi-electron emission
    if beamline == 'ES2':
        ppD_APE_MOA = [0, 0, 1.0, 2, 0, 8.0, 4.0, 4.0, 8.0, 0, 0, 0]  # settings for single electron emission
        # ppD_APE_MOA = [ 0,  0, 1.0,  2,  0, 8.0, 2.0, 2.0, 8.0,  0,  0,  0 ]   #settings for multi-electron emission

    ppMOAT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_MOA_HFM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_APE_HFM = ppD_APE_MOA
    ppHFML = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppHFMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_HFM_VFM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppVFML = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppVFMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppVMT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_VFM_SSA = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_VFM_VM = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_VM_SSA = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppSSA = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_SSA_CRL = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppApCRL = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppCRL = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_CRL_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_SSA_ES1 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_ES1_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    # downstream of ES1
    ppD_ES1_HKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_SSA_HKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppHKB = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppHKBT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_HKB_VKB = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppVKB = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppVKBT = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]
    ppD_VKB_ES2 = [0, 0, 1.0, 2, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

    # [ 0]: Auto-Resize (1) or not (0) Before propagation
    # [ 1]: Auto-Resize (1) or not (0) After propagation
    # [ 2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
    # [ 3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation (2) semi-analytical treatment + autoresizing
    # [ 4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
    # [ 5]: Horizontal Range modification factor at Resizing (1. means no modification)
    # [ 6]: Horizontal Resolution modification factor at Resizing
    # [ 7]: Vertical Range modification factor at Resizing
    # [ 8]: Vertical Resolution modification factor at Resizing
    # [ 9]: Type of wavefront Shift before Resizing (not yet implemented)
    # [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
    # [11]: New Vertical wavefront Center position after Shift (not yet implemented)

    # with transission for error shape
    # BML         = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA ],
    #                        [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA ] )
    # without error shape
    # if beamline=='ES1' and not bump:
    #   BML     = SRWLOptC( [   D_APE_HFM,   HFML,   D_HFM_VFM,   VFML,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
    #                      [ ppD_APE_HFM, ppHFML, ppD_HFM_VFM, ppVFML, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])
    # if beamline=='ES1' and bump:
    #   BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   D_HFM_VFM,   VFML,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
    #                      [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppD_HFM_VFM, ppVFML, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])
    if beamline == 'ES1' and bump and BMmode == 'Norm':
        # BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA,   SSA,   D_SSA_ES1 ],
        #                    [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA, ppSSA, ppD_SSA_ES1 ])
        BML = SRWLOptC(
            [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_ES1],
            [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT,
             ppD_VM_SSA, ppSSA, ppD_SSA_ES1])

        # if beamline=='ES2' and not bump:
        #   BML     = SRWLOptC( [   D_APE_HFM,   HFML,   HFMT,   D_HFM_VFM,   VFML,   VFMT,   D_VFM_SSA,   SSA,   D_SSA_HKB,   HKB,   D_HKB_VKB,   VKB,   D_VKB_ES2 ],
        #                      [ ppD_APE_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_SSA, ppSSA, ppD_SSA_HKB, ppHKB, ppD_HKB_VKB, ppVKB, ppD_VKB_ES2 ] )

        # if beamline =='ES2' and bump and BMmode == 'LowDiv': #focus ES2 without Kb with low divergence
        #   BML     = SRWLOptC( [   D_APE_MOA,   MOAT,   D_MOA_HFM,   HFML,   HFMT,  D_HFM_VFM,   VFML,    VFMT,   D_VFM_VM,   VMT,   D_VM_SSA,   SSA,   D_SSA_CRL,   ApCRL,   CRL,   D_CRL_ES2 ],
        #                      [ ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT, ppD_VM_SSA, ppSSA, ppD_SSA_CRL, ppApCRL, ppCRL, ppD_CRL_ES2 ] )
    if beamline == 'ES2' and bump and BMmode == 'LowDiv':  # focus ES2 without Kb with low divergence
        BML = SRWLOptC(
            [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_ES1,
             D_ES1_ES2],
            [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT,
             ppD_VM_SSA, ppSSA, ppD_SSA_ES1, ppD_ES1_ES2])

    if beamline == 'ES2' and bump and BMmode == 'Norm':
        BML = SRWLOptC(
            [D_APE_MOA, MOAT, D_MOA_HFM, HFML, HFMT, D_HFM_VFM, VFML, VFMT, D_VFM_VM, VMT, D_VM_SSA, SSA, D_SSA_CRL,
             ApCRL, CRL, D_CRL_ES2],
            [ppD_APE_MOA, ppMOAT, ppD_MOA_HFM, ppHFML, ppHFMT, ppD_HFM_VFM, ppVFML, ppVFMT, ppD_VFM_VM, ppVMT,
             ppD_VM_SSA, ppSSA, ppD_SSA_CRL, ppApCRL, ppCRL, ppD_CRL_ES2])


    # CALCULATION

    # **********************Calculation (SRWLIB function calls)

    # ***Electron Trajectory
    partTraj = SRWLPrtTrj()
    partTraj.partInitCond = elecBeam.partStatMom1
    partTraj.allocate(npTraj)
    partTraj.ctStart = elecBeam.partStatMom1.z  # Start Time for the calculation
    # partTraj.ctEnd = (numPer + 2)*per + magFldCnt.arMagFld[0].rz + magFldCnt.arMagFld[2].rz #End Time
    partTraj.ctEnd = 4

    # print('   Performing trajectory calculation ... ', end='')
    # partTraj = srwl.CalcPartTraj(partTraj, magFldCnt, [1])
    # AuxSaveTrajData(partTraj, os.path.join(os.getcwd(), strExDataFolderName, strTrajOutFileName))
    # print('done')

    # sys.exit(0)

    sys.stdout.write('   Performing Single Electron calculation ... ');
    sys.stdout.flush()
    srwl.CalcElecFieldSR(wfr, 0, magFldCnt, arPrecParSpec)
    sys.stdout.write('done\n')

    sys.stdout.write('   Saving Single Electron UR Intensity ... ');
    sys.stdout.flush()
    arI = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    # srwl.CalcIntFromElecField(arI, wfr, 6, 1, 3, wfr.mesh.eStart, 0, 0)

    srwl_uti_save_intens_ascii(arI, wfr.mesh, strIntSourSE_OutFileName)
    sys.stdout.write('done\n')

    # sys.exit(0)

    sys.stdout.write('   Performing Single Electron Radiation Propagation ... ');
    sys.stdout.flush()
    srwl.PropagElecField(wfr, BML)
    sys.stdout.write('done\n')

    sys.stdout.write('   Saving Single Electron Propagated Intensity ... ');
    sys.stdout.flush()
    arI = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
    srwl_uti_save_intens_ascii(arI, wfr.mesh, strIntPropSE_OutFileName)
    sys.stdout.write('done\n')

    print('Switching from Coordinate to Angular Representation ... ', end='')
    srwl.SetRepresElecField(wfr, 'a');
    print('done')

    print('Extracting Intensity from the Propagated Electric Field in Angular Representation  ... ', end='')
    arIa = array('f', [0] * wfr.mesh.nx * wfr.mesh.ny)  # "flat" 2D array to take intensity data
    srwl.CalcIntFromElecField(arIa, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)

    # Converting intensity in angular representation to [photons/s/.1%bw/mrad^2] vs rad
    meshAng = deepcopy(wfr.mesh)
    wavelength = 1.239842e-06 / meshAng.eStart  # wavelength in [m] - converiosn / multiplier for argument
    invSqWavelength = 1 / (wavelength * wavelength)  # converiosn / multiplier for intensity values
    meshAng.xStart *= wavelength
    meshAng.xFin *= wavelength
    meshAng.yStart *= wavelength
    meshAng.yFin *= wavelength
    for i in range(len(arIa)):
        arIa[i] *= invSqWavelength

    srwl_uti_save_intens_ascii(arIa, meshAng, os.path.join(os.getcwd(), strExDataFolderName, "prop_se_ang.dat"))
    print('done')

    '''
    sys.exit(0)

    sys.stdout.write('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... '); sys.stdout.flush()
    nMacroElec = 800000 #Total number of Macro-Electrons (Wavefronts)
    nMacroElecAvgOneProc = 5 #Number of Macro-Electrons (Wavefronts) to average on each node (for MPI calculations)
    nMacroElecSavePer = 5 #Saving periodicity (in terms of Macro-Electrons) for the Resulting Intensity
    #arPrecParSpec[6] = sampFactNxNyForProp
    radStokesProp = srwl_wfr_emit_prop_multi_e(elecBeam, magFldCnt, meshInitPartCoh, 1, 0.01, nMacroElec, nMacroElecAvgOneProc, nMacroElecSavePer,
                                               strIntPropME_OutFileName, sampFactNxNyForProp, BML)
    sys.stdout.write('done\n')
    '''
    return os.path.join(os.getcwd(), strExDataFolderName, "prop_se_ang.dat")


# -------------------------------------------------------------------------------


# GUI:
# TODO: check http://stackoverflow.com/a/24728212/4143531 for Python 3.* solution:
from ttk import Frame, Button, Label
from Tkinter import Tk, N, S, W, E, BOTH, RIGHT, DoubleVar, BooleanVar, StringVar, Entry, Checkbutton, OptionMenu
import tkMessageBox as box
import tkFileDialog
from uti_plot import *

from PIL import ImageTk, Image


class Example(Frame):
    def __init__(self, parent):
        Frame.__init__(self, parent)
        self.parent = parent
        self.w = 700
        self.h = 350

        script_dir = os.path.dirname(os.path.realpath(__file__))
        os.chdir(script_dir)
        self.images_dir = os.path.join(script_dir, 'images')
        self.elements_list = [
            'APE',
            'DCM',
            'HFM',
            'VFM',
            'VM',
            'SSA',
            'CRL',
            'Sample',
        ]
        self.default_dist = {
            'APE': 29.5,
            'DCM': 2.44,
            'HFM': 3.023696,
            'VFM': 3.42,
            'VM': 0.7,
            'SSA': 8.0,
            'CRL': 10.4052,
            'Sample': 3.9,
        }

        self.counter = 0
        self.images = []
        self.entries = []
        self.mirror = None
        self.plot_graphs = None
        self.beamline = None

        self.initUI()
        self.centerWindow()

        self.parent.minsize(width=self.w, height=self.h)
        self.parent.minsize(height=self.h)

        self.parent.maxsize(width=self.w, height=self.h)
        self.parent.maxsize(height=self.h)

    def centerWindow(self):
        sw = self.parent.winfo_screenwidth()
        sh = self.parent.winfo_screenheight()

        x = sw * 3 / 4 - self.w / 2
        y = sh * 1 / 2 - self.h / 2
        self.parent.geometry('%dx%d+%d+%d' % (self.w, self.h, x, y))

    def initUI(self):
        self.parent.title("SRW GUI")
        # self.style = Style()
        # self.style.theme_use("default")

        self.pack(fill=BOTH, expand=1)

        # Just empty row to make layout better:
        label = Label(self, text='').grid(row=0)

        # Add descriptions for the rows:
        label = Label(self, text='Beamline\nelements:').grid(row=1, column=0)
        label = Label(self, text='Distance (m):').grid(row=4, column=0, sticky=E)

        # Add buttons for the elements:
        for i in range(len(self.elements_list)):
            self.add_button(self.elements_list[i], i + 1)

        # Just empty row to make layout better:
        label = Label(self, text='').grid(row=5)

        # Add static image of the undulator:
        self.static_element('Undulator', 0)

        # Print button to print gathered values:
        print_button = Button(self, text="Print", command=self.on_print).grid(row=6, column=0, sticky=N + E)
        # Delete last added element:
        del_button = Button(self, text="Delete last", command=self.on_delete).grid(row=7, column=0, sticky=E)
        # Delete all elements:
        del_button = Button(self, text="Clear all", command=self.on_clear).grid(row=8, column=0, sticky=S + E)

        # Check button for plotting MPL plots after execution:
        self.plot_graphs = BooleanVar()
        check_button = Checkbutton(self, text="Plot graphs?", variable=self.plot_graphs)
        check_button.select()
        check_button.grid(row=7, column=3, columnspan=2)

        # Add mirror file:
        self.mirror = Button(self, text="Browse mirror file", command=self.on_browse).grid(row=7, column=5,
                                                                                           columnspan=2, sticky=E + W)

        # Add selection of beamline:
        options = ['ES1', 'ES2']
        self.beamline = StringVar()
        self.beamline.set(options[-1])
        beamline = apply(OptionMenu, (self, self.beamline) + tuple(options))
        beamline.grid(row=1, column=i + 2)
        Label(self, text='Beamline:').grid(row=1, column=i + 2, sticky=N)

        # Close and OK buttons:
        close_button = Button(self, text="Close", command=self.parent.destroy).grid(row=9, column=0, sticky=S + E)
        ok_button = Button(self, text="OK", command=self.on_ok).grid(row=9, column=i + 2, sticky=S + W)

        # Configure sizes of rows and columns:
        for j in range(i):
            self.columnconfigure(j + 1, minsize=66)

        self.rowconfigure(3, minsize=80)  # row with images of elements appearing after clicking buttons
        self.columnconfigure(0, minsize=90)  # First column with all buttons
        self.columnconfigure(i + 2, minsize=80)  # Last column with 'Cancel' button
        self.rowconfigure(9, minsize=50)  # last row with 'OK' button

    # Events processing:
    def on_click(self, name):
        self.counter += 1
        self.add_element(name, self.counter)

    def on_delete(self):
        if len(self.images) > 0:
            self.images[-1].destroy()
            del self.images[-1]
            self.counter -= 1

        if len(self.entries) > 0:
            self.entries[-1][0].destroy()
            del self.entries[-1]

    def on_clear(self):
        for i in range(len(self.images)):
            self.images[i].destroy()
            self.entries[i][0].destroy()

        self.images = []
        self.entries = []
        self.counter = 0

    def add_button(self, name, col):
        img = self.resize_image(name)

        button = Button(self, image=img, command=lambda: self.on_click(name))
        button.image = img
        button.grid(row=1, column=col)

        label = Label(self, text='Add ' + name)
        label.grid(row=2, column=col)

        return button

    def add_element(self, name, col):
        img = self.resize_image(name)
        image = Label(self, image=img)
        image.image = img
        image.grid(row=3, column=col)

        var = DoubleVar()
        entry = Entry(self, textvariable=var, width=8, justify=RIGHT)
        var.set(self.default_dist[name])
        entry.selection_adjust(len(str(self.default_dist[name])))
        entry.focus_set()
        entry.grid(row=4, column=col)

        self.images.append(image)
        self.entries.append((entry, name))

    def static_element(self, name, col):
        img = self.resize_image(name, height=80)
        image = Label(self, image=img)
        image.image = img
        image.grid(row=3, column=col, sticky=N + S)

    def resize_image(self, name, width=None, height=50):
        img = Image.open(self.images_dir + '/' + name + '.png')
        if not width:
            scale_coef = (height / float(img.size[1]))
            width = int(float(img.size[0]) * float(scale_coef))
        img = img.resize((width, height), Image.ANTIALIAS)
        img = ImageTk.PhotoImage(img)
        return img

    def on_print(self):
        print('Distances  :', ', '.join([x[0].get() for x in self.entries]))
        print('Elements   :', ', '.join([x[1] for x in self.entries]))
        print('Plot graphs:', self.plot_graphs.get())
        print('Mirror file:', self.mirror)
        print('Beamline   :', self.beamline.get())

    def on_scale(self, val):
        self.var1.set(int(float(self.scale.get())))

    def on_browse(self):
        ftypes = [('DAT files', '*.dat'), ('All files', '*')]
        dlg = tkFileDialog.Open(self, filetypes=ftypes)
        fl = dlg.show()
        if fl != '':
            self.mirror = os.path.abspath(fl)
            label = Label(self, text=os.path.basename(self.mirror))
            label.grid(row=7, column=7, columnspan=2, sticky=W)

    def on_ok(self):
        GUI_data = {}
        for i in range(len(self.entries)):
            GUI_data[self.entries[i][1]] = float(self.entries[i][0].get())

        '''
        GUI_data = {
            'npTraj': self.var1.get(),
            'phBy': self.var2.get(),
            'numPer': self.var3.get(),
            'BMmode': self.var4.get(),
            'plot_graphs': self.var5.get(),
            'mirror_file': self.var6.get(),
            'spinbox': self.var7.get(),
            'radiobutton': self.var8.get(),
        }
        '''
        GUI_data['mirror_file'] = self.mirror
        GUI_data['beamline'] = self.beamline.get()
        GUI_data['defaults'] = self.default_dist

        a = '\n' + '=' * 80 + '\n'
        a += '%40s' % ('Parameters:') + '\n'
        a += '-' * 80 + '\n'
        keys = sorted(GUI_data.keys())
        for key in keys:
            a += '%-15s : %30s\n' % (key, GUI_data[key])
        a += '=' * 80 + '\n'
        print(a)

        # Run the function for beamline calculation:
        outfile = parm_func(GUI_data)

        # Show a window with information about completion:
        msg = "Calculation completed, data is saved in:\n" + outfile
        box.showinfo("Information", msg)
        print(msg)

        # Close the window:
        self.parent.destroy()

        # Plot output information:
        if self.plot_graphs.get():
            try:
                uti_data_file_plot(outfile, 0, 0, 0, 0, 0)
                uti_plot_show()
            except:
                print('Outfile is not found or plotting program met an error.')


def main():
    root = Tk()
    ex = Example(root)
    # set the window background to hex code '#a1dbcd'
    # root.configure(background="#a1dbcd")
    # root.configure(background="black")
    # root.geometry("500x300+300+300")
    root.mainloop()


# End of functions part
# -------------------------------------------------------------------------------

if __name__ == '__main__':
    main()
