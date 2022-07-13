from scipy.ndimage import binary_dilation, binary_erosion
from radiomics import featureextractor
from mahotas.labeled import bwperim
from skimage.measure import label
from skimage import exposure
from pywt import dwtn, idwtn
import scipy.ndimage as ndimage
import numpy.ma as ma
import pandas as pd
import numpy as np
import tempfile
import mahotas
import nrrd

def getDifferentialUniformity(rawCT):

    rd = np.zeros(4)
    cdiff = np.zeros(4)
    sd = np.zeros(4)
    maxcd = 0
    maxrd = 0
    maxsd = 0
    mincd = 10001
    minrd = 10001
    minsd = 10001
    k = 1


    dimRows, dimColumns, dimSlices = np.shape(rawCT)

    for s in np.arange(dimSlices):
        for r in np.arange(dimRows):
            for c in np.arange(dimColumns):
                for m in np.arange(4):
                    if r+m < dimRows-1:
                        rd[m] = np.abs(rawCT[r,c,s] - rawCT[r+m+1,c,s])
                    if c+m < dimColumns-1:
                        cdiff[m] = np.abs(rawCT[r,c,s] - rawCT[r,c+m+1,s])
                    if s+m < dimSlices-1:
                        sd[m] = np.abs(rawCT[r,c,s] - rawCT[r,c,s+m+1])
                maxcd = np.maximum(maxcd,np.max(rd))
                maxrd = np.maximum(maxrd,np.max(rd))
                maxsd = np.maximum(maxsd,np.max(sd))
                mincd = np.minimum(mincd,np.min(cdiff, where=cdiff>0, initial=10000))
                minrd = np.minimum(minrd,np.min(rd, where=rd>0, initial=10000))
                minsd = np.minimum(minsd,np.min(sd, where=sd>0, initial=10000))
                rd = np.zeros(4)
                cdiff = np.zeros(4)
                sd = np.zeros(4)
                print(r)

    maxUni = np.max(maxcd,maxrd,maxsd)
    minUni = np.max(mincd, minrd, minsd)

    DifferentialUniformity = (maxUni-minUni)/(maxUni+minUni)

    return DifferentialUniformity


def getEccentricity(rawCT, metaCT, mask):
    mask = mask.astype(bool).astype('uint8')
    mask_bool = np.array(mask, dtype=bool)
    mask_bool = mask_bool.astype('uint8')

    ROI = rawCT.copy()
    ROI[mask_bool == 0] = 0
    ROI[mask_bool != 0] = rawCT[mask_bool != 0]

    pixelW = metaCT["space directions"][0,0]
    sliceS = metaCT["space directions"][2,2]

    # ISOTROPIC RESAMPLING
    sFactor = sliceS/pixelW
    resampledMask = ndimage.interpolation.zoom(ROI, zoom=[1,1,sFactor])

    # INITIALIZATION OF ECCENTRICITY COMPUTATION
    perimeter = mahotas.bwperim(resampledMask)
    nPoints = np.shape(np.where(perimeter==1))[1]
    X = np.zeros([nPoints,1])
    Y = np.zeros([nPoints,1])
    Z = np.zeros([nPoints,1])

    counter = 0
    for i in np.arange(np.shape(perimeter)[2]):
        row, col = np.where(perimeter[:,:,i] == 1)
        p = len(row)
        if p > 0:
            X[counter:(counter+p),0] = col
            Y[counter:(counter+p),0] = row
            Z[counter:(counter+p),0] = i
            counter = counter + p

    # LI'S ELLIPSOID FITTING ALGORITHM
    dx = X
    dy = Y
    dz = Z
    n = len(dx)
    D = np.array([dx*dx,dy*dy,dz*dz,2*dy*dz,2*dx*dz,2*dx*dy,2*dx,2*dy,2*dz,np.ones([n,1])])
    D = D[:,:,0]
    S = np.matmul(D,D.T)

    # Create constraint matrix C:
    C = np.zeros([6,6])
    C[0,0] = -1
    C[1,1] = -1
    C[2,2] = -1
    C[3,3] = -4
    C[4,4] = -4
    C[5,5] = -4
    C[0,1] = 1
    C[1,0] = 1
    C[0,2] = 1
    C[2,0] = 1
    C[1,2] = 1
    C[2,1] = 1

    # Solve generalized eigensystem
    S11 = S[0:6, 0:6]
    S12 = S[0:6, 6:10]
    S22 = S[6:10, 6:10]

    A = S11-np.matmul(np.matmul(S12,np.linalg.pinv(S22)),S12.T)
    CA = np.matmul(np.linalg.inv(C),A)
    geval, gevec = np.linalg.eig(CA)

    # Find the largest eigenvalue
    In = 0
    maxVal = geval[0]
    for i in np.array([1,2,3,4,5]):
        if (geval[i]>maxVal):
            maxVal = geval[i]
            In = i

    v1 = gevec[:,In]
    v2 = np.matmul(np.matmul(-np.linalg.pinv(S22),S12.T),v1)
    v = np.concatenate((v1,v2))

    # Algebraic form of the ellipsoid
    A = np.array([[v[0], v[5], v[4], v[6]], [v[5], v[1],v[3],v[7]], [v[4],v[3],v[2], v[8]],[v[6], v[7], v[8], -1]])

    # Center of the ellipsoid
    center = np.linalg.solve(-A[0:3,0:3],np.array([v[6], v[7], v[8]]))
    center = np.sort(center)[::-1]

    # Corresponding translation matrix
    T = np.eye(4)
    T[3,0:3] = center.T

    # Translating to center
    R = np.matmul(np.matmul(T,A),T.T)

    # Solving eigenproblem 
    evals, evecs = np.linalg.eig((R[0:3,0:3]/-R[3,3]))
    radii = np.sqrt(1/evals)
    radii = np.sort(radii)[::-1]

    # ECCENTRICY COMPUTATION
    eccentricity = np.sqrt(1 - (radii[1]*radii[2])/(radii[0]**2))

    return eccentricity


def getPercentInactive(rawCT, mask, thresh=0.6):
    
    mask = mask.astype(bool).astype('uint8')
    mask_bool = np.array(mask, dtype=bool)
    mask_bool = mask_bool.astype('uint8')

    ROI = rawCT.copy()
    ROI[mask_bool == 0] = 0
    ROI[mask_bool != 0] = rawCT[mask_bool != 0]

    maskInactive = np.zeros(np.shape(ROI))
    maskInactive[ROI > (thresh*np.max(ROI)**2)] = 1

    # if (maskInactive).any() != 0:
    #     conn = np.zeros([5,5,5])
    #     conn[:,:,0] = np.array([[0,0,0,0,0], [0,0,0,0,0], [0,0,1,0,0], [0,0,0,0,0], [0,0,0,0,0]])
    #     conn[:,:,4] = conn[:,:,0]
    #     conn[:,:,1] = np.array([[0,0,0,0,0], [0,0,1,0,0], [0,1,1,1,0], [0,0,1,0,0], [0,0,0,0,0]])
    #     conn[:,:,2] = np.array([[0,0,1,0,0], [0,1,1,1,0], [1,1,1,1,1], [0,1,1,1,0], [0,0,1,0,0]])
    #     conn[:,:,3] = conn[:,:,1]

    #     maskInactive = cv2.morphologyEx(maskInactive, cv2.MORPH_CLOSE, conn)
    #     maskInactive = cv2.morphologyEx(maskInactive, cv2.MORPH_OPEN, conn)

    perim = bwperim(mask, 8)
    newMask = maskInactive + perim
    newMask[mask == 0] = 10
    newMask[newMask == 1] = 10
    newMask[newMask == 2] = 10
    newMask[newMask == 0] = 1
    newMask[newMask == 10] = 0

    all_labels = label(newMask)
    numLabels = np.shape(np.unique(all_labels))[0]

    b = np.zeros([1, numLabels-1])
    for l in np.arange(numLabels-1):
        if np.shape(np.where(all_labels==l))[1]<15:
            b[0,l] = 0
        else:
            b[0,l] = 1

    [row, col] = np.where(b>0)

    sumInactive = 0
    for i in np.arange(np.size(col)):
        sumInactive = sumInactive + np.shape(np.where(all_labels==(i+1)))[1]

    sumVolume = np.sum(mask)
    percentInactive = sumInactive/sumVolume * 100
    
    return percentInactive


def getSUVMetricP(rawCT, mask):

    nbins = 1000
    mask_bool = np.array(mask, dtype=bool)
    ROI = rawCT[mask_bool]

    HUmax = np.round(np.max(ROI))
    HUmean = np.round(np.mean(ROI))
    HUmin = np.round(np.min(ROI))
    HUsd = np.round(np.std(ROI))
    range = np.round(HUmax-HUmin)
    HUmedian = np.round(np.median(ROI))
    cov = 100*HUsd/HUmean
    cov = np.round(np.abs(cov))

    # SUVpeak
    radius = 2
    rawCTMasked = ma.masked_array(rawCT, mask=np.invert(mask_bool))
    x0,y0,z0 = np.unravel_index(np.argmax(rawCTMasked), np.shape(rawCTMasked))
    SUVPeakMask = np.zeros(np.shape(rawCTMasked))

    for x in np.arange(x0-radius, x0+radius+1):
        for y in np.arange(y0-radius, y0+radius+1):
            for z in np.arange(z0-radius, z0+radius+1):
                deb = radius - abs(x0-x) - abs(y0-y) - abs(z0-z) 
                if (deb)>=0: 
                    SUVPeakMask[x,y,z] = 1
    
    SUVPeakMask = np.array(SUVPeakMask, dtype=bool)
    rawCTMaskedTest = ma.masked_array(rawCT, mask=np.invert(SUVPeakMask))
    HUpeak = np.round(np.mean(rawCTMaskedTest))
    
    # AUC
    H, X1 = np.histogram(ROI, bins = nbins, normed=True)
    dx = X1[1]-X1[0]
    aucCSH = np.round(np.sum(np.cumsum(H)*dx / nbins))

    return HUmax, HUmean, HUmin, HUsd, range, HUmedian, cov, HUpeak, aucCSH


def waveletFilter(rawCT, R):
    
    ratio = R
    coeffs = dwtn(rawCT, 'sym8')

    ratio = 0.5

    weightHHH_LLL = 8/(2*(3*ratio+1))
    weight_Rest = 8*ratio/(2*(3*ratio+1))
    nbcell = len(coeffs)

    contador = 1
    for key in coeffs.keys():
        if contador is 1:
            coeffs[key] = coeffs[key] * weightHHH_LLL 
            contador = contador+1
        elif contador is not nbcell:
            coeffs[key] = coeffs[key] * weight_Rest
            contador = contador+1
        else:
            coeffs[key] = coeffs[key] * weightHHH_LLL

    filteredVol = idwtn(coeffs, 'sym8')

    return filteredVol


def featureExtractionWavelet(rawCT, metaCT, maskPath, R):

    filteredVol = waveletFilter(rawCT, R)

    temp = tempfile.NamedTemporaryFile()
    nrrd.write(temp.name, filteredVol, header=metaCT)
    extractor = featureextractor.RadiomicsFeatureExtractor()
    result = extractor.execute(temp.name, maskPath)
    temp.close()

    dfWavelet = pd.DataFrame.from_dict(dict(result), orient='index')
    dfWavelet.transpose()

    return dfWavelet


def featureExtractionEqHist(rawCT, metaCT, maskPath):
    
    volumeEq = exposure.equalize_hist(rawCT)

    temp = tempfile.NamedTemporaryFile()
    nrrd.write(temp.name, volumeEq, header=metaCT)
    extractor = featureextractor.RadiomicsFeatureExtractor()
    result = extractor.execute(temp.name, maskPath)
    temp.close()

    dfEqHist = pd.DataFrame.from_dict(dict(result), orient='index')
    dfEqHist.transpose()

    return dfEqHist

def volumeErosion(mask, kernel):
    x,y,z = np.shape(mask)
    M_I = np.zeros([x,y,z])

    for slice in np.arange(z):
        M_I[:,:,slice] = binary_erosion(mask[:,:,slice], kernel, iterations=1)
    
    return M_I

def volumeDilation(mask, kernel):
    x,y,z = np.shape(mask)
    M_E = np.zeros([x,y,z])

    for slice in np.arange(z):
        M_E[:,:,slice] = binary_dilation(mask[:,:,slice], kernel, iterations=1)
    
    return M_E

def EstructuringElementDisk(radius):
    kernel = np.zeros([radius*2-1,radius*2-1])
    rest = np.floor(radius/2).astype('int8')
    kernel[rest:-rest,:] = 1 

    counter = 0
    for i in np.arange(rest):
        kernel[i, rest-counter:-(rest-counter)] = 1
        counter = counter+1

    counter = 0
    for i in np.arange(rest):
        kernel[-(i+1), rest-counter:-(rest-counter)] = 1
        counter = counter+1

    return kernel

def decay(rawCT, metaCT, mask, ii2):
    # Parameters
    se_length = np.ceil(np.shape(rawCT)[0]/256*2.5)
    se_length_5 = np.ceil(np.shape(rawCT)[0]/256)
    porc=0.1
    pixelWX = metaCT["space directions"][0,0]
    pixelWY = metaCT["space directions"][1,1]

    # Decay of maximum or minimum intensity
    ini_mask_5 = np.ceil(3.4/pixelWX)
    fin_mask_5 = np.ceil(17/pixelWX)
    ini_mask = np.ceil(23.8/pixelWX)
    fin_mask = np.ceil(37.4/pixelWX)

    kernel = EstructuringElementDisk(se_length.astype('int8'))

    # Inner and outer contour calculation
    M_E = volumeDilation(mask, kernel)
    M_I = volumeErosion(mask, kernel)
    M = M_E - M_I

    # Intensity difference init and mask pixel list
    Diff_C = np.zeros(np.shape(M))
    CC_b_d_r = label(M)
    idLabels = np.argwhere(CC_b_d_r == 1)

    # Saving variables of position (x,y,z) where there is label
    posX = idLabels[:,0]
    posY = idLabels[:,1]
    posZ = idLabels[:,2]

    # Inner and outer contour calculation
    a_aux = rawCT

    # M_E_aux = volumeDilation(mask, EstructuringElementDisk(np.ceil(se_length/2).astype('int8')))
    M_E_aux = volumeDilation(mask, np.ones([np.ceil(se_length).astype('int8'), np.ceil(se_length).astype('int8')]))
    M_E = volumeDilation(mask, EstructuringElementDisk(np.ceil(se_length*1.5).astype('int8'))) - M_E_aux

    if ii2 != 5:
        M_I = volumeErosion(mask, EstructuringElementDisk(se_length.astype('int8')))
    else:
        M_I = volumeErosion(mask, EstructuringElementDisk(se_length_5.astype('int8')))

    Int_value = np.mean(a_aux[np.where(M_I==1)])
    Ext_value = np.mean(a_aux[np.where(M_E==1)])

    # Set the coordinates for profile calculation
    CC_Mask = label(mask)

    I1, I2, I3 = np.where(mask==1)
    inix = np.min(I1)
    finx = np.max(I1)
    meanx = np.round(np.mean(I1)) + 1
    iniy = np.min(I2)
    finy = np.max(I2)
    meany = np.round(np.mean(I2)) + 1
    z = np.unique(I3)

    # meany and meanx have been added one to match the results of matlab

    # Defining reversed volumes for later use
    reversed_a_aux = a_aux[::-1,:,:]
    reversed_a_aux_y = a_aux[:,::-1,:]

    # Define 4 profile at 0º, 90º, 180º and 270º for each slice
    Dist_borde = list([])
    for lin in np.arange(4):
        line = list([])
        for aux in np.arange(len(z)):
            # Checking if insertion is dense bone (little one)
            if ii2 == 5:
                if lin == 0:
                    line.append(a_aux[(inix + ini_mask_5).astype('int32'):(inix + fin_mask_5).astype('int32')+1, meany.astype('int32'), z[aux].astype('int32')])
                elif lin == 1:
                    line.append(reversed_a_aux[(len(reversed_a_aux)-(finx-ini_mask_5)-1).astype('int32'):(len(reversed_a_aux)-(finx-fin_mask_5)).astype('int32')+1, meany.astype('int32'), z[aux].astype('int32')])
                elif lin == 2:
                    line.append(a_aux[meanx.astype('int32'), (iniy+ini_mask_5).astype('int32'):(iniy+fin_mask_5).astype('int32')+1, z[aux].astype('int32')])
                else:
                    line.append(reversed_a_aux_y[meanx.astype('int32'), (len(reversed_a_aux_y)-(finy-ini_mask_5)).astype('int32')-1:(len(reversed_a_aux_y)-(finy-fin_mask_5)).astype('int32'), z[aux].astype('int32')])
            else:
                if lin == 0:
                    line.append(a_aux[(inix + ini_mask).astype('int32'):(inix + fin_mask).astype('int32')+1, meany.astype('int32'), z[aux].astype('int32')])
                elif lin == 1:
                    line.append(reversed_a_aux[(len(reversed_a_aux)-(finx-ini_mask)-1).astype('int32'):(len(reversed_a_aux)-(finx-fin_mask)).astype('int32')+1, meany.astype('int32'), z[aux].astype('int32')])
                elif lin == 2:
                    line.append(a_aux[meanx.astype('int32'), (iniy+ini_mask).astype('int32'):(iniy+fin_mask).astype('int32')+1, z[aux].astype('int32')])
                else:
                    line.append(reversed_a_aux_y[meanx.astype('int32'), (len(reversed_a_aux_y)-(finy-ini_mask)).astype('int32')-1:(len(reversed_a_aux_y)-(finy-fin_mask)).astype('int32'), z[aux].astype('int32')])

        # Mean of all of slices for one direction
        mean_curve = np.mean(line, axis=0)

        # Set maximum and minimum interval intensity value
        Inter = porc*np.abs(Int_value-Ext_value)
        if Int_value > Ext_value:
            Imax = Int_value - Inter
            Imin = Ext_value + Inter
        else:
            Imax = Ext_value - Inter
            Imin = Int_value + Inter

        # Detect position of voxel whitin interval
        pos = list([])
        for aux2 in np.arange(len(mean_curve)):
            if (mean_curve[aux2] > Imin) and (mean_curve[aux2] < Imax):
                pos.append(aux2)

        # Detect maximum number of contiguous voxel positions whithin interval
        cont_tmp = list([])
        cont_final = list([])
        if len(pos) > 1:
            for i in np.arange(len(pos)-1):
                if pos[i+1]-pos[i] == 1:
                    cont_tmp.append(pos[i])
                else:
                    cont_tmp.append(pos[i])
                    if len(cont_tmp) > len(cont_final):
                        cont_final = cont_tmp
                    cont_tmp = list([])
            if len(cont_tmp) > len(cont_final):
                        cont_final = cont_tmp
        elif len(pos) == 1:
            cont_final = pos
        elif len(pos) == 0:
            deriv_pos = list([])
            for i in np.arange(len(mean_curve)-1):
                deriv_pos.append(mean_curve[i+1] - mean_curve[i])
            cont_final = np.argmax(deriv_pos)
            


        return mean_curve, cont_final, Int_value, Ext_value

# Calculates decay of the mean of decays
def contCalc(Int_value, Ext_value, mean_curve):
    porc = 0.1
# Set maximum and minimum interval intensity value
    Inter = porc*np.abs(Int_value-Ext_value)
    if Int_value > Ext_value:
        Imax = Int_value - Inter
        Imin = Ext_value + Inter
    else:
        Imax = Ext_value - Inter
        Imin = Int_value + Inter

    # Detect position of voxel whitin interval
    pos = list([])
    for aux2 in np.arange(len(mean_curve)):
        if (mean_curve[aux2] > Imin) and (mean_curve[aux2] < Imax):
            pos.append(aux2)

    # Detect maximum number of contiguous voxel positions whithin interval
    cont_tmp = list([])
    cont_final = list([])
    if len(pos) > 1:
        for i in np.arange(len(pos)-1):
            if pos[i+1]-pos[i] == 1:
                cont_tmp.append(pos[i])
            else:
                cont_tmp.append(pos[i])
                if len(cont_tmp) > len(cont_final):
                    cont_final = cont_tmp
                cont_tmp = list([])
        if len(cont_tmp) > len(cont_final):
                    cont_final = cont_tmp
    elif len(pos) == 1:
        cont_final = pos
    elif len(pos) == 0:
        deriv_pos = list([])
        for i in np.arange(len(mean_curve)-1):
            deriv_pos.append(mean_curve[i+1] - mean_curve[i])
        cont_final = np.argmax(deriv_pos)
    
    return cont_final
