from utilities import getSUVMetricP, getPercentInactive, getEccentricity, EstructuringElementDisk, volumeDilation, volumeErosion
from radiomics import featureextractor
from RandRauto import Register
from skimage.measure import label
import pandas as pd
import numpy as np
import nrrd
import os

def doCalibration(rawCT, metaCT, GTV, protocol, mask):

    pixelW = metaCT["space directions"][0,0]
    slice_S = metaCT["space directions"][2,2]
    
    # nbins = 64
    # volume = rawCT
    # R = 1
    # scale = 1
    
    HUmax, HUmean, HUmin, HUsd, range, HUmedian, cov, HUpeak, aucCSH = getSUVMetricP(rawCT, mask) 

    ParametersCT = {
    'protocol':protocol, 
    'mean':HUmean,
    'SD':HUsd,
    'CoV':cov,
    'max':HUmax,
    'min':HUmin,
    'range':range,
    'median':HUmedian,
    'SUVpeak':HUpeak
    }

    return GTV, ParametersCT

def doRadiomics(pathCT, rawCT, metaCT, protocol, GTV, pathFixed, pathMoving, pathSegmentation):

    mask, pathMask_tmp = Register(pathFixed, pathMoving, pathSegmentation, 1)
    pathMask = os.path.join(pathMask_tmp.name, 'mask.nrrd')

    # Prueba tamaño
    testMask, testMeta = nrrd.read(pathMask)
    print(np.shape(testMask))
    print(np.shape(rawCT))

    HUmax, HUmean, HUmin, HUsd, range, HUmedian, cov, HUpeak, aucCSH = getSUVMetricP(rawCT, mask) 
    
    extractor = featureextractor.RadiomicsFeatureExtractor()
    result = extractor.execute(pathCT, pathMask)

    df = pd.DataFrame.from_dict(dict(result), orient='index')
    df.transpose()

    R = 0.5

    # dfWavelet = featureExtractionWavelet(rawCT, metaCT, pathMask, R)
    # dfEqHist = featureExtractionEqHist(rawCT, metaCT, pathMask)

    # Habría que ponerlo df.loc['original_shapeMeshVolume'][0]
    volume = df.loc['original_shape_MeshVolume'][0]
    HUxV = volume*HUmean
    if HUmax + HUmin == 0:
        integralUniformity = 0
    else:
        integralUniformity = (HUmax-HUmin)/(HUmax+HUmin)
    # print('Getting differential uniformity...')   De momento lo quito porque me quema el ordenador
    # differentialUniformity = getDifferentialUniformity(rawCT)
    # print('Finished getting diferential uniformity')
    skewness = df.loc['original_firstorder_Skewness'][0]
    kurtosis = df.loc['original_firstorder_Kurtosis'][0]
    entropy = df.loc['original_firstorder_Entropy'][0]
    energy = df.loc['original_firstorder_Energy'][0]
    # Necesario chequear solidity, no da lo mismo que en matlab, sacado de paper https://arxiv.org/pdf/1612.07003.pdf
    solidity = df.loc['original_shape_VoxelVolume'][0] / df.loc['original_shape_MeshVolume'][0]
    percentInactive = getPercentInactive(rawCT, mask) # Default thresh = 0.6
    eccentricity = getEccentricity(rawCT, metaCT, mask)
    glcmEnergy = df.loc['original_glcm_JointEnergy'][0]
    glcmContrast = df.loc['original_glcm_Contrast'][0]
    glcmEntropy = df.loc['original_glcm_JointEnergy'][0]
    glcmHomogeinity = df.loc['original_glcm_Idm'][0]
    glcmCorrelation = df.loc['original_glcm_Correlation'][0]
    glcmVariance = df.loc['original_glcm_DifferenceVariance'][0] #Chequear
    glcmDissimilarity = df.loc['original_glcm_DifferenceAverage'][0]
    glcmAutocorrelation = df.loc['original_glcm_Autocorrelation'][0]

    glszmSAE = df.loc['original_glszm_SmallAreaEmphasis'][0]
    glszmLAE = df.loc['original_glszm_LargeAreaEmphasis'][0]
    glszmGLN = df.loc['original_glszm_GrayLevelNonUniformity'][0]
    glszmSZN = df.loc['original_glszm_SizeZoneNonUniformity'][0]
    glszmZP = df.loc['original_glszm_ZonePercentage'][0]
    glszmLGLZE = df.loc['original_glszm_LowGrayLevelZoneEmphasis'][0]
    glszmHGLZE = df.loc['original_glszm_HighGrayLevelZoneEmphasis'][0]
    glszmSALGLE = df.loc['original_glszm_SmallAreaLowGrayLevelEmphasis'][0]
    glszmSAHGLE = df.loc['original_glszm_SmallAreaHighGrayLevelEmphasis'][0]
    glszmLALGLE = df.loc['original_glszm_LargeAreaLowGrayLevelEmphasis'][0]
    glszmLAHGLE = df.loc['original_glszm_LargeAreaHighGrayLevelEmphasis'][0]
    glszmGLV = df.loc['original_glszm_GrayLevelVariance'][0]
    glszmZV = df.loc['original_glszm_ZoneVariance'][0]

    glrlmSRE = df.loc['original_glrlm_ShortRunEmphasis'][0]
    glrlmLRE = df.loc['original_glrlm_LongRunEmphasis'][0]
    glrlmGLN = df.loc['original_glrlm_GrayLevelNonUniformity'][0]
    glrlmRLN = df.loc['original_glrlm_RunLengthNonUniformity'][0]
    glrlmRP = df.loc['original_glrlm_RunPercentage'][0]
    glrlmLGRE = df.loc['original_glrlm_LowGrayLevelRunEmphasis'][0]
    glrlmHGRE = df.loc['original_glrlm_HighGrayLevelRunEmphasis'][0]
    glrlmSRLGLE = df.loc['original_glrlm_ShortRunLowGrayLevelEmphasis'][0]
    glrlmSRHGLE = df.loc['original_glrlm_ShortRunHighGrayLevelEmphasis'][0]
    glrlmLRLGLE = df.loc['original_glrlm_LongRunLowGrayLevelEmphasis'][0]
    glrlmLRHGLE = df.loc['original_glrlm_LongRunHighGrayLevelEmphasis'][0]
    glrlmGLV = df.loc['original_glrlm_GrayLevelVariance'][0]
    glrlmRLV = df.loc['original_glrlm_RunVariance'][0]

    ngtdmCoarseness = df.loc['original_ngtdm_Coarseness'][0]
    ngtdmContrast = df.loc['original_ngtdm_Contrast'][0]
    ngtdmBusyness = df.loc['original_ngtdm_Busyness'][0]
    ngtdmComplexity = df.loc['original_ngtdm_Complexity'][0]
    ngtdmStrength = df.loc['original_ngtdm_Strength'][0]

    radiomics = {'Protocol': protocol, 'GTV':GTV, 'volume':volume, 'HUxV':HUxV, 'IntergalUniformity':integralUniformity, 'Skewness':skewness, 'Kurtosis':kurtosis, 'Entropy':entropy, 
    'Energy':energy, 'Solidity':solidity, 'PercentInactive':percentInactive, 'Eccentricity':eccentricity, 'GLCMEnergy':glcmEnergy, 'GLCMContrast':glcmContrast, 'GLCMEntropy':glcmEntropy,
    'GLCMHomogeinity':glcmHomogeinity, 'GLCMCorrelation':glcmCorrelation, 'GLCMVariance':glcmVariance, 'GLCMDissimilarity':glcmDissimilarity, 'GLCMAutocorrelation':glcmAutocorrelation, 
    'GLSZMSAE':glszmSAE, 'GLSZMLAE':glszmLAE, 'GLSZMGLN':glszmGLN, 'GLSZMSZN':glszmSZN, 'GLSZMZP':glszmZP, 'GLSZMLGLZE':glszmLGLZE, 'GLSZMHGLZE': glszmHGLZE, 'GLSZMSALGLE':glszmSALGLE, 
    'GLSZMSAHGLE':glszmSAHGLE, 'GLSZMLALGLE':glszmLALGLE, 'GLSZMLAHGLE':glszmLAHGLE, 'GLSZMGLV':glszmGLV, 'GLSZMZV':glszmZV, 'GLRLMSRE':glrlmSRE, 'GLRLMLRE':glrlmLRE, 'GLRLMGLN': glrlmGLN,
    'GLRLMRLN': glrlmRLN, 'GLRLMRP':glrlmRP, 'GLRLMLGRE':glrlmLGRE, 'GLRLMHGRE': glrlmHGRE, 'GLRLMSRLGLE':glrlmSRLGLE, 'GLRLMSRHGLE':glrlmSRHGLE, 'GLRLMLRLGLE':glrlmLRLGLE, 'GLRLMLRHGLE':glrlmLRHGLE,
    'GLRLMGLV':glrlmGLV, 'GLRLMRLV':glrlmRLV, 'NGTDMCoarseness':ngtdmCoarseness, 'NGTDMContrast':ngtdmContrast, 'NGTDMBusyness':ngtdmBusyness, 'NGTDMComplexity':ngtdmComplexity, 'NGTDMStrength':ngtdmStrength}

    # Radiomics for wavelet filtered image

    # wglcmEnergy = dfWavelet.loc['original_glcm_JointEnergy']
    # wglcmContrast = dfWavelet.loc['original_glcm_Contrast']
    # wglcmEntropy = dfWavelet.loc['original_glcm_JointEnergy']
    # wglcmHomogeinity = dfWavelet.loc['original_glcm_Idm']
    # wglcmCorrelation = dfWavelet.loc['original_glcm_Correlation']
    # wglcmVariance = dfWavelet.loc['original_glcm_DifferenceVariance'] #Chequear
    # wglcmDissimilarity = dfWavelet.loc['original_glcm_DifferenceAverage']
    # wglcmAutocorrelation = dfWavelet.loc['original_glcm_Autocorrelation']

    # wglszmSAE = dfWavelet.loc['original_glszm_SmallAreaEmphasis']
    # wglszmLAE = dfWavelet.loc['original_glszm_LargeAreaEmphasis']
    # wglszmGLN = dfWavelet.loc['original_glszm_GrayLevelNonUniformity']
    # wglszmSZN = dfWavelet.loc['original_glszm_SizeZoneNonUniformity']
    # wglszmZP = dfWavelet.loc['original_glszm_ZonePercentage']
    # wglszmLGLZE = dfWavelet.loc['original_glszm_LowGrayLevelZoneEmphasis']
    # wglszmHGLZE = dfWavelet.loc['original_glszm_HighGrayLevelZoneEmphasis']
    # wglszmSALGLE = dfWavelet.loc['original_glszm_SmallAreaLowGrayLevelEmphasis']
    # wglszmSAHGLE = dfWavelet.loc['original_glszm_SmallAreaHighGrayLevelEmphasis']
    # wglszmLALGLE = dfWavelet.loc['original_glszm_LargeAreaLowGrayLevelEmphasis']
    # wglszmLAHGLE = dfWavelet.loc['original_glszm_LargeAreaHighGrayLevelEmphasis']
    # wglszmGLV = dfWavelet.loc['original_glszm_GrayLevelVariance'] 
    # wglszmZV = dfWavelet.loc['original_glszm_ZoneVariance']   

    # wglrlmSRE = dfWavelet.loc['original_glrlm_ShortRunEmphasis']
    # wglrlmLRE = dfWavelet.loc['original_glrlm_LongRunEmphasis']
    # wglrlmGLN = dfWavelet.loc['original_glrlm_GrayLevelNonUniformity']
    # wglrlmRLN = dfWavelet.loc['original_glrlm_RunLengthNonUniformity']
    # wglrlmRP = dfWavelet.loc['original_glrlm_RunPercentage']
    # wglrlmLGRE = dfWavelet.loc['original_glrlm_LowGrayLevelRunEmphasis']
    # wglrlmHGRE = dfWavelet.loc['original_glrlm_HighGrayLevelRunEmphasis']
    # wglrlmSRLGLE = dfWavelet.loc['original_glrlm_ShortRunLowGrayLevelEmphasis']
    # wglrlmSRHGLE = dfWavelet.loc['original_glrlm_ShortRunHighGrayLevelEmphasis']
    # wglrlmLRLGLE = dfWavelet.loc['original_glrlm_LongRunLowGrayLevelEmphasis']
    # wglrlmLRHGLE = dfWavelet.loc['original_glrlm_LongRunHighGrayLevelEmphasis']
    # wglrlmGLV = dfWavelet.loc['original_glrlm_GrayLevelVariance']
    # wglrlmRLV = dfWavelet.loc['original_glrlm_RunVariance']

    # wngtdmCoarseness = dfWavelet.loc['original_ngtdm_Coarseness']
    # wngtdmContrast = dfWavelet.loc['original_ngtdm_Contrast']
    # wngtdmBusyness = dfWavelet.loc['original_ngtdm_Busyness']
    # wngtdmComplexity = dfWavelet.loc['original_ngtdm_Complexity']
    # wngtdmStrength = dfWavelet.loc['original_ngtdm_Strength']

    # # Radiomics for equalized histogram image

    # qglcmEnergy = dfEqHist.loc['original_glcm_JointEnergy']
    # qglcmContrast = dfEqHist.loc['original_glcm_Contrast']
    # qglcmEntropy = dfEqHist.loc['original_glcm_JointEnergy']
    # qglcmHomogeinity = dfEqHist.loc['original_glcm_Idm']
    # qglcmCorrelation = dfEqHist.loc['original_glcm_Correlation']
    # qglcmVariance = dfEqHist.loc['original_glcm_DifferenceVariance'] #Chequear
    # qglcmDissimilarity = dfEqHist.loc['original_glcm_DifferenceAverage']
    # qglcmAutocorrelation = dfEqHist.loc['original_glcm_Autocorrelation']

    # qglszmSAE = dfEqHist.loc['original_glszm_SmallAreaEmphasis']
    # qglszmLAE = dfEqHist.loc['original_glszm_LargeAreaEmphasis']
    # qglszmGLN = dfEqHist.loc['original_glszm_GrayLevelNonUniformity']
    # qglszmSZN = dfEqHist.loc['original_glszm_SizeZoneNonUniformity']
    # qglszmZP = dfEqHist.loc['original_glszm_ZonePercentage']
    # qglszmLGLZE = dfEqHist.loc['original_glszm_LowGrayLevelZoneEmphasis']
    # qglszmHGLZE = dfEqHist.loc['original_glszm_HighGrayLevelZoneEmphasis']
    # qglszmSALGLE = dfEqHist.loc['original_glszm_SmallAreaLowGrayLevelEmphasis']
    # qglszmSAHGLE = dfEqHist.loc['original_glszm_SmallAreaHighGrayLevelEmphasis']
    # qglszmLALGLE = dfEqHist.loc['original_glszm_LargeAreaLowGrayLevelEmphasis']
    # qglszmLAHGLE = dfEqHist.loc['original_glszm_LargeAreaHighGrayLevelEmphasis']
    # qglszmGLV = dfEqHist.loc['original_glszm_GrayLevelVariance'] 
    # qglszmZV = dfEqHist.loc['original_glszm_ZoneVariance']   

    # qglrlmSRE = dfEqHist.loc['original_glrlm_ShortRunEmphasis']
    # qglrlmLRE = dfEqHist.loc['original_glrlm_LongRunEmphasis']
    # qglrlmGLN = dfEqHist.loc['original_glrlm_GrayLevelNonUniformity']
    # qglrlmRLN = dfEqHist.loc['original_glrlm_RunLengthNonUniformity']
    # qglrlmRP = dfEqHist.loc['original_glrlm_RunPercentage']
    # qglrlmLGRE = dfEqHist.loc['original_glrlm_LowGrayLevelRunEmphasis']
    # qglrlmHGRE = dfEqHist.loc['original_glrlm_HighGrayLevelRunEmphasis']
    # qglrlmSRLGLE = dfEqHist.loc['original_glrlm_ShortRunLowGrayLevelEmphasis']
    # qglrlmSRHGLE = dfEqHist.loc['original_glrlm_ShortRunHighGrayLevelEmphasis']
    # qglrlmLRLGLE = dfEqHist.loc['original_glrlm_LongRunLowGrayLevelEmphasis']
    # qglrlmLRHGLE = dfEqHist.loc['original_glrlm_LongRunHighGrayLevelEmphasis']
    # qglrlmGLV = dfEqHist.loc['original_glrlm_GrayLevelVariance']
    # qglrlmRLV = dfEqHist.loc['original_glrlm_RunVariance']

    # qngtdmCoarseness = dfEqHist.loc['original_ngtdm_Coarseness']
    # qngtdmContrast = dfEqHist.loc['original_ngtdm_Contrast']
    # qngtdmBusyness = dfEqHist.loc['original_ngtdm_Busyness']
    # qngtdmComplexity = dfEqHist.loc['original_ngtdm_Complexity']
    # qngtdmStrength = dfEqHist.loc['original_ngtdm_Strength']

    pathMask_tmp.cleanup()

    return GTV, radiomics
    
def doEdges(rawCT, metaCT, NameContour, NamePatient, mask, ii2):

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

    # Intensity difference calculation
    for i in np.arange(np.shape(idLabels)[0]):
        x = posX[i]
        y = posY[i]
        z = posZ[i]

        if rawCT[x+1,y,z] == 0:
            a1 = 0
        else:
            a1 = (rawCT[x,y,z] - rawCT[x+1,y,z])/pixelWX
        
        if rawCT[x-1,y,z] == 0:
            a2 = 0
        else:
            a2 = (rawCT[x,y,z] - rawCT[x-1,y,z])/pixelWX
        
        if rawCT[x,y+1,z] == 0:
            a3 = 0
        else:
            a3 = (rawCT[x,y,z] - rawCT[x,y+1,z])/pixelWY
        
        if rawCT[x,y-1,z] == 0:
            a4 = 0
        else:
            a4 = (rawCT[x,y,z] - rawCT[x,y-1,z])/pixelWY
        
        if rawCT[x+1,y+1,z] == 0:
            a5 = 0
        else:
            a5 = (rawCT[x,y,z] - rawCT[x+1,y+1,z])/np.sqrt(2)*pixelWY
        
        if rawCT[x-1,y+1,z] == 0:
            a6 = 0
        else:
            a6 = (rawCT[x,y,z] - rawCT[x-1,y+1,z])/np.sqrt(2)*pixelWY
        
        if rawCT[x-1,y-1,z] == 0:
            a7 = 0
        else:
            a7 = (rawCT[x,y,z] - rawCT[x-1,y-1,z])/np.sqrt(2)*pixelWY
        
        if rawCT[x+1,y-1,z] == 0:
            a8 = 0
        else:
            a8 = (rawCT[x,y,z] - rawCT[x+1,y-1,z])/np.sqrt(2)*pixelWY

        Diff_C[x,y,z] = (np.mean(np.abs([a1,a2,a3,a4,a5,a6,a7,a8])))**2

    # Final Tissue Contrast metrics
    TC_mean = np.round(np.mean(Diff_C[(M==1)]))

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
        cont_final = list([])
        if len(pos) > 1:
            for aux in np.arange(len(pos)-1):
                cont = 1
                for aux2 in np.arange(aux+1,len(pos),1):
                    if pos[aux2]-pos[aux2-1] != 1:
                        break
                    cont = cont+1
                cont_final.append(cont)
            num_voxels = np.max(cont_final)
        elif len(pos) == 1:
            num_voxels = 1
        else:
            num_voxels = 0

        if num_voxels != 0:
            Dist_borde.append(1/(pixelWX * num_voxels))
        else:
            Dist_borde.append(0)

    # Final Drop Rage metrics
    DR_mean = np.round(np.mean(Dist_borde), decimals=2)

    ParametersEdges = {
    'protocol':NamePatient, 
    'TC_mean':TC_mean,
    'DR_mean':DR_mean
    }

    return NameContour, ParametersEdges
