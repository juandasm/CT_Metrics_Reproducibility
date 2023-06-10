from featureExtraction import doCalibration, doRadiomics, doEdges
from RandRauto import Register
import pandas as pd
import numpy as np
import nrrd
import glob
import os

def Main(pathSave, pathData, position, flag, pathReference):

    if position == 'Center':
        idPosition = 'C'
    elif position == 'Edge':
        idPosition = 'E'

    dirData = os.listdir(pathData)
    os.chdir(pathData)

    for folder in dirData:
        dirData = os.listdir(pathData)
        os.chdir(pathData)
        if os.path.isdir(folder) and folder[0][0] is 'D':
            nameData = folder
            currDir = os.path.join(pathData,nameData)
            print('Current data ' + nameData)
            os.chdir(currDir)
            pathTypes = currDir
            dirTypes = os.listdir(pathTypes)

            for type in dirTypes: 
                if os.path.isdir(type) and type[0][0] is idPosition:
                    nameType = type
                    currDir = os.path.join(pathTypes, type)
                    pathProtocols = currDir
                    print('Current image type ' + nameType)
                    os.chdir(currDir)
                    pathProtocolsRyR = currDir
                    dirProtocols = os.listdir(pathProtocolsRyR)
        
                    # Reading folder of protocols
                    for folder in dirProtocols:
                        if os.path.isdir(os.path.join(pathProtocolsRyR,folder)) and folder[0] is 'P':
                            nameProtocol = folder
                            currDir = os.path.join(pathProtocolsRyR,nameProtocol)
                            print('Current Protocol ' + nameProtocol)
                            os.chdir(currDir)
                            currPathCT = currDir
                            currCTnrrd = glob.glob('*nrrd')
                            
                            # Loop for calibration features extraction
                            if flag == 1 or flag == 4:
                                writer = pd.ExcelWriter(pathSave + nameProtocol + nameData + '_calibration_' + position + '.xlsx', engine = 'xlsxwriter')
                                for CTnrrd in currCTnrrd:
                                    os.chdir(currPathCT)
                                    print('Current CT: ' + CTnrrd)
                                    rawCT, metaCT = nrrd.read(CTnrrd)
                                    
                                    # Reading segmentations
                                    currPathGTV = os.path.join(pathProtocols, 'SegmentationCalibration')
                                    os.chdir(currPathGTV)
                                    dirData = os.listdir()
                                    currGTVs = glob.glob('*.nrrd')
                                    listCalibration = list([])
                                    listContours = list([])
                                    for GTV in currGTVs:
                                        os.chdir(currPathGTV)
                                        print('Current Segmentation: ' + GTV)
                                        # mask, meta_mask = nrrd.read(GTV)
                                        mask = Register(os.path.join(currPathCT, CTnrrd), pathReference, GTV)
                                        contour, calibration = doCalibration(rawCT, metaCT, GTV, nameProtocol, mask)
                                        listCalibration.append(calibration)
                                        listContours.append(contour)
                                    
                                    parametersCalibrationDf = pd.DataFrame()
                                    for i in np.arange(len(currGTVs)):
                                        df = pd.DataFrame.from_records(listCalibration[i], index = [listContours[i]])
                                        df = df[['mean'] + [col for col in df.columns if col != 'mean']]
                                        parametersCalibrationDf = pd.concat([parametersCalibrationDf, df])
                                    parametersCalibrationDf.sort_index
                                    parametersCalibrationDf.to_excel(writer, sheet_name=CTnrrd)
                                writer.save()

                            # Loop for edges features extraction
                            if flag == 2 or flag == 4 and position == 'Center':
                                writer = pd.ExcelWriter(pathSave + nameProtocol + nameData + '_edges.xlsx', engine = 'xlsxwriter')
                                for CTnrrd in currCTnrrd:
                                    os.chdir(currPathCT)
                                    print('Current CT: ' + CTnrrd)
                                    rawCT, metaCT = nrrd.read(CTnrrd)
                                    
                                    # Reading segmentations
                                    currPathGTV = os.path.join(pathProtocols, 'SegmentationEdges')
                                    os.chdir(currPathGTV)
                                    dirData = os.listdir()
                                    currGTVs = glob.glob('*.nrrd')
                                    currGTVs = sorted(currGTVs)
                                    parametersEdgesDf = list([])
                                    counter = 1
                                    listEdges = list([])
                                    listContours = list([])
                                    for GTV in currGTVs:
                                        os.chdir(currPathGTV)
                                        print('Current Segmentation: ' + GTV)
                                        # mask, meta_mask = nrrd.read(GTV)
                                        mask = Register(os.path.join(currPathCT, CTnrrd), pathReference, GTV)
                                        contour, edgesDicc = doEdges(rawCT, metaCT, GTV, nameProtocol, mask, counter)
                                        listEdges.append(edgesDicc)
                                        listContours.append(contour)
                                        counter = counter + 1
                                    
                                    parametersEdgesDf = pd.DataFrame()
                                    for i in np.arange(len(currGTVs)):
                                        df = pd.DataFrame.from_records(listEdges[i], index = [listContours[i]])
                                        parametersEdgesDf = pd.concat([parametersEdgesDf, df])
                                    parametersEdgesDf.sort_index
                                    parametersEdgesDf.to_excel(writer, sheet_name=CTnrrd)
                                writer.save()
                            
                            # Loop for edges features extraction
                            if flag == 3 or flag == 4 and position == 'Center':
                                writer = pd.ExcelWriter(pathSave + nameProtocol + nameData + '_radiomics.xlsx', engine = 'xlsxwriter')
                                for CTnrrd in currCTnrrd:
                                    os.chdir(currPathCT)
                                    print('Current CT: ' + CTnrrd)
                                    rawCT, metaCT = nrrd.read(CTnrrd)
                                    
                                    # Reading segmentations
                                    currPathGTV = os.path.join(pathProtocols, 'SegmentationRadiomics')
                                    os.chdir(currPathGTV)
                                    dirData = os.listdir()
                                    currGTVs = glob.glob('*.nrrd')
                                    currGTVs = sorted(currGTVs)
                                    parametersRadiomicsDf = list([])
                                    counter = 1
                                    listRadiomics = list([])
                                    listContours = list([])
                                    pathCT = os.path.join(currPathCT, CTnrrd)
                                    for GTV in currGTVs:
                                        os.chdir(currPathGTV)
                                        print('Current Segmentation: ' + GTV)
                                        # mask, meta_mask = nrrd.read(GTV)
                                        # pathMask = os.path.join(currPathGTV, GTV)
                                        contour, radiomicsDicc = doRadiomics(pathCT, rawCT, metaCT, nameProtocol, GTV, os.path.join(currPathCT, CTnrrd), pathReference, GTV)
                                        listRadiomics.append(radiomicsDicc)
                                        listContours.append(contour)
                                        counter = counter + 1
                                    
                                    parametersRadiomicsDf = pd.DataFrame()
                                    for i in np.arange(len(currGTVs)):
                                        df = pd.DataFrame.from_records(listRadiomics[i], index = [listContours[i]])
                                        parametersRadiomicsDf = pd.concat([parametersRadiomicsDf, df])
                                    parametersRadiomicsDf.sort_index
                                    parametersRadiomicsDf.to_excel(writer, sheet_name=CTnrrd)
                                writer.save()
