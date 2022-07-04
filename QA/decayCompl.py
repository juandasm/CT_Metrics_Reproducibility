from utilities import decay, contCalc
import matplotlib.pyplot as plt
import numpy as np
import nrrd
import glob
import os

def main(pathSave, pathData, position):

    if position == 'Center':
        idPosition = 'C'
    elif position == 'Edge':
        idPosition = 'E'

    dirData = os.listdir(pathData)
    os.chdir(pathData)
    # listProtocols = ['P04_Radiocirugia', 'P03_Fina', 'P02_Mejorada', 'P01_Normal']
    listProtocols = ['P04_Radiocirugia']
    listGTVs = ['I01In.nrrd', 'I02In.nrrd', 'I03In.nrrd', 'I04In.nrrd', 'I05In.nrrd', 'I06In.nrrd', 'I07In.nrrd', 'I08In.nrrd']
    ii2 = 1

    for protocol in listProtocols:
        for GTV in listGTVs:
            for folder in dirData:
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
                            currDir = os.path.join(pathProtocols, 'Output_R&R')
                            os.chdir(currDir)
                            pathProtocolsRyR = currDir
                            dirProtocols = os.listdir(pathProtocolsRyR)
                
                            # Reading folder of protocols
                            for folderProtocol in dirProtocols:
                                os.chdir(currDir)
                                if os.path.isdir(os.path.join(pathProtocolsRyR,folderProtocol)) and folderProtocol == protocol:
                                    nameProtocol = folderProtocol
                                    currDir = os.path.join(pathProtocolsRyR,nameProtocol)
                                    print('Current Protocol ' + nameProtocol)
                                    os.chdir(currDir)
                                    currPathCT = currDir
                                    currCTnrrd = glob.glob('*nrrd')

                                    currPathGTV = os.path.join(pathProtocols, 'SegmentationEdges')
                                    os.chdir(currPathGTV)
                                    print('Current Segmentation: ' + GTV)
                                    mask, meta_mask = nrrd.read(GTV)

                                    list_means = np.array([])
                                    list_int = np.array([])
                                    list_ext = np.array([])
                                    contador = 1 # Counter to check the number of images
                                    for CTnrrd in currCTnrrd:
                                        os.chdir(currPathCT)
                                        print('Current CT: ' + CTnrrd)
                                        rawCT, metaCT = nrrd.read(CTnrrd)
                                        media, nada, int_value, ext_value = decay(rawCT, metaCT, mask, ii2)
                                        if contador == 1:
                                            list_means = media
                                        else:
                                            list_means = np.append(np.reshape(list_means, [contador-1,len(media)]), np.reshape(media,[1,len(media)]), axis = 0)
                                        list_int = np.append(list_int, int_value)
                                        list_ext = np.append(list_ext, ext_value)
                                        if contador == 4:
                                            break
                                        contador = contador + 1
            mean_int_value = np.mean(list_int)
            mean_ext_value = np.mean(list_ext)
            mean_means = np.mean(list_means, axis = 0)
            cont = contCalc(mean_int_value, mean_ext_value, mean_means)

            plt.plot(mean_means)
            plt.xlabel('Pixels')
            plt.ylabel('Intensity')
            plt.title(GTV + ' Protocol: ' + folderProtocol)
            for i in cont:
                plt.axvspan(i, i+1, facecolor='b', alpha = 0.2)
            plt.savefig(os.path.join(pathSave, folderProtocol, GTV + '.jpg'),dpi=600, format='jpg', bbox_inches = 'tight')
            plt.close()
            ii2 = ii2 + 1


                            # # Loop for edges features extraction
                            # print('Proceeding with edge feature extraction...')
                            # for CTnrrd in currCTnrrd:
                            #     os.chdir(currPathCT)
                            #     print('Current CT: ' + CTnrrd)
                            #     rawCT, metaCT = nrrd.read(CTnrrd)
                                
                            #     # Reading segmentations
                            #     currPathGTV = os.path.join(pathProtocols, 'SegmentationEdges')
                            #     os.chdir(currPathGTV)
                            #     dirData = os.listdir()
                            #     currGTVs = glob.glob('*.nrrd')
                            #     currGTVs = sorted(currGTVs)
                            #     parametersEdgesDf = list([])
                            #     contador = 1
                            #     list_means = np.array([])
                            #     list_int = np.array([])
                            #     list_ext = np.array([])
                            #     for GTV in currGTVs:
                            #         os.chdir(currPathGTV)
                            #         print('Current Segmentation: ' + GTV)
                            #         mask, meta_mask = nrrd.read(GTV)
                            #         media, cont, int_value, ext_value = decay(rawCT, metaCT, mask, contador)
                            #         contador = contador + 1
                            #         list_means = np.append(list_means, media)
                            #         list_int = np.append(list_int, int_value)
                            #         list_ext = np.append(list_ext)


                                # plt.plot(media)
                                # plt.xlabel('Pixels')
                                # plt.ylabel('Intensity')
                                # plt.title(GTV + ' Protocol: ' + protocol)
                                # for i in cont:
                                #     plt.axvspan(i, i+1, facecolor='b', alpha = 0.2)
                                # plt.savefig(os.path.join(pathSave, protocol, GTV + '.jpg'),dpi=600, format='jpg', bbox_inches = 'tight')
                                # plt.close()
