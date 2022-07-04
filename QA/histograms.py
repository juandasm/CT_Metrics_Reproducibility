from sklearn.metrics import r2_score
from scipy.stats import shapiro, ranksums
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np 
import glob
import os

def Histograms(pathOutput, pathSave, position, mode, flagMeasure):

    if mode == 'CalibrationOutput':
        pathOutput = os.path.join(pathOutput, mode, position)

        material_index = {'In_Trabecular_Bone':0, 
                    'Out_Trabecular_Bone':1, 
                    'In_Breast':2, 
                    'Out_Breast':3,
                    'In_Muscle':4,
                    'Out_Muscle':5,
                    'In_Adipose':6,
                    'Out_Adipose':7,
                    'In_Dense_Bone':8,
                    'Out_Dense_Bone':9,
                    'In_Lung_Exhale':10,
                    'Out_Lung_Exhale':11,
                    'In_Lung_Inhale':12,
                    'Out_Lung_Inhale':13,
                    'In_Liver':14,
                    'Out_Liver':15,
                    'Water':16}
            
    if mode == 'EdgesOutput':

        pathOutput = os.path.join(pathOutput, mode, position)

        material_index = {
            'I01In.nrrd':0,
            'I02In.nrrd':1,
            'I03In.nrrd':2,
            'I04In.nrrd':3,
            'I05In.nrrd':4,
            'I06In.nrrd':5,
            'I07In.nrrd':6,
            'I08In.nrrd':7,
        }

    listProtocols = ['Radiocirugia', 'Fina', 'Mejorada', 'Normal']
    listMeasures = ['mean', 'CoV', 'SD', 'SUVpeak', 'max', 'min', 'median', 'range', 'DR_mean', 'TC_mean']

    flagMeasure = flagMeasure - 1
    measure = listMeasures[flagMeasure]

    for protocol in listProtocols:
        listDfsRadiomics = list([])

        pathOutputCal = os.path.join(pathOutput, protocol)
        os.chdir(pathOutputCal)
        currOutputs = glob.glob('*xlsx')
        for excel in currOutputs:
            for i in np.arange(5):
                df = pd.read_excel(os.path.join(pathOutputCal, excel), index_col=0, sheet_name=i)
                df = df.sort_index()
                listDfsRadiomics.append(df)

        mean_material = np.array([])
        normalFit = {}
        normalTest = {}

        for material, index in material_index.items():
            for df in listDfsRadiomics:
                mean_material = np.append(mean_material,(df[measure].iloc[material_index[material]]))
            sorted_material = np.sort(mean_material)

            stat, p = shapiro(mean_material)
            # print('Statistics=%.3f, p=%.3f' % (stat, p))
            # interpret
            alpha = 0.05
            if p > alpha:
                # print('Sample ' + material + ' looks Gaussian (fail to reject H0)')
                normalFit.update({material:p})
                normalTest.update({material:'looks gaussian'})
            else:
                # print('Sample ' + material + ' does not look Gaussian (reject H0)')
                normalFit.update({material:p})
                normalTest.update({material:'does not look gaussian'})
            mean_value = np.round(np.mean(sorted_material))
            plt.hist(sorted_material, bins=8)
            if measure is 'mean':
                plt.xlabel('Hounsfield units')
            else:
                plt.xlabel(measure)
            plt.ylabel('Frecuency')
            plt.title('Material: ' + material + ' Protocol: ' + protocol + '. Mean = ' + str(mean_value))
            plt.savefig(os.path.join(pathSave,position,measure,protocol, material + '.jpg'),dpi=600, format='jpg', bbox_inches = 'tight')
            plt.close()

            mean_material = np.array([])


def Calibracion(pathOutput, pathSave, position, loc):
    # Representación HU vs Densidad

    material_index = {'In_Trabecular_Bone':0, 
                    'Out_Trabecular_Bone':1, 
                    'In_Breast':2, 
                    'Out_Breast':3,
                    'In_Muscle':4,
                    'Out_Muscle':5,
                    'In_Adipose':6,
                    'Out_Adipose':7,
                    'In_Dense_Bone':8,
                    'Out_Dense_Bone':9,
                    'In_Lung_Exhale':10,
                    'Out_Lung_Exhale':11,
                    'In_Lung_Inhale':12,
                    'Out_Lung_Inhale':13,
                    'In_Liver':14,
                    'Out_Liver':15,
                    'Water':16}

    material_index_in = ['In_Trabecular_Bone',
                'In_Breast',
                'In_Muscle',
                'In_Adipose',
                'In_Dense_Bone',
                'In_Lung_Exhale',
                'In_Lung_Inhale',
                'In_Liver',
                'Water']

    material_index_out = ['Out_Trabecular_Bone',
                    'Out_Breast',
                    'Out_Muscle',
                    'Out_Adipose',
                    'Out_Dense_Bone',
                    'Out_Lung_Exhale',
                    'Out_Lung_Inhale',
                    'Out_Liver']
    
    if loc == 'In':
        material_index_pos = material_index_in
    elif loc == 'Out':
        material_index_pos = material_index_out

    RED = {'In_Trabecular_Bone':1.117, 
        'Out_Trabecular_Bone':1.117, 
        'In_Breast':0.976, 
        'Out_Breast':0.976,
        'In_Muscle':1.043,
        'Out_Muscle':1.043,
        'In_Adipose':0.949,
        'Out_Adipose':0.949,
        'In_Dense_Bone':1.695,
        'Out_Dense_Bone':1.695,
        'In_Lung_Exhale':0.496,
        'Out_Lung_Exhale':0.496,
        'In_Lung_Inhale':0.200,
        'Out_Lung_Inhale':0.200,
        'In_Liver':1.052,
        'Out_Liver':1.052,
        'Water':0.998}

    # RED = {'In_Trabecular_Bone':1.16, 
    #     'Out_Trabecular_Bone':1.16, 
    #     'In_Breast':0.99, 
    #     'Out_Breast':0.99,
    #     'In_Muscle':1.06,
    #     'Out_Muscle':1.06,
    #     'In_Adipose':0.97,
    #     'Out_Adipose':0.97,
    #     'In_Dense_Bone':1.61,
    #     'Out_Dense_Bone':1.61,
    #     'In_Lung_Exhale':0.5,
    #     'Out_Lung_Exhale':0.5,
    #     'In_Lung_Inhale':0.200,
    #     'Out_Lung_Inhale':0.200,
    #     'In_Liver':1.07,
    #     'Out_Liver':1.07,
    #     'Water':0.998}

    listProtocols = ['Radiocirugia', 'Fina', 'Mejorada', 'Normal']
    pathOutput = os.path.join('/Volumes/T7/iMac/Output_full/CalibrationOutput', position)
    measure = 'mean'
    for protocol in listProtocols:
        listDfsRadiomics = list([])

        pathOutputCal = os.path.join(pathOutput, protocol)
        os.chdir(pathOutputCal)
        currOutputs = glob.glob('*xlsx')
        for excel in currOutputs:
            for i in np.arange(5):
                df = pd.read_excel(os.path.join(pathOutputCal, excel), index_col=0, sheet_name=i)
                df = df.sort_index()
                listDfsRadiomics.append(df)

        mean_material = np.array([])
        normalFit = {}
        normalTest = {}

        list_means = np.array([])
        list_densities = np.array([])
        for material in material_index_pos:
            for df in listDfsRadiomics:
                mean_material = np.append(mean_material,(df[measure].iloc[material_index[material]]))
            list_means = np.append(list_means,np.mean(mean_material))
            list_densities = np.append(list_densities, RED[material])
            mean_material = np.array([])
        
        # Regression equation
        regresion = np.polyfit(np.sort(list_densities), np.sort(list_means), 1)
        x = np.linspace((np.min(list_densities)), (np.max(list_densities)), 2000)
        y = regresion[0]*x + regresion[1]

        # R2
        GT = regresion[0]*np.sort(list_densities) + regresion[1]
        R2 = r2_score(np.sort(list_means), GT)
        R2 = np.round(R2, decimals=4)
        
        plt.plot(np.sort(list_densities), np.sort(list_means), 'bo')
        plt.plot(x,y, 'g--', alpha=0.5)
        if regresion[1] > 0:
            plt.title('Regresión lineal: ' + str(np.round(regresion[0]).astype('int16')) + 'x + ' + str(np.round(regresion[1]).astype('int16')) + '. R square value = ' + str(R2))
        else:
            plt.title('Regresión lineal: ' + str(np.round(regresion[0]).astype('int16'))  + 'x - ' + str(-np.round(regresion[1]).astype('int16'))  + '. R square value = ' + str(R2))
        plt.xlabel('RED (Relative electron density to Water)')
        plt.ylabel('Hounsfield Units')
        # plt.savefig('/Users/juansaboridomoral/Desktop/' + protocol + '.jpg', dpi=600, format='jpg', bbox_inches = 'tight')
        plt.savefig(os.path.join(pathSave, position, 'Calibration',loc, protocol + '.jpg'), dpi = 600, format = 'jpg', bbox_inches = 'tight')
        plt.close()


def NormalityTest(pathOutput, pathSave, position, mode):
    ## Código para analizar si los datos se ajustan a una distribución normal

    listProtocols = ['Radiocirugia', 'Fina', 'Mejorada', 'Normal']

    if mode == 'CalibrationOutput':
        # En caso de calibracion
        listMeasures = ['mean', 'CoV', 'SD', 'SUVpeak', 'max', 'min', 'median', 'range']
        material_index = {'In_Trabecular_Bone':0, 
                    'Out_Trabecular_Bone':1, 
                    'In_Breast':2, 
                    'Out_Breast':3,
                    'In_Muscle':4,
                    'Out_Muscle':5,
                    'In_Adipose':6,
                    'Out_Adipose':7,
                    'In_Dense_Bone':8,
                    'Out_Dense_Bone':9,
                    'In_Lung_Exhale':10,
                    'Out_Lung_Exhale':11,
                    'In_Lung_Inhale':12,
                    'Out_Lung_Inhale':13,
                    'In_Liver':14,
                    'Out_Liver':15,
                    'Water':16}
        saveFolder = 'Calibration'

    if mode == 'EdgesOutput':
        # En caso de edge
        listMeasures = ['DR_mean', 'TC_mean']
        material_index = {
        'I01In.nrrd':0,
        'I02In.nrrd':1,
        'I03In.nrrd':2,
        'I04In.nrrd':3,
        'I05In.nrrd':4,
        'I06In.nrrd':5,
        'I07In.nrrd':6,
        'I08In.nrrd':7,
        }
        saveFolder = 'Edges'
    
    if mode == 'RadiomicsOutput':
        # En caso de radiómica
        listMeasures = ['volume', 'HUxV', 'IntergalUniformity', 'Skewness', 'Kurtosis', 'Entropy', 'Energy', 'Solidity', 'PercentInactive', 'Eccentricity', 'GLCMEnergy', 'GLCMContrast', 'GLCMEntropy', 'GLCMHomogeinity', 'GLCMCorrelation', 'GLCMVariance', 'GLCMDissimilarity', 'GLCMAutocorrelation', 'GLSZMSAE', 'GLSZMLAE', 'GLSZMGLN', 'GLSZMSZN', 'GLSZMZP', 'GLSZMLGLZE', 'GLSZMHGLZE', 'GLSZMSALGLE', 'GLSZMSAHGLE', 'GLSZMLALGLE', 'GLSZMLAHGLE', 'GLSZMGLV', 'GLSZMZV', 'GLRLMSRE', 'GLRLMLRE', 'GLRLMGLN', 'GLRLMRLN', 'GLRLMRP', 'GLRLMLGRE', 'GLRLMHGRE', 'GLRLMSRLGLE', 'GLRLMSRHGLE', 'GLRLMLRLGLE', 'GLRLMLRHGLE', 'GLRLMGLV', 'GLRLMRLV', 'NGTDMCoarseness', 'NGTDMContrast', 'NGTDMBusyness', 'NGTDMComplexity', 'NGTDMStrength']
        material_index = {
            'H1.nrrd':0,
            'H2.nrrd':1,
            'H3.nrrd':2,
            'H4.nrrd':3,
            'H5.nrrd':4,
            'H6.nrrd':5,
            'H7.nrrd':6,
            'H8.nrrd':7,
            'H9.nrrd':8
        }
        saveFolder = 'Radiomics'

    for measure in listMeasures:
        listTests = list([])
        for protocol in listProtocols:
            listDfsRadiomics = list([])
            pathOutputCal = os.path.join(pathOutput, mode, position, protocol)
            os.chdir(pathOutputCal)
            currOutputs = glob.glob('*xlsx')
            for excel in currOutputs:
                for i in np.arange(5):
                    df = pd.read_excel(os.path.join(pathOutputCal, excel), index_col=0, sheet_name=i)
                    df = df.sort_index()
                    listDfsRadiomics.append(df)

            mean_material = np.array([])
            normalFit = {}
            for material, index in material_index.items():
                for df in listDfsRadiomics:
                    mean_material = np.append(mean_material,(df[measure].iloc[material_index[material]]))
                sorted_material = np.sort(mean_material)
                stat, p = shapiro(mean_material)
                # print('Statistics=%.3f, p=%.3f' % (stat, p))
                # interpret
                alpha = 0.05
                if p > alpha:
                    # print('Sample ' + material + ' looks Gaussian (fail to reject H0)')
                    normalFit.update({material:p})
                else:
                    # print('Sample ' + material + ' does not look Gaussian (reject H0)')
                    normalFit.update({material:p})

                mean_material = np.array([])

            listTests.append(normalFit)

        testDf = pd.DataFrame()
        for i in np.arange(len(listTests)):
            df = pd.DataFrame.from_dict(listTests[i], orient='index', columns = [listProtocols[i]])
            testDf = pd.concat([testDf, df], axis = 1)
        testDf.sort_index
        testDf.to_excel(os.path.join(pathSave, saveFolder, position, measure + '_' + position + '.xlsx'))


def WRST(pathOutput, pathSave, mode, position):

    if mode == 'CalibrationOutput':
        material_index = {'In_Trabecular_Bone':0, 
                        'Out_Trabecular_Bone':1, 
                        'In_Breast':2, 
                        'Out_Breast':3,
                        'In_Muscle':4,
                        'Out_Muscle':5,
                        'In_Adipose':6,
                        'Out_Adipose':7,
                        'In_Dense_Bone':8,
                        'Out_Dense_Bone':9,
                        'In_Lung_Exhale':10,
                        'Out_Lung_Exhale':11,
                        'In_Lung_Inhale':12,
                        'Out_Lung_Inhale':13,
                        'In_Liver':14,
                        'Out_Liver':15,
                        'Water':16}

        listMeasures = ['mean']
        saveFolder = 'Calibration'

    if mode == 'EdgesOutput':

        material_index = {
            'I01In.nrrd':0,
            'I02In.nrrd':1,
            'I03In.nrrd':2,
            'I04In.nrrd':3,
            'I05In.nrrd':4,
            'I06In.nrrd':5,
            'I07In.nrrd':6,
            'I08In.nrrd':7,
        }

        listMeasures = ['DR_mean', 'TC_mean']
        saveFolder = 'Edges'

    if mode == 'RadiomicsOutput':
        # En caso de radiómica
        listMeasures = ['volume', 'HUxV', 'IntergalUniformity', 'Skewness', 'Kurtosis', 'Entropy', 'Energy', 'Solidity', 'PercentInactive', 'Eccentricity', 'GLCMEnergy', 'GLCMContrast', 'GLCMEntropy', 'GLCMHomogeinity', 'GLCMCorrelation', 'GLCMVariance', 'GLCMDissimilarity', 'GLCMAutocorrelation', 'GLSZMSAE', 'GLSZMLAE', 'GLSZMGLN', 'GLSZMSZN', 'GLSZMZP', 'GLSZMLGLZE', 'GLSZMHGLZE', 'GLSZMSALGLE', 'GLSZMSAHGLE', 'GLSZMLALGLE', 'GLSZMLAHGLE', 'GLSZMGLV', 'GLSZMZV', 'GLRLMSRE', 'GLRLMLRE', 'GLRLMGLN', 'GLRLMRLN', 'GLRLMRP', 'GLRLMLGRE', 'GLRLMHGRE', 'GLRLMSRLGLE', 'GLRLMSRHGLE', 'GLRLMLRLGLE', 'GLRLMLRHGLE', 'GLRLMGLV', 'GLRLMRLV', 'NGTDMCoarseness', 'NGTDMContrast', 'NGTDMBusyness', 'NGTDMComplexity', 'NGTDMStrength']
        material_index = {
            'H1.nrrd':0,
            'H2.nrrd':1,
            'H3.nrrd':2,
            'H4.nrrd':3,
            'H5.nrrd':4,
            'H6.nrrd':5,
            'H7.nrrd':6,
            'H8.nrrd':7,
            'H9.nrrd':8
        }
        saveFolder = 'Radiomics'

    listProtocols = ['Radiocirugia', 'Fina', 'Mejorada', 'Normal']


    # Primero se realiza para el protocolo de referencia (mejorada)

    for measure in listMeasures:
        testNormal = np.array([])
        testFina = np.array([])
        testRadio = np.array([])

        meanRef = np.array([])
        meanNormal = np.array([])
        meanFina = np.array([])
        meanRadio = np.array([])

        for material, index in material_index.items():

            protocol_values = np.array([])
            reference_values = np.array([])

            protocol = 'Mejorada'
            listDfsRadiomics = list([])
            pathOutputCal = os.path.join(pathOutput, mode, position, protocol)
            os.chdir(pathOutputCal)
            currOutputs = glob.glob('*xlsx')
            for excel in currOutputs:
                for i in np.arange(4):
                    df = pd.read_excel(os.path.join(pathOutputCal, excel), index_col=0, sheet_name=i)
                    df = df.sort_index()
                    listDfsRadiomics.append(df)

            reference_values = np.array([])

            for df in listDfsRadiomics:
                reference_values = np.append(reference_values, (df[measure].iloc[material_index[material]]))
            
            meanRef = np.append(meanRef, np.mean(reference_values))
            
            
            # print(material)
            # print(reference_values)

            for protocol in listProtocols:
                if protocol != 'Mejorada':
                    listDfsRadiomics = list([])

                    pathOutputCal = os.path.join(pathOutput,mode, position, protocol)
                    os.chdir(pathOutputCal)
                    currOutputs = glob.glob('*xlsx')
                    for excel in currOutputs:
                        for i in np.arange(4):
                            df = pd.read_excel(os.path.join(pathOutputCal, excel), index_col=0, sheet_name=i)
                            df = df.sort_index()
                            listDfsRadiomics.append(df)

                    protocol_values = np.array([])

                    for df in listDfsRadiomics:
                        protocol_values = np.append(protocol_values,(df[measure].iloc[material_index[material]]))
                
                    # print(protocol_values)

                    if protocol == 'Normal':
                        meanNormal = np.append(meanNormal, np.mean(protocol_values))

                    if protocol == 'Fina':
                        meanFina = np.append(meanFina, np.mean(protocol_values))
                    
                    if protocol == 'Radiocirugia':
                        meanRadio = np.append(meanRadio, np.mean(protocol_values))

        alternatives = ['two-sided', 'less', 'greater']
        listNormal = np.array([])
        listFina = np.array([])
        listRadio = np.array([])

        if mode == 'CalibrationOutput':

            for alt in alternatives:

                testNormal = ranksums(meanRef, meanNormal, alternative = alt)
                listNormal = np.append(listNormal, testNormal.pvalue)
                testFina = ranksums(meanRef, meanFina, alternative = alt)
                listFina = np.append(listFina, testFina.pvalue)
                testRadio = ranksums(meanRef, meanRadio, alternative = alt)
                listRadio = np.append(listRadio, testRadio.pvalue)
        
        if mode == 'EdgesOutput':

            for alt in alternatives:

                testNormal = ranksums(np.ediff1d(meanRef), np.ediff1d(meanNormal), alternative = alt)
                listNormal = np.append(listNormal, testNormal.pvalue)
                testFina = ranksums(np.ediff1d(meanRef), np.ediff1d(meanFina), alternative = alt)
                listFina = np.append(listFina, testFina.pvalue)
                testRadio = ranksums(np.ediff1d(meanRef), np.ediff1d(meanRadio), alternative = alt)
                listRadio = np.append(listRadio, testRadio.pvalue)

        listTests = np.vstack((listNormal, listFina, listRadio))

        listProtocolsDf = ['Normal', 'Fina', 'Radiocirugia']

        df = pd.DataFrame.from_records(listTests.T, columns = [listProtocolsDf], index = alternatives)
        df.to_excel(os.path.join(pathSave, saveFolder, position, measure + '_' + position + '.xlsx'))




