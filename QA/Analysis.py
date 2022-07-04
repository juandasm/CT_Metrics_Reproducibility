# from MainDecay import main
from decayCompl import main
from histograms import Histograms, Calibracion, NormalityTest, WRST

# Figure selection
print('Which type of figure do you want to calculate: Intensity drop (1), Histograms (2), Calibration (3), Normal Test (4), WRST (5)')
typeFigure = int(input())

if typeFigure == 1:
    # Asking for path
    print('Introduce DDBB path:')
    pathDDBB = input()

    print('Select: [C]enter or [E]dge')
    posLetter = input()
    if posLetter == 'C' or posLetter == 'c' or posLetter == 'Center' or posLetter == 'center':
        position = 'Center'
    elif posLetter == 'E' or posLetter == 'e' or posLetter == 'Edge' or posLetter == 'edge':
        position = 'Edge'

    print('Printing intensity drop figures, introduce path to save them:')
    pathSave = input()
    main(pathSave, pathDDBB, position)

if typeFigure == 2:
    # Asking for path data
    print('From calibration [C] or from edges [E]?')
    letterMode = input()
    if letterMode == 'C' or letterMode == 'c' or letterMode == 'Calibration' or letterMode == 'calibration':
        mode = 'CalibrationOutput'
    elif letterMode == 'E' or letterMode == 'e' or letterMode == 'Edges' or letterMode == 'edges':
        position = 'EdgesOutput'

    print('Introduce path to output folder:')
    pathData = input()

    print('Select: [C]enter or [E]dge')
    posLetter = input()
    if posLetter == 'C' or posLetter == 'c' or posLetter == 'Center' or posLetter == 'center':
        position = 'Center'
    elif posLetter == 'E' or posLetter == 'e' or posLetter == 'Edge' or posLetter == 'edge':
        position = 'Edge'
    print('Select which measure are you interested in: [1] Mean, [2] CoV, [3] SD, [4] SUVpeak, [9] DR_mean, [10] TC_mean')
    flagMeasure = int(input())
    print('Printing intensity drop figures, introduce path to save them:')
    pathSave = input()
    Histograms(pathData, pathSave, position, mode, flagMeasure)

if typeFigure == 3:
    # Asking for path data
    print('Introduce path to calibration output folder:')
    pathData = input()
    print('Select: [C]enter or [E]dge')
    posLetter = input()
    if posLetter == 'C' or posLetter == 'c' or posLetter == 'Center' or posLetter == 'center':
        position = 'Center'
    elif posLetter == 'E' or posLetter == 'e' or posLetter == 'Edge' or posLetter == 'edge':
        position = 'Edge'
    print('Select: [I]n or [O]ut')
    locLetter = input()
    if locLetter == 'I' or locLetter == 'i' or locLetter == 'In' or locLetter == 'in':
        loc = 'In'
    elif locLetter == 'O' or locLetter == 'o' or locLetter == 'Out' or locLetter == 'out':
        loc = 'Out'
    print('Printing calibration graphs, introduce path to save them:')
    pathSave = input()
    Calibracion(pathData, pathSave, position, loc)

if typeFigure == 4:
    # Asking for path data
    print('Introduce path to output folder:')
    pathData = input()
    print('Which metrics do you want to check if they fit a normal distribution: (1) Calibration, (2) Edge, (3) Radiomics?')
    flagMode = int(input())
    if flagMode == 1:
        mode = 'CalibrationOutput'
    if flagMode == 2:
        mode = 'EdgesOutput'
    if flagMode == 3:
        mode = 'RadiomicsOutput'
    print('Select: [C]enter or [E]dge')
    posLetter = input()
    if posLetter == 'C' or posLetter == 'c' or posLetter == 'Center' or posLetter == 'center':
        position = 'Center'
    elif posLetter == 'E' or posLetter == 'e' or posLetter == 'Edge' or posLetter == 'edge':
        position = 'Edge'
    print('Introduce path to save results:')
    pathSave = input()
    NormalityTest(pathData, pathSave, position, mode)

if typeFigure == 5:
    # Asking for path data
    print('Introduce path to output folder:')
    pathData = input()
    print('Which metrics do you want to check if they fit a normal distribution: (1) Calibration, (2) Edge, (3) Radiomics?')
    flagMode = int(input())
    if flagMode == 1:
        mode = 'CalibrationOutput'
    if flagMode == 2:
        mode = 'EdgesOutput'
    if flagMode == 3:
        mode = 'RadiomicsOutput'
    print('Select: [C]enter or [E]dge')
    posLetter = input()
    if posLetter == 'C' or posLetter == 'c' or posLetter == 'Center' or posLetter == 'center':
        position = 'Center'
    elif posLetter == 'E' or posLetter == 'e' or posLetter == 'Edge' or posLetter == 'edge':
        position = 'Edge'
    print('Introduce path to save results:')
    pathSave = input()
    WRST(pathData, pathSave, mode, position)
