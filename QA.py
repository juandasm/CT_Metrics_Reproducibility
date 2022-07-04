from main import Main
import os

# Asking for path
print('Introduce DDBB path:')
pathDDBB = input()

# Cd to path
os.chdir(pathDDBB)
dataFolders = os.listdir(pathDDBB)

# Readind txt file folders
analysedFolders = list([])
with open(".info/dataHistory.txt", "r") as file:
    analysedFolders=file.read().splitlines()

# Readind txt file original folder
with open(".info/dataOriginal.txt", "r") as file:
    originalFolder=file.read().splitlines()
    originalFolder = originalFolder[0]

# Reading existing folders
listFolders = list([])
for folder in dataFolders:
    if os.path.isdir(os.path.join(pathDDBB,folder)) and folder[0] is 'D':
        listFolders.append(folder)
listFolders.sort()

# Seeing if there are new folders
if listFolders == analysedFolders:
    print('DDBB has not sufferd any change since last time, neither resampling nor registering necessary')
else:
    newFolder = list([])
    for folder in listFolders:
        if folder not in analysedFolders:
            newFolder.append(folder)
    
    print('New data folder found: ' + str(newFolder))
    print('Do you want to register to ' + originalFolder + ', or use a dedicated segmentation? [O]ld, [N]ew')
    answer_seg = input()
    if answer_seg == 'O' or answer_seg == 'Old' or answer_seg == 'o':
        print('Registering using old segmentations')
        print('\n')
        print('Center [1] or Edge [2]?')
        answerPos = input()
        if answerPos == 1:
            position = 'Center'
        elif answerPos == 2:
            position = 'Edge'
        else:
            print('error')
            
        print('Which segmentation do you want to use?')
        counter = 1
        for folder in analysedFolders:
            print(folder + '[' + str(counter) + '], ')
            counter = counter + 1
        answer = input()
        print('\n')
        print('Using folder ' + analysedFolders[int(answer)-1])
    elif answer_seg == 'N' or answer_seg == 'New' or answer_seg == 'n':
        print('Using new segmentations')
    else:
        print('Learn to write')
    
print('Select which calculations do you want to do: Calibration (1), Edges (2), Radiomics (3), All (4)')
flag = int(input())

print('Which position are you interested in? [C]enter or [E]dge')
answerPos = input()

if answerPos == 'Center' or answerPos == 'center' or answerPos == 'C' or answerPos == 'c':
    position = 'Center'
elif answerPos == 'Edge' or answerPos == 'edge' or  answerPos == 'E' or answerPos == 'e':
    position = 'Edge'

print('Set path to save results:')
pathSave = input()

Main(pathSave, pathDDBB, position, flag)
