from main import Main
import os

# Asking for path
print('Introduce path of the Data Folder to be analyzed:')
pathDDBB = input()
    
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

print('Introduce the path to the image used as a reference to perform the segmentations')
pathReference = input()

Main(pathSave, pathDDBB, position, flag, pathReference)
print('Finalizado')
