import tellurium as te

import SyntheticGRN as syn
modelStr = syn.GetModel(3, 'sauro', 'c:\\tmp\\sauro\\', [1, 0,0])

if modelStr == '':
    print 'No viable model found'
else:
    r = te.loada(modelStr)

r.draw()