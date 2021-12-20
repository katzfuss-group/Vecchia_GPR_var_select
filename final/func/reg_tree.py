import numpy as np
from sklearn.tree import DecisionTreeRegressor


def reg_tree_wrap(XTrn, yTrn, XTst, yTst, pIn):
    # dataTrn = np.genfromtxt(dataFn + "_train.csv", delimiter=",")
    # dataTst = np.genfromtxt(dataFn + "_test.csv", delimiter=",")
    # fkMatch = re.search(r'f([0-9]+)', dataFn)
    # if fkMatch:
    #     dF = int(fkMatch.group(1))
    # else:
    #     dF = 0

    nTrn = XTrn.shape[0]
    nOOS = nTrn - int(nTrn * pIn)
    yTrn = np.array(yTrn)
    yTst = np.array(yTst)


    np.random.seed(123)
    idxOOS = np.random.choice(nTrn, nOOS, replace=False)
    XOOS = XTrn[idxOOS, :]
    XTrn = np.delete(XTrn, idxOOS, axis=0)
    yOOS = yTrn[idxOOS]
    yTrn = np.delete(yTrn, idxOOS, axis=0)
    # nTrn = nTrn - nOOS

    depths = range(3, 20)
    OOSScores = []
    idx = 0
    for i in depths:
        regTree = DecisionTreeRegressor(max_depth=i, max_features=min(100, XTrn.shape[1]), random_state=123)
        regTree.fit(XTrn, yTrn)
        yOOSPred = regTree.predict(XOOS)
        OOSScores.append(np.sqrt(np.square(yOOSPred - yOOS).mean()))

        if idx > 0 and OOSScores[idx] > OOSScores[idx - 1] * 0.99:
            break
        else:
            regTreeBst: DecisionTreeRegressor = regTree
        idx = idx + 1

    yTstPred = regTreeBst.predict(XTst)
    rmseScr = np.sqrt(np.square(yTst - yTstPred).mean() / np.var(yTst))
    idxSel = (regTreeBst.feature_importances_ > 0).nonzero()[0] + 1

    return dict({"idxSel": idxSel, "yTstPred": yTstPred, "rmseScr": rmseScr})
