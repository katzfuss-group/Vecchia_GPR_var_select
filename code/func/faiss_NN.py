import numpy as np
import faiss

def get_NN_py(x, nNbor, y=None):
    n, d = x.shape
    x = np.ascontiguousarray(x.astype(np.float32))
    if y is not None:
        y = np.ascontiguousarray(y.astype(np.float32))
    quantizer = faiss.IndexFlatL2(d)
    index = faiss.IndexIVFFlat(quantizer, d, 1024)
    index.train(x)
    index.add(x)
    index.nprobe = 256
    if y is not None:
        _, IAprx = index.search(y, int(nNbor))
    else:
        _, IAprx = index.search(x, int(nNbor))
    return IAprx
