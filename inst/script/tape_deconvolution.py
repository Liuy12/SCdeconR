from TAPE import Deconvolution
def tape_deconvolution(sc_data, bulk_data):
    SignatureMatrix, CellFractionPrediction = \
    Deconvolution(sc_data, bulk_data, sep='\t', scaler='mms',
                  datatype='counts',
                  mode='overall', adaptive=True, variance_threshold=0.98,
                  save_model_name=None,
                  batch_size=32, epochs=128, seed=1)
    return CellFractionPrediction
