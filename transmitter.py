import matlab.engine
import time
import os

if __name__ == '__main__':

    # Configure experiment
    numTxAntennas = 4  # Select between 1, 2, 3 and 4
    gain = -10.0  # in dB
    fileName = 'data/channelEstimation.bin'  # Store channel estimation measures

    start1 = time.time()
    myMatlab = matlab.engine.start_matlab()
    end1 = time.time()
    print "matlab startup time ", end1 - start1, " seconds"

    # Add current directory to Matlab path
    cwd = os.getcwd()
    myMatlab.addpath(cwd, nargout=0)

    start2 = time.time()
    [transmitter, bits, symbols, payload, trainingSig, chEst_prelim] = myMatlab.f_loadConfig(numTxAntennas, gain, fileName, nargout=6)
    end2 = time.time()

    # print"t(Load config) = ", end2 - start2, " seconds"
    # print"Type of Bits: ", type(bits)
    # print"Length of Bits: ", len(bits)
    # print"Length of Bits: ", len(bits[0])
    # for i in range(0, 10):
    #     print"sample bits[0][", i, "] = ", str(bits[0][i][0])
    #     print"type bits[0][", i, "] = ", type(bits[0][i][0])
    #     print"sample trainingSig[", i, "][0] = ", str(trainingSig[i][0])
    #     print"type trainingSig[", i, "][0] = ", type(trainingSig[i][0])

    start2 = time.time()
    [txSig, chEst_new] = myMatlab.f_BFPayload(numTxAntennas, payload, trainingSig, fileName, chEst_prelim, nargout=2)
    end2 = time.time()

    # print"t(Compute weights) = ", end2 - start2, " seconds"
    # print"Type of channelEst_new: ", type(chEst_new)
    # print"Length of channelEst_new: ", len(chEst_new)
    # print"Length of channelEst_new: ", len(chEst_new[0])
    # for i in range(0, 4):
    #     print"type channelEst_prelim[0][", i, "] = ", chEst_prelim[0][i].real
    #     print"type channelEst_new[0][", i, "] = ", chEst_new[0][i].real
    #     print"type channelEst_prelim[0][", i, "] = ", chEst_prelim[0][i].imag
    #     print"type channelEst_new[0][", i, "] = ", chEst_new[0][i].imag

    start2 = time.time()s
    eng.f_transmit(transmitter,nargout=0)
    end2 = time.time()
    print"t(Transmit signal) = ", end2 - start2, " seconds"