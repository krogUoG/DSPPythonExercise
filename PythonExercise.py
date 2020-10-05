import numpy as np
from matplotlib import pyplot as plt

""" Functions """

""" Generate na array of n samples between between two values
    this method initalizes list first then appends all elements 
    using the for loop and then converts it into numpy array and returns """
def GenerateNumpyArrayFromList(startValue, endValue, nSamples):
    # Initialize empty list
    listToArray = []
    # Calculate Step size between samples ( range / num of samples)
    stepSize = (endValue - startValue) / nSamples
    # Populate an Array with values
    print( "Populting an array of " + str(nSamples) + " samples between " + str(startValue) + " and " + str(endValue) + " using list and for loop")
    for i in range(nSamples):
        listToArray.append( (i*stepSize) )
        
    # Form and return array
    return np.asarray(listToArray)

""" Generate na array of n samples between between two values """
def GenerateNumpyArray(startValue, endValue, nSamples):
    print( "Generating an array of " + str(nSamples) + " samples between " + str(startValue) + " and " + str(endValue) )
    return np.linspace( startValue, endValue, num=nSamples, endpoint=True, retstep=False, dtype=None, axis=0 )

""" Generate a sine wave using an array of samples, sample by sample """
def GenerateSineWaveFromSamples(SampleArray):
    print("Generating a sine signal sample by smaple")
    AmplitudeArray = SampleArray
    
    for i in len(AmplitudeArray):
        AmplitudeArray[i] = np.sin(AmplitudeArray[i])
    return AmplitudeArray

""" Generate a sine wave using an array of samples """
def GenerateSineWaveFromArray( SampleArray ):
    print("Generating a sine signal from set of samples")
    
    return np.sin( SampleArray )

def GenerateNSinWaveCycles(nCycles, nSamples):
    # Calculate time steps 
    timeSteps = np.linspace( 0, (2*np.pi), num=nSamples, endpoint=True, retstep=False, dtype=None, axis=0 )
    sineWaveAmpitude = np.sin( timeSteps * nCycles )
    
    return sineWaveAmpitude
    

""" Perform Frequency modulation on set of samples """
def ModulationFM(inputSampleArray, signalFreq, carrierFreq, nSamples ):
    print("Frequency Modulating the signal")

    # Calculate delta frequency
    deltaFreq = np.abs(carrierFreq - signalFreq)
    # Calculate angular frequencies
    signalAngularFreq = signalFreq * 2 * np.pi
    carierAngularFreq = carrierFreq * 2 * np.pi
    
    # Generate time variable array 
    timeSteps = np.linspace( 0, (2*np.pi), num=8000, endpoint=True, retstep=False, dtype=None, axis=0 )
    
    modulatedOutput = np.cos( carierAngularFreq*timeSteps + ( (deltaFreq/signalFreq) * np.sin( signalAngularFreq * timeSteps) ) )
    return modulatedOutput
    
""" Plot x,y graph using two sets of values """
def Plot2dGraph(title, xcords, ycords, xlabel, ylabel):
    plt.title(title) 
    plt.xlabel(xlabel) 
    plt.ylabel(ylabel) 
    plt.plot(xcords,ycords) 
    plt.show()
    
""" Perform a half-wave rectification on a set of samples """
def HalfWaveRectification(SignalAmplitude):
    outputAmplitude = SignalAmplitude
    print("Half-wave rectifing the signal")
    # For each sample
    for i in range(len(outputAmplitude)):
        # If the amplitude is negative
        if outputAmplitude[i] < 0.0:
            # Set amplitude to 0
            outputAmplitude[i] = 0.0
            
    return outputAmplitude

""" Main function """

# Generate number arrays
startValue = 0
endValue = (2 * np.pi)
nSamples = 8000
nCycles =  200
signalFreq = 10
carrierFreq = 20

#sampleNumpyArray = GenerateNumpyArray( startValue, endValue, nSamples)
sampleNumpyArray = GenerateNumpyArrayFromList(startValue, endValue, nSamples)
amplitudeNumpyArray = GenerateSineWaveFromArray(sampleNumpyArray)
#Plot2dGraph( "Sin (x) betwen 0 and 2pi ", sampleNumpyArray, amplitudeNumpyArray, "Angle [rad]", "Sin(x)")

SineAmplitudeArray = GenerateNSinWaveCycles(nCycles,nSamples)
#Plot2dGraph( str(nCycles) + "Hz Sin (x) betwen 0 and 2pi ", amplitudeNumpyArray, SineAmplitudeArray, "Angle [rad]", "Sin(wx)")

FMModulatedOutput = ModulationFM(sampleNumpyArray, signalFreq, carrierFreq, nSamples )
Plot2dGraph( "FM modulated Sin Wave, Signal Freq: " + str(signalFreq)  + "Hz Carrier Freq: " + str(carrierFreq) + "Hz, betwen 0 and 2pi", sampleNumpyArray, FMModulatedOutput, "Angle [rad]", "Amplitude")

rectifiedOutput = HalfWaveRectification(FMModulatedOutput)
Plot2dGraph( "Rectified FM modulated Sin Wave, Signal Freq: " + str(signalFreq)  + "Hz Carrier Freq: " + str(carrierFreq) + "Hz, betwen 0 and 2pi", sampleNumpyArray, rectifiedOutput, "Angle [rad]", "Amplitude")

#numpyArrayFromList = GenerateNumpyArrayFromList( startValue, endValue, nSamples)