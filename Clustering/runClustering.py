from build.clustering import clustering
import numpy as np
from scipy.spatial import cKDTree
import os
from pathlib import Path
import sys
import inspect
import uproot
import uproot.models.TTree as TTreeType
from typing import Union, Tuple, List


currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
builddir = currentdir+r"/build"

def loadSugarFile(sugarPosFilename: str) -> Tuple[cKDTree, cKDTree]:
    """Function to read in sugar positions file

    Args:
        sugarPosFilename (str): Path to sugar positions file

    Returns:
        Tuple[cKDTree, cKDTree]: KDTrees containing the positions of the sugars on stand 0 and strand 1
    """
    # Read sugar file
    with open(sugarPosFilename, "rb") as f:
        data = np.fromfile(f,np.float32)

    data = data.reshape(int(len(data)/12),12)
    
    sugar0 = data[:,0:3]
    sugar1 = data[:,3:6]

    sugar0 = np.array(sugar0)
    sugar1 = np.array(sugar1)

    T0 = cKDTree(sugar0)
    T1 = cKDTree(sugar1)

    return T0, T1

def doseCalculation(input_tree: TTreeType, pathLength: dict, dosePerEvent: dict, meanKEperEvent: dict, eventType: str, chromatinVolume: float, copyCheck: bool = False, copy: int = 0, part1_particleSource: Union[bool, list] = False, part1_particleTime: Union[float, bool] = False) -> Tuple[Union[str,dict], dict, dict]:
    """ Calculates the dose to the whole box from the simulation output. Selecting box copy number and particle as required. If the required fields are available, also calculates the mean KE and path length per event, and mean KE of the primary particle.

    Args:
        input_tree (TTreeType): TTree from simulation containing energy deposits in box per event
        pathLength (dict): dictionary to save total path length per event
        dosePerEvent (dict): dictionary to save total dose per event
        meanKEperEvent (dict): dictionary to save mean KE per event
        eventType (str): type of simulation, either standalone or PS
        chromatinVolume (float): volume of simulation box
        copyCheck (bool, optional): True if simulation results to be returned only from one box number
        copy (int, optional): Box number to use
        part1_particleSource (bool or list): list of primary particles to use, default: false
        part1_particleTime (bool or float): maximum time to include particles, default: false


    Returns:
        Tuple[Union[str,dict], dict, dict, dict, dict]: Mean primary energy (if available, event mapping from primary event to PS event (empty if standalone), dose per event), mean KE per event, path length per event
    """

    eventMapping = {}

    if eventType == "standalone":
        event = input_tree["EventNum"].array(library = "np")
    elif eventType == "photon":
        event = input_tree["part1_EventNum"].array(library = "np")
    elif eventType == "decay":
        event = input_tree["EventNum"].array(library = "np")
    
    if copyCheck:
        part1_CopyNum = input_tree["part1_CopyNum"].array(library = "np")
        mask = part1_CopyNum==copy
    else:
        mask = np.array([True]*len(event))

    if part1_particleSource:
        PrimaryParticle = input_tree["part1_particleSource"].array(library = "np")

        tempMask = np.zeros_like(mask)
        for i in part1_particleSource:
            tempMask+=(PrimaryParticle==i) # returns a truth array for each particle in requested particles. Temp mask changed to true if particle in requested list.
        mask *=tempMask

    if part1_particleTime:
        particleTime = input_tree["time"].array(library = "np")
        tempMask = (particleTime<= part1_particleTime ) 
        mask *=tempMask

    if ("PathLengthChromatin" in input_tree) and ("PrimaryKEEntrance" in input_tree):
        data = np.vstack([event, input_tree["Edep_J"].array(library = "np"), input_tree["PathLengthChromatin"].array(library = "np"), input_tree["PrimaryKEEntrance"].array(library = "np"), input_tree["PrimaryKEExit"].array(library = "np")]).T
    elif ("PathLengthChromatin" in input_tree):
        data = np.vstack([event, input_tree["Edep_J"].array(library = "np"), input_tree["PathLengthChromatin"].array(library = "np")]).T
    else:
        data = np.vstack([event, input_tree["Edep_J"].array(library = "np")]).T

    data = data[mask]
    volumeFactor = 1/ (1000 * chromatinVolume)

    # remove events with no dose
    data = data[data[:,1]>0]
    # apply mask to events
    event = event[mask]

    if eventType == "standalone":
        dosePerEvent = {e: volumeFactor * np.sum((data[data[:,0]==e])[:,1]) for e in event if (np.sum((data[data[:,0]==e])[:,1])>0)} #only include if dose>0
        
        pathLength = {e: (data[data[:,0]==e])[0,2] for e in event if ("PathLengthChromatin" in input_tree) and (np.sum((data[data[:,0]==e])[:,1])>0) }  #only include if dose>0

        meanKEperEvent = {e: np.sum((data[data[:,0]==e])[0,3:5])/2 for e in event if ("PrimaryKEEntrance" in input_tree) and (np.sum((data[data[:,0]==e])[:,1])>0)} #only include if dose>0



    else:
        dosePerEventNew = {e: volumeFactor * np.sum((data[data[:,0]==e])[:,1]) for e in event if np.sum((data[data[:,0]==e])[:,1])>0}
        dosePerEvent = {e: np.sum([dosePerEventNew[e], dosePerEvent[e]]) if (e in dosePerEvent) and (e in dosePerEventNew) else dosePerEventNew[e] if (e in dosePerEventNew) else dosePerEvent[e] for e in set(list(dosePerEvent.keys())+list(dosePerEventNew.keys()))}
    

    if eventType == "photon":
        eventMapping = {a:b for a,b in zip(input_tree["EventNum"].array(library = "np"), input_tree["part1_EventNum"].array(library = "np"))} 


    if meanKEperEvent and list(meanKEperEvent.values())[0]!="N/A":
        energy = sum(meanKEperEvent.values())/len(meanKEperEvent.values())
    else:
        energy = "N/A"

    return energy, eventMapping, dosePerEvent, meanKEperEvent, pathLength

def getNumIntersecting(fInfo: TTreeType) -> float:
    """Returns the number of particles crossing the target box

    Args:
        fInfo (TTreeType): TTree containing information about the simulaiton

    Returns:
        float: number of particles crossing the target box.
    """
    return fInfo["NumIntersecting"].array()[0]

def AccumulateEdep(input_tree: TTreeType, cumulatedEnergyDep: dict, T0: cKDTree, T1: cKDTree, eventMapping: Union[bool,dict] = False, copyCheck: bool = False, copy: int = 0, part1_particleSource: Union[bool, list]  = False, part1_particleTime: Union[float, bool] = False)->dict:
    """Accumulate energy deposits from all steps if near to sugar molecule

    Args:
        input_tree (TTreeType): TTree containing energy deposits per step
        cumulatedEnergyDep (dict): dictionary to accumulate energy deposits in sugar volumes per event
        T0 (cKDTree): positions of strand 0
        T1 (cKDTree): positions of strand 1
        eventMapping (dict, optional): if required mapping between primary simulation event number and PS simulation event number. Defaults to False.
        copyCheck (bool, optional): Whether to only accumulate energy from a specific box number. Defaults to False.
        copy (int, optional): box number to use. Defaults to 0.
        part1_particleSource (bool or list): list of primary particles to use, default: false
        part1_particleTime (bool or float): maximum time to include particles, default: false

    Returns:
        dict: dictionary containing all energy deposits per sugar per event
    """

    if eventMapping:
        event = np.array([eventMapping[a] for a in input_tree["EventNum"].array(library = "np")])
    else:
        event = input_tree["EventNum"].array(library = "np")

    if copyCheck:
        part1_CopyNum = input_tree["part1_CopyNum"].array(library = "np")
        mask = part1_CopyNum==copy
    else:
        mask = np.array([True]*len(event))

    if part1_particleSource:
        PrimaryParticle = input_tree["part1_particleSource"].array(library = "np")

        tempMask = np.zeros_like(mask)
        for i in part1_particleSource:
            tempMask+=(PrimaryParticle==i)
        mask *=tempMask

    if part1_particleTime:
        particleTime = input_tree["time"].array(library = "np")
        tempMask = (particleTime<= part1_particleTime ) 
        mask *=tempMask

    data = np.vstack([event, input_tree["x"].array(library = "np"),  input_tree["y"].array(library = "np"),  input_tree["z"].array(library = "np"), input_tree["eDep_eV"].array(library = "np")]).T

    data = data[mask]

    cumulatedEnergyDep = {}

    results = checkPoints(data[:,1:4], T0,T1)

    for i in range(len(data)):
        if (results[i][0] != -1):
            key = (int(data[i,0]),results[i][0],results[i][1])
            if key in cumulatedEnergyDep:
                cumulatedEnergyDep[key] += data[i,4]
            else:
                cumulatedEnergyDep[key] = data[i,4]
    
    return cumulatedEnergyDep

def Direct(cumulatedEnergyDep: dict, fEMinDamage: float, fEMaxDamage: float)-> Tuple[list, list, list]:
    """Convert accumulated energy deposition per sugar into strand breaks using the linear probability model. Exports lists of event number, sugar copy number and strand number for input into clustering algorithm.

    Args:
        cumulatedEnergyDep (dict): dictionary containing all energy deposits per sugar per event
        fEMinDamage (float): minimimum energy deposit for a strand break
        fEMaxDamage (float): maximum energy above which a strand break always occurs

    Returns:
        eventsListDirect, copyListDirect, strandListDirect Tuple(list, list, list): lists of event number, sugar copy number and strand number for input into clustering algorithm
    """

    isBreak = [IsEdepSufficient(cumulatedEnergyDep[key],fEMinDamage, fEMaxDamage) for key in cumulatedEnergyDep]

    eventsListDirect = [a[0] for a,b in zip(cumulatedEnergyDep.keys(), isBreak) if b]
    strandListDirect = [a[1] for a,b in zip(cumulatedEnergyDep.keys(), isBreak) if b]
    copyListDirect = [a[2] for a,b in zip(cumulatedEnergyDep.keys(), isBreak) if b]

    return eventsListDirect, copyListDirect, strandListDirect

def Indirect(input_tree: TTreeType, probIndirect: float, T0: cKDTree, T1: cKDTree, eventMapping: Union[dict, bool] = False, copyCheck: bool = False, copy:int = 0,part1_particleSource: Union[bool, list] = False, part1_particleTime: Union[float, bool] = False)->Tuple[list, list, list]:
    """ Calculate the number of indirect stand breaks (OH + sugar->damaged sugar) and returns lists of the events, sugar copy number and strand number for the clustering algorithm

    Args:
        input_tree (TTreeType): TTree from simulation containing water radical reactions with DNA per event
        probIndirect (float): probability that a OH + sugar reaction results in a strand break
        T0 (cKDTree):  positions of strand 0
        T1 (cKDTree):  positions of strand 1
        eventMapping (dict, optional): if required mapping between primary simulation event number and PS simulation event number. Defaults to False.
        copyCheck (bool, optional): True if simulation results to be returned only from one box number
        copy (int, optional): Box number to use
        part1_particleSource (Union[bool, list], optional): Whether to accumulate from a specific particle, numbering of particles 1 indexed. Defaults to False.
        part1_particleTime (bool or float): maximum time to include particles, default: false

    Returns:
        eventsListIndirect, copyListIndirect, strandListIndirect Tuple[list, list, list]: lists of event number, sugar copy number and strand number for input into clustering algorithm
    """
    eventsListIndirect=[]
    copyListIndirect=[]
    strandListIndirect=[]

    if eventMapping:
        event = np.array([eventMapping[a] for a in input_tree["EventNum"].array(library = "np")])
    else:
        event = input_tree["EventNum"].array(library = "np")

    if len(event)==0: #RBE can be run with indirect off
        return eventsListIndirect, copyListIndirect, strandListIndirect 

    if copyCheck:
        part1_CopyNum = input_tree["part1_CopyNum"].array(library = "np")
        mask = part1_CopyNum==copy
    else:
        mask = np.array([True]*len(event))
        
    if part1_particleSource:
        PrimaryParticle = input_tree["part1_particleSource"].array(library = "np")

        tempMask = np.zeros_like(mask)
        for i in part1_particleSource:
            tempMask+=(PrimaryParticle==i)
        mask *=tempMask

    if part1_particleTime:
        particleTime = input_tree["time"].array(library = "np")
        tempMask = (particleTime<= part1_particleTime ) 
        mask *=tempMask

    data = np.vstack([event, input_tree["DNAmolecule"].array(library = "np"),  input_tree["radical"].array(library = "np"), input_tree["x"].array(library = "np"),  input_tree["y"].array(library = "np"),  input_tree["z"].array(library = "np")]).T

    data = data[mask]

    # extract only deoxyribose and OH reactions
    data = data[data[:,1]=="Deoxyribose^0"]
    data = data[data[:,2]=="OH^0"]

    results = getIndex(data[:,3:6], T0,T1)

    for i in range(len(data)):
        if np.random.rand() <= probIndirect:
            eventsListIndirect.append(data[i,0])
            copyListIndirect.append(results[i][0])
            strandListIndirect.append(results[i][1])

    return eventsListIndirect, copyListIndirect, strandListIndirect

def IsEdepSufficient(pEdep:float, fEMinDamage:float, fEMaxDamage:float) -> bool:
    """Apply the linear probability model to determine if a direct strand break occurs.

    Args:
        pEdep (float): energy deposit in a sugar molecule
        fEMinDamage (float): minimimum energy deposit for a strand break
        fEMaxDamage (float): maximum energy above which a strand break always occurs

    Returns:
        bool: whether strand break occurs
    """
    if (pEdep<fEMinDamage):
        return False
    if(pEdep>fEMaxDamage):
        return True
    else:
        proba = (pEdep - fEMinDamage)/(fEMaxDamage-fEMinDamage)
        return (proba>np.random.rand())

def checkPoints(points: List[List], T0: cKDTree, T1: cKDTree) -> List[List]:
    """ Check points are within Rdirect (0.35nm) of a sugar molecule centre

    Args:
        points (List[List]): x,y,z positions of energy depositions
        T0 (cKDTree): sugar molecules on strand 0
        T1 (cKDTree): sugar molecules on strand 1

    Returns:
        result (list): [strand, sugar index] if within Rdirect [-1,-1] otherwise
    """
    result = [-1,-1]
    Rdirect = 0.35 #nm
    d0, idx0 = T0.query(points, k=1)
    d1, idx1 = T1.query(points, k=1)

    result = [[0, c] if (a<b and a< Rdirect) else [1,d] if (b < a and b<Rdirect) else [-1,-1] for a,b,c,d in zip(d0,d1,idx0,idx1) ]

    return result

def getIndex(points: list, T0: cKDTree, T1: cKDTree) -> List[List]:
    """Get sugar copy number and strand number from x,y,z positions of reactions

    Args:
        points (list): position of the reaction
        T0 (cKDTree): positions of sugar molecules on strand 0
        T1 (cKDTree): positions of sugar molecules on strand 1

    Raises:
        ValueError: Error if no sugar molecule nearby

    Returns:
        List[List]: [sugar copy number, sugar strand]
    """

    d0, idx0 = T0.query(points, k=1)
    d1, idx1 = T1.query(points, k=1)

    result = [[c, 0] if (a<b) else [d,1] for a,b,c,d in zip(d0,d1,idx0,idx1)]

    d = [min(a,b) for a,b in zip(d0, d1)]

    if any(i > 1e-5 for i in d):
        raise ValueError("d > 1e-10")

    return result

def readClusteringGitHash(builddir: str)->list:
    """Read text file containing git hash when clustering algorithm last compiled

    Args:
        builddir (str): path to text file
    Returns:
        list[str]: git hash
    """
    import subprocess
    result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=currentdir,
                            capture_output=True, text=True)
    return result.stdout.strip() or "unversioned"

def runClustering(filename: str, outputFilename: str, fEMinDamage: float, fEMaxDamage: float, probIndirect: float, sugarPosFilename: str, simulationType: str="standalone", filenamePhoton: Union[str,bool]=False, continuous: bool = True, part1_CopyNum: Union[int,bool] = False, part1_particleSource: Union[List[str],bool]=False, part1_particleTime: Union[float,bool]=False):
    """ Run all stages to interpret the root file from the simulation and output the number of strand breaks.

    Args:
        filename (str): root file to analyse
        outputFilename (str): filename to save results
        fEMinDamage (float): minimimum energy deposit for a strand break
        fEMaxDamage (float): maximum energy above which a strand break always occurs
        probIndirect (float): probability that OH+sugar reaction leads to strand damage
        sugarPosFilename (str): text file containing the sugar molecule positions, must match the one used for simulation
        simulationType (str, optional): Type of simulation standalone, decay, photon. Defaults to "standalone".
        filenamePhoton (Union[str,bool], optional): If photon, path to root file from photon simulation. Defaults to False.
        continuous (bool, optional): Continuous DNA structure from fractal DNA - True, original DNA structure with discontinuities - False. Defaults to True.
        part1_CopyNum (Union[int,bool], optional): If decay, which copy number to use for clustering. Defaults to False.
        part1_particleSource (Union[str,bool], optional): If decay, list of particles to use for clustering, if False all. Defaults to False.
        part1_particleTime (Union[float, bool], optional): If decay, maximum time to include in clustering, if False all. Defaults to False.
        

    Raises:
        ValueError: Incompatible simulationType given
    """
    if simulationType == "ps":
        simulationType = "decay"
    particleMap = {
      "Ra224": 1,
      "Rn220": 2,
      "Po216": 3,
      "Pb212": 4,
      "Bi212": 5,
      "Tl208": 6,
      "Po212": 7,
      "Pb208": 8,
      "alphaRa224": 9,
      "alphaRn220": 10,
      "alphaPo216": 11,
      "alphaBi212": 12,
      "alphaPo212": 13,
      "e-Rn220": 14,
      "e-Po216": 15,
      "e-Pb212": 16,
      "e-Bi212": 17,
      "e-Tl208": 18,
      "e-Po212": 19,
      "e-Pb208": 20,
      "gammaRn220": 21,
      "gammaPo216": 22,
      "gammaPb212": 23,
      "gammaBi212": 24,
      "gammaTl208": 25,
      "gammaPo212": 26,
      "gammaPb208": 27,
      "e+": 28,

      "At211": 29,
      "Po211": 30,
      "Bi207": 31,
      "Pb207": 32,
      "alphaAt211": 33,
      "alphaPo211": 34,
      "e-At211": 35,
      "e-Po211": 36,
      "e-Bi207": 37,
      "e-Pb207": 38,

      "gammaAt211": 39,
      "gammaPo211": 40,
      "gammaBi207": 41,
      "gammaPb207": 42,


      "Ra226": 43,
      "Rn222": 44,
      "Po218": 45,
      "Pb214": 46,
      "Bi214": 47,
      "Po214": 48,
      "Pb210": 49,
      "Bi210": 50,
      "Po210": 51,
      "Pb206": 52,
      "Tl210": 81,
      "Tl206": 82,
      "Hg206": 83,
      "At218": 90,
      "Rn218": 91,




      "alphaRa226": 53,
      "alphaRn222": 54,
      "alphaPo218": 55,
      "alphaPo214": 56,
      "alphaPo210": 57,

      "e-Ra226": 58,
      "e-Rn222": 59,
      "e-Po218": 60,
      "e-Pb214": 61,
      "e-Bi214": 62,
      "e-Po214": 63,
      "e-Pb210": 64,
      "e-Bi210": 65,
      "e-Po210": 66,
      "e-Pb206": 67,
      "e-Tl210": 84,
      "e-Tl206": 85,
      "e-Hg206": 86,
      "e-At218": 92,
      "e-Rn218": 93,


      "gammaRa226": 68,
      "gammaRn222": 69,
      "gammaPo218": 70,
      "gammaPb214": 71,
      "gammaBi214": 72,
      "gammaPo214": 73,
      "gammaPb210": 74,
      "gammaBi210": 75,
      "gammaPo210": 76,
      "gammaPb206": 77,
      "gammaTl210": 87,
      "gammaTl206": 88,
      "gammaHg206": 89,
      "gammaAt218": 94,
      "gammaRn218": 95,


      
      "Co60": 78,            
      "gamma": 79,         
      "e-Co60": 80,

      "Ra223": 181,
      "Rn219": 182,
      "Po215": 183,
      "Pb211": 184,
      "Bi211": 185,
      "Tl207": 186,
    #   "Pb207_Ra223": 87,  # لتجنب التضارب مع Pb207 من At211

      "alphaRa223": 187,
      "alphaRn219": 188,
      "alphaPo215": 189,
      "alphaBi211": 190,

      "e-Ra223": 191,
      "e-Rn219": 192,
      "e-Po215": 193,
      "e-Pb211": 194,
      "e-Bi211": 195,
      "e-Tl207": 196,

      "gammaRa223": 197,
      "gammaRn219": 198,
      "gammaPo215": 199,
      "gammaPb211": 100,
      "gammaBi211": 101,
      "gammaTl207": 102,


      "Ra225": 103,
      "Rn221": 104,
      "Po217": 105,
      "Pb213": 106,
      "Ac225": 107,
      "Fr221": 108,
      "Ra221": 109,
      "Rn217": 110,
      "Po213": 111,
      "Pb209": 112,
      "Bi209": 113,
      "Tl205": 114,
      "At217": 115,
      "Bi213": 116,
      "Tl209": 117,



      "alphaRa225": 118,
      "alphaRn221": 119,
      "alphaPo217": 120,
      "alphaAc225": 121,
      "alphaFr221": 124,
      



      "alphaRa221": 125,
      "alphaRn217": 126,
      "alphaPo213": 127,
      "alphaBi209": 128,
      "alphaAt217": 129,
      "alphaBi213": 130,


      "e-Fr221": 131,
      "e-Tl209": 132,
      "e-Ra225": 133,
      "e-Pb209": 134,
      "e-Pb213": 135,
      "e-Bi213": 136,
      "e-Po217": 137,
      "e-Rn221": 138,






      }
    
    if part1_particleSource:
        part1_particleSource = [particleMap[a] for a in part1_particleSource]

    T0, T1 = loadSugarFile(sugarPosFilename)

    # Load root data to analyse
    file = uproot.open(filename)
    
    fDirectory = file['ntuple']

    tInfo = fDirectory["Info"]
    tEdep = fDirectory["EventEdep"]
    input_tree = fDirectory["Direct"]
    input_tree_indirect = fDirectory["Indirect"]

    # Simulation Parameters
    chromatinVolume=tInfo["ChromatinVolume_m3"].array()[0] # in m3
    numBP = tInfo["NumBasepairs"].array()[0]
    gitHash = tInfo["GitHash"].array()[0]

    assert T0.n==numBP, "Sugar file does not match simulation sugar file"

    pathLength={}
    dosePerEvent={}
    meanKEperEvent={}
    cumulatedEnergyDep = {}

    if simulationType == "standalone":
        energy, _, dosePerEvent, meanKEperEvent, pathLength = doseCalculation(tEdep, pathLength, dosePerEvent,meanKEperEvent, "standalone", chromatinVolume) 
        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1)
        eventsListDirect, copyListDirect, strandListDirect = Direct(cumulatedEnergyDep,fEMinDamage, fEMaxDamage)
        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree_indirect, probIndirect, T0, T1)

        LET = np.mean(tInfo["MeanLET"])

        NumIntersecting = len(pathLength.keys())

    elif simulationType == "decay":
        energy, _, dosePerEvent, meanKEperEvent, pathLength = doseCalculation(tEdep, pathLength, dosePerEvent,meanKEperEvent, "decay", chromatinVolume, copyCheck = True, copy = part1_CopyNum, part1_particleSource = part1_particleSource, part1_particleTime = part1_particleTime) 
        

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, copyCheck = True, copy = part1_CopyNum, part1_particleSource = part1_particleSource, part1_particleTime = part1_particleTime)

        eventsListDirect, copyListDirect, strandListDirect = Direct(cumulatedEnergyDep,fEMinDamage, fEMaxDamage)

        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree_indirect, probIndirect, T0, T1, copyCheck = True, copy = part1_CopyNum, part1_particleSource = part1_particleSource, part1_particleTime = part1_particleTime)
        LET = "N/A"

        NumIntersecting = len(list(set(eventsListDirect+eventsListIndirect)))
        # Path length and KE are not recorded for decay simulations
        pathLength = {e: "N/A" for e in dosePerEvent.keys()}
        meanKEperEvent = {e: "N/A" for e in dosePerEvent.keys()}

    elif simulationType == "photon":
        fFilePhoton = uproot.open(filenamePhoton) 
        fDirectoryPhoton = fFilePhoton["ntuple"]
        tEdepPhoton = fDirectoryPhoton["EventEdep"]
        input_tree_photon = fDirectoryPhoton["Direct"]
        NumIntersecting = getNumIntersecting(fDirectoryPhoton["Info"])

        energy, eventMapping, dosePerEvent, _, _ = doseCalculation(tEdepPhoton, pathLength, dosePerEvent,meanKEperEvent, "standalone", chromatinVolume) 
        energy, eventMapping, dosePerEvent, _, _  = doseCalculation(tEdep, pathLength, dosePerEvent,meanKEperEvent, "photon", chromatinVolume) 

        cumulatedEnergyDep = AccumulateEdep(input_tree_photon, cumulatedEnergyDep, T0, T1) #read in photon edep first

        cumulatedEnergyDep = AccumulateEdep(input_tree, cumulatedEnergyDep, T0, T1, eventMapping = eventMapping)

        eventsListDirect, copyListDirect, strandListDirect = Direct(cumulatedEnergyDep, fEMinDamage, fEMaxDamage)
        eventsListIndirect, copyListIndirect, strandListIndirect = Indirect(input_tree_indirect, probIndirect, T0, T1, eventMapping = eventMapping, copyCheck = False)
        LET = "N/A"
        # Path length and KE are not recorded for photon simulations
        pathLength = {e: "N/A" for e in dosePerEvent.keys()}
        meanKEperEvent = {e: "N/A" for e in dosePerEvent.keys()}
    else:
        raise ValueError("Incompatible simulation type")

    numEvt = list(set(eventsListDirect+eventsListIndirect))

    file.close()     
    # clustering
    tempResults = clustering(numEvt,eventsListDirect,copyListDirect,strandListDirect,eventsListIndirect,copyListIndirect, strandListIndirect, continuous)

    clusteringResults = tempResults[0] # strand break number results
    clusterSize = tempResults[1] # DSB cluster size, cluster of size 1-10
    clusterDistance = tempResults[2] # distance in bp between DSB clusters (if more 2 or more)

    events = list(dosePerEvent.keys())
    events.sort()
    eventsWithClusteringResults = [a[0] for a in clusteringResults]

    clusteringGitHash = readClusteringGitHash(builddir)


    # write out results to txt file
    with open(outputFilename, "w") as f:
        f.write("Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,GitHash,ClusteringGitHash\n")
        f.write("{},{},{},{},{},{},{},{},{}\n".format(filename,energy,LET,NumIntersecting,chromatinVolume,numBP,sugarPosFilename, gitHash, clusteringGitHash)) 

        f.write("EventNum,DoseGy,PathLength_nm,MeanKEMeV,TotalSBdirect,SSBdirect,cSSBdirect,DSBdirect,TotalSBindirect,SSBindirect,cSSBindirect,DSBindirect,TotalSBtotal,SSBtotal,cSSBtotal,DSBtotal\n")
        for event in events:
            if event in eventsWithClusteringResults:
                text = str([clusteringResults[eventsWithClusteringResults.index(event)][1:]]).strip("[").strip("]")
            else:
                text = "0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0"           
            f.write("{},{},{},{},{}".format(event,dosePerEvent[event], pathLength[event],meanKEperEvent[event],text))
            f.write('\n')



    with open(str(Path(outputFilename).with_name("DSBclusterSize_" + Path(outputFilename).name)), "w") as f:
        f.write("Filename,EnergyMeV,LET,numEvtIntersectingVolume,chromatinVolume,numBP,sugarPosFilename,SimGitHash,ClusteringGitHash\n")
        f.write("{},{},{},{},{},{},{},{},{}\n".format(filename,energy,LET,NumIntersecting,chromatinVolume,numBP,sugarPosFilename, gitHash, clusteringGitHash)) 
        f.write("EventNum,DoseGy,PathLength_nm,Direct_1_SBperDSBcluster,Direct_2_SBperDSBcluster,Direct_3_SBperDSBcluster,Direct_4_SBperDSBcluster,Direct_5_SBperDSBcluster,Direct_6_SBperDSBcluster,Direct_7_SBperDSBcluster,Direct_8_SBperDSBcluster,Direct_9_SBperDSBcluster,Direct_10_SBperDSBcluster,Indirect_1_SBperDSBcluster,Indirect_2_SBperDSBcluster,Indirect_3_SBperDSBcluster,Indirect_4_SBperDSBcluster,Indirect_5_SBperDSBcluster,Indirect_6_SBperDSBcluster,Indirect_7_SBperDSBcluster,Indirect_8_SBperDSBcluster,Indirect_9_SBperDSBcluster,Indirect_10_SBperDSBcluster,Hybrid_1_SBperDSBcluster,Hybrid_2_SBperDSBcluster,Hybrid_3_SBperDSBcluster,Hybrid_4_SBperDSBcluster,Hybrid_5_SBperDSBcluster,Hybrid_6_SBperDSBcluster,Hybrid_7_SBperDSBcluster,Hybrid_8_SBperDSBcluster,Hybrid_9_SBperDSBcluster,Hybrid_10_SBperDSBcluster,Mixed_1_SBperDSBcluster,Mixed_2_SBperDSBcluster,Mixed_3_SBperDSBcluster,Mixed_4_SBperDSBcluster,Mixed_5_SBperDSBcluster,Mixed_6_SBperDSBcluster,Mixed_7_SBperDSBcluster,Mixed_8_SBperDSBcluster,Mixed_9_SBperDSBcluster,Mixed_10_SBperDSBcluster,Total_1_SBperDSBcluster,Total_2_SBperDSBcluster,Total_3_SBperDSBcluster,Total_4_SBperDSBcluster,Total_5_SBperDSBcluster,Total_6_SBperDSBcluster,Total_7_SBperDSBcluster,Total_8_SBperDSBcluster,Total_9_SBperDSBcluster,Total_10_SBperDSBcluster,Total_DSBdistances\n")
        for event in events:
            if event in eventsWithClusteringResults:
                text = str([clusterSize[eventsWithClusteringResults.index(event)][1:]]).strip("[").strip("]")+", "+ str([d for d in clusterDistance[eventsWithClusteringResults.index(event)]]).replace(",","\t")
            else:
                text = "0, "*49+"0" + ", []"

            f.write("{},{},{},{}".format(event,dosePerEvent[event], pathLength[event], text))
            f.write('\n')

    print("Finished: {}".format(filename))

if __name__ == "__main__":
    pass

