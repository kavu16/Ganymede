# Ganymede - Liquid Chromatography 
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt

@dataclass
class InjectionData:
    dataVault: str
    injection: str
    injectionNum: int
    position: str
    comment: str
    processMethod: str
    instrumentMethod: str
    type: str
    status: str
    injectionDate: str
    injectionTime: str
    injectionVolume: float
    dilutionFactor: float
    weight: float

@dataclass
class ChromatogramData:
    timeMin: float
    timeMax: float
    dataPoints: int
    detector: str
    genDataSystem: str
    expDataSystem: str
    operator: str
    signalQuantity: str
    sigUnit: str
    sigMin: float
    sigMax: float
    channel: str
    driverName: str
    channelType: str
    minStep: float
    maxStep: float
    avgStep: float

@dataclass
class SignalData:
    signalInfo: str

class ChromatogramRun:
    injData: InjectionData
    chromData: ChromatogramData
    sigData: SignalData

    rawChromData: list[list[float]] = []
    peakData: list[list[float]] = []

    def __init__(self) -> None:
        pass
    
    def find_peaks(self, threshold, influence, lag):
        if len(self.rawChromData) < lag + 2:
            print("Not Enough Data")
            return
        filtered_chrom = np.array(self.rawChromData)
        running_avg = np.mean(filtered_chrom[0:lag, 2:])
        running_std = np.std(filtered_chrom[0:lag, 2:])
        pcount = 0
        ncount = 0

        i = lag
        while i < len(filtered_chrom):
            if abs(self.rawChromData[i][2] - running_avg) <= threshold * running_std:
                running_avg = np.mean(filtered_chrom[(i-lag+1):i+1, 2:])
                running_std = np.std(filtered_chrom[(i-lag+1):i+1, 2:])
                i += 1
                filtered_chrom[i][2] = self.rawChromData[i][2]
                ncount += 1
            else:
                start = self.rawChromData[i][0]
                peak_max = self.rawChromData[i][2]
                while i < len(filtered_chrom) and abs(self.rawChromData[i][2] - running_avg) > threshold * running_std:
                    peak_max = max(peak_max, self.rawChromData[i][2])
                    filtered_chrom[i][2] = influence * self.rawChromData[i][2] + (1 - influence) * filtered_chrom[i-1][2]
                    running_avg = np.mean(filtered_chrom[(i-lag+1):i+1, 2:])
                    running_std = np.std(filtered_chrom[(i-lag+1):i+1, 2:])
                    i += 1
                    pcount += 1

                end = self.rawChromData[i-1][0]
                self.peakData.append([start, end, peak_max])

        print("Number of points in peaks =", pcount)
        print("Number of points not in peaks =", ncount)
        return
        


    def elutionVolume(self) -> float:
        # vol = 0
        # for peak in self.peaks:
        #     vol += (left + 2 * height + right) * step / 2

        # return vol
        pass

if __name__ == '__main__':
    file_name = "IgG Vtag 1_ACQUITY FLR ChA.txt"
    with open(file_name) as file:
        data = file.read().splitlines()


    inj_data, chrom_data, sig_param = data[3:18], data[19:37], data[38:40]
    
    raw_data = data[43:]

    inj_data = InjectionData(
        inj_data[1].split("Data Vault")[1].strip(),
        inj_data[2].split("Injection")[1].strip(),
        int(inj_data[3].split("Injection Number")[1].strip()),
        inj_data[4].split("Position")[1].strip(),
        inj_data[5].split("Comment")[1].strip(),
        inj_data[6].split("Processing Method")[1].strip(),
        inj_data[7].split("Instrument Method")[1].strip(),
        inj_data[8].split("Type")[1].strip(),
        inj_data[9].split("Status")[1].strip(),
        inj_data[10].split("Injection Date")[1].strip(),
        inj_data[11].split("Injection Time")[1].strip(),
        float(inj_data[12].split("Injection Volume")[1].strip().split()[1]),
        float(inj_data[13].split("Dilution Factor")[1].strip()),
        float(inj_data[14].split("Weight")[1].strip())
    )

    chrom_data = ChromatogramData(
        float(chrom_data[1].split("Time Min.")[1].strip().split('\t')[1]),
        float(chrom_data[2].split("Time Max.")[1].strip().split('\t')[1]),
        int(chrom_data[3].split("Data Points")[1].strip()),
        chrom_data[4].split("Detector")[1].strip(),
        chrom_data[5].split("Generating Data System")[1].strip(),
        chrom_data[6].split("Exporting Data System")[1].strip(),
        chrom_data[7].split("Operator")[1].strip(),
        chrom_data[8].split("Signal Quantity")[1].strip(),
        chrom_data[9].split("Signal Unit")[1].strip(),
        float(chrom_data[10].split("Signal Min.")[1].strip()),
        float(chrom_data[11].split("Signal Max.")[1].strip()),
        chrom_data[12].split("Channel")[1].strip(),
        chrom_data[13].split("Driver Name")[1].strip(),
        chrom_data[14].split("Channel Type")[1].strip(),
        float(chrom_data[15].split("Min. Step")[1].strip().split()[1]),
        float(chrom_data[16].split("Max. Step")[1].strip().split()[1]),
        float(chrom_data[17].split("Average Step")[1].strip().split()[1])
    )

    sig_param = SignalData(
        sig_param[1].split("Signal Info")[1].strip()
    )

    chrom_run = ChromatogramRun
    chrom_run.injData = inj_data
    chrom_run.chromData = chrom_data
    chrom_run.sigData = sig_param

    chrom_run.rawChromData.append([0.0, 0.0, 0.0])
    for tsv in raw_data[1:]:
        time, step, val = tsv.split()
        chrom_run.rawChromData.append([float(time), float(step), float(val)])

    chrom_run.find_peaks(chrom_run, 3.75, 0, 100)
    avgs = []

    raw = np.array(chrom_run.rawChromData)
    for i in range(60, len(raw)):
        avgs.append([raw[i][0], np.median(raw[i-60+1:i+1, 2:])])
    avgs = np.array(avgs)
    signal = np.array(chrom_run.peakData)
    plt.plot(raw[:, :1], raw[:, 2:])
    plt.plot(avgs[:, :1], avgs[:, 1:])
    plt.scatter(signal[:, :1], signal[:, 2:])
    plt.show()