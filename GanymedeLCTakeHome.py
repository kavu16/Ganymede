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
    peakVolumes: list[float] = []

    def __init__(self, file_name: str):
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

        return
    
    def find_peaks(self, threshold, influence, lag):
        # Find peaks using z-score smoothing
        # If current point within *threshold* times the rolling average of the *lag* then we have detected a peak
        if len(self.rawChromData) < lag + 2:
            print("Not Enough Data")
            return
        
        # Initialize Rolling average and standard deviations
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
            else:
                # Using index to store left/right threshold (can also store time value if needed)
                start = i
                peak_max = self.rawChromData[i][2]
                # Loop to find left and right threshold of the peak
                while i < len(filtered_chrom) and abs(self.rawChromData[i][2] - running_avg) > threshold * running_std:
                    peak_max = max(peak_max, self.rawChromData[i][2])
                    filtered_chrom[i][2] = influence * self.rawChromData[i][2] + (1 - influence) * filtered_chrom[i-1][2]
                    running_avg = np.mean(filtered_chrom[(i-lag+1):i+1, 2:])
                    running_std = np.std(filtered_chrom[(i-lag+1):i+1, 2:])
                    i += 1

                end = i - 1
                self.peakData.append([start, end, peak_max])
        return
        


    def elutionVolume(self) -> float:
        # using Simpson's Quadrature rule: Int(f(x)) a -> b = (b-a)/6[f(a) + 4f(a+b/2) + f(b)]
        # Assuming constant step size of 0.5 s
        # if step size varies, other quadrature methods need to be implemented
        
        total = 0
        for start, end, _ in self.peakData:
            sum = 0
            for i in range(start, end, 2):
                a = self.rawChromData[i][2] if i < len(self.rawChromData) else 0
                m = self.rawChromData[i+1][2] if i + 1 < len(self.rawChromData) else 0
                b = self.rawChromData[i+2][2] if i + 2 < len(self.rawChromData) else 0
                sum += (1/6) * (a + 4 * m + b) # each step is 0.5 in the data, so b-a = 1

            self.peakVolumes.append(sum)
            total += sum

        return total
    
    def visualize(self):
        pass


if __name__ == '__main__':
    file_name = "IgG Vtag 4_ACQUITY FLR ChA.txt"
    chrom_run = ChromatogramRun(file_name)

    chrom_run.find_peaks(3.75, 0, 100)
    print(chrom_run.elutionVolume())