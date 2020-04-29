import csv
import vary_vol
import numpy as np

# FCO2 = np.linspace(1501.1, 97.7, 100)
# FCO2 = (1501.1, 1106.1, 840.6, 525.4, 467, 404.1, 350.3, 344.5, 289.9, 254.7, 227.2, 195.5, 161.7, 131.3, 110.6, 97.7)
# FCO2 = (1501.1, 1106.1, 840.6, 525.4, 467, 404.1, 350.3, 344.5, 289.9, 254.7, 227.2, 195.5, 161.7, 97.7)

FCO2 = np.linspace(0.5, 5, 100)
Z = np.zeros((len(FCO2), 9))
for i in range(len(FCO2)):
    Z[i] = vary_vol.main(FCO2[i])

# WoFCO2, W, FCO2, XCO2, Ntubes, heat, Pfinal, Re, COfinal = Z

with open('rwgs_results2.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["index", "WoFCO2", "W", "L", "XCO2", "Ntubes", "heat", "deltaP", "Re", "COfinal"])
    for i in range(len(FCO2)):
        writer.writerow([str(i + 1), Z[i][0], Z[i][1], Z[i][2], Z[i][3], Z[i][4], Z[i][5], Z[i][6], Z[i][7], Z[i][8]])
