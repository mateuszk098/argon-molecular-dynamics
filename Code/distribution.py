# Momentum distribution
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("../Out/p0_init.txt", skiprows = 2, header=None, sep = '\t')
p0_abs = data[4]
mean_p0_abs = int(p0_abs.mean())

plt.style.use('bmh')
plt.figure(figsize=(8, 5.5), dpi = 100)
plt.subplots_adjust(bottom=0.17, top=0.95)

plt.hist(p0_abs, density = True, color = 'navy', edgecolor = 'black', bins = 50, rwidth = 0.75, alpha = 0.75)
# plt.axvline(mean_p0_abs, color="orange", label="Mean Momentum")
# plt.grid(color = 'gainsboro', linestyle = '-.', linewidth = 0.5)
plt.xlabel("Absolute Momentum", fontname = "Times New Roman", fontsize=16)
plt.ylabel("Momentum Probability", fontname = "Times New Roman", fontsize=16)
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
plt.tick_params(which = 'both', direction='in', top=True, right=True, length = 5, width = 1)
plt.yticks(fontname = "Times New Roman")
plt.xticks(ticks = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22], fontname = "Times New Roman")
plt.title("Absolute Momentum of 15625 Molecules at Temperature 100K", fontname = "Times New Roman", fontsize=16)

plt.show()