import matplotlib.pyplot as plt

with open("iso.txt","r",encoding="utf8") as file :
        lines = file.readlines()

# print(lines[16].split)
logg = []
logTeff = []

for i in lines[14:-1]:
        logg.append(float(i.split()[8]))
        logTeff.append(float(i.split()[7]))

print(logg[0], logg[-1], logTeff[0], logTeff[-1])


f = plt.figure(figsize=(10, 5))
gs = f.add_gridspec(1)
ax = gs.subplots(sharex=False, sharey=True)

ax.plot(logTeff, logg)
ax.set_xlabel("log T$_{eff}$")
ax.set_ylabel("log g")

plt.show()