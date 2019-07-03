import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

def generatePlot():
  with open("scoredPeps.csv") as f:
    peptides = f.read()
  peptides = peptides.split("\n")
  peptides = [x.split(",") for x in peptides]

  cat1 = []
  cat2 = []
  cat3 = []
  cat4 = []

  for pep in peptides:
    if pep[0] == "1":
      cat1.append((float(pep[2]), float(pep[3])))
    elif pep[0] == "2":
      cat2.append((float(pep[2]), float(pep[3])))
    elif pep[0] == "3":
      cat3.append((float(pep[2]), float(pep[3])))
    elif pep[0] == "4":
      cat4.append((float(pep[2]), float(pep[3])))

  cat3set = set(cat3)
  cat4set = set(cat4)
  cat3 = list(cat3set - cat4set)
  print(len(cat3))
  print(len(cat4))

  # abovePrior = 0
  # aboveGlobal = 0
  # aboveBoth = 0
  # for pep in cat3:
  #   #print(str(pep[0]) + " " + str(pep[1]))
  #   #input()
  #   if pep[0] >= 0.1 and pep[1] >= 0.1:
  #     aboveBoth += 1
  #     abovePrior += 1
  #     aboveGlobal += 1
  #   elif pep[0] >= 0.1:
  #     aboveGlobal +=1
  #   elif pep[1] >= 0.1:
  #     abovePrior += 1

  # print(abovePrior)
  # print(aboveGlobal)
  # print(aboveBoth)

  # print()
  # print(len(cat1))
  # print(len(cat2))
  # print(len(cat3))

  plt.xlabel("Prior Probabilities")
  plt.ylabel("Global Matching Score")
  #plt.plot([x[1] for x in cat3], [x[0] for x in cat3], "ro", markersize=2, marker="x", label="Constrained De Novo Sequencing Only")
  #plt.plot([x[1] for x in cat1], [x[0] for x in cat1], "bo", markersize=2, marker="x", label="Overlap")
  #plt.plot([x[1] for x in cat2], [x[0] for x in cat2], "go", markersize=2, marker="x", label="Database Search Only")
  plt.scatter([x[1] for x in cat4], [x[0] for x in cat4], c="#55aaffff", marker='x', s=1)
  plt.scatter([x[1] for x in cat3], [x[0] for x in cat3], c="#000000ff", marker='x', s=1)
  plt.scatter([x[1] for x in cat1], [x[0] for x in cat1], c="#ffaa00ff", marker='x', s=1)
  plt.scatter([x[1] for x in cat2], [x[0] for x in cat2], c="#aa0000ff", marker='x', s=1)

  blue_patch = mpatches.Patch(color="#aa0000ff", label="Database Search Only")
  red_patch = mpatches.Patch(color="#000000ff", label="Constrained De Novo Sequencing Only")
  green_patch = mpatches.Patch(color="#ffaa00ff", label="Overlap")
  black_patch = mpatches.Patch(color="#55aaffff", label="Constrained De Novo Similar to Human Proteins")

  plt.legend(handles=[blue_patch, red_patch, green_patch, black_patch])

  #plt.show()
  plt.savefig('scatterDatabaseDeNovoOverlapWith4.png', dpi=300)


if __name__ == "__main__":
  generatePlot()
