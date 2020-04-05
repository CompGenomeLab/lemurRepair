import scipy.stats as stats
import statistics


human = []
lemur = []
with open("tcr.csv", "r") as tcr:
    for line in tcr:
        h, l = line.strip().split(',')
        human.append(float(h))
        lemur.append(float(l))
    print("Human Median: ", statistics.median(human), "Lemur Median: ", statistics.median(lemur))
    print(stats.mannwhitneyu(human, lemur, alternative='two-sided'))