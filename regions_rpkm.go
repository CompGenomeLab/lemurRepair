package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

// Region -- Matching regions of human&lemur
type Region struct {
	Identity      float64
	Genic         bool
	LemurRNACount int
	HumanRNACount int
	LemurRNARPKM  float64
	HumanRNARPKM  float64
	LemurSimRPKM  float64
	HumanSimRPKM  float64
	LemurSimCount int
	HumanSimCount int
	LemurXrRPKM   float64
	HumanXrRPKM   float64
	LemurXrCount  int
	HumanXrCount  int
	LemurChr      string
	LemurStart    int
	LemurEnd      int
	LemurLength   int
	HumanChr      string
	HumanStart    int
	HumanEnd      int
	HumanLength   int
}

// Range -- Region start & end coordinates
type Range struct {
	Min, Max int
	Value    string
}

// Ranges -- Range slice
type Ranges []Range

func (r Ranges) Len() int {
	return len(r)
}

func (r Ranges) Less(i, j int) bool {
	return r[i].Min < r[j].Min
}

func (r Ranges) Swap(i, j int) {
	r[i], r[j] = r[j], r[i]
}

//Sort -- Sort Range slice
func (r Ranges) Sort() {
	sort.Sort(r)
}

//Search -- find if given integer in which range
func (r Ranges) Search(v int) string {
	ln := r.Len()
	if i := sort.Search(ln, func(i int) bool {
		return v <= r[i].Max
	}); i < ln {
		if it := &r[i]; v >= it.Min && v <= it.Max {
			return it.Value
		}
	}
	return ""
}

func main() {
	allRegions, lemurRegions, humanRegions := parseRegions("/cta/users/umitakkose/lemur/overlap_wGenes.bed")

	start := time.Now()
	allRegions, lemurSimLineCount := parseLemurSimBed("/cta/users/umitakkose/lemur/random_shuffled_lemur_cutadapt_sorted_plus.bed", allRegions, lemurRegions)
	elapsed := time.Since(start)
	fmt.Printf("Finished Parsing Lemur Sim Bed file in : %s\n", elapsed)

	start = time.Now()
	allRegions, humanSimLineCount := parseHumanSimBed("/cta/users/umitakkose/lemur/random_shuffled_human_cutadapt_sorted_plus.bed", allRegions, humanRegions)
	elapsed = time.Since(start)
	fmt.Printf("Finished Parsing Human Sim Bed file in : %s\n", elapsed)

	start = time.Now()
	allRegions, lemurXrLineCount := parseLemurBed("/cta/groups/adebali/data/lemur/res/CPD_60m_A/Lemur_CPD_GTCCGC_S5_L002_R1_001_cutadapt_sorted.bed", allRegions, lemurRegions)
	elapsed = time.Since(start)
	fmt.Printf("Finished Parsing Lemur Bed file in : %s\n", elapsed)

	start = time.Now()
	allRegions, humanXrLineCount := parseHumanBed("/cta/groups/adebali/data/human/GSM1659156_cutadapt_sorted.bed", allRegions, humanRegions)
	elapsed = time.Since(start)
	fmt.Printf("Finished Parsing Human Bed file in : %s\n", elapsed)

	// filter 0 reads
	start = time.Now()
	allRegions, lemurRNALineCount := parseLemurRNABed("/cta/groups/adebali/data/lemur/reads/Lemur1_cutadapt_sorted.bed", allRegions, lemurRegions)
	elapsed = time.Since(start)
	fmt.Printf("Finished Parsing Lemur RNA Bed file in : %s\n", elapsed)

	start = time.Now()
	allRegions, humanRNALineCount := parseHumanRNABed("/cta/groups/adebali/data/human/rnaseq/SRR3192539_cutadapt_sorted.bed", allRegions, humanRegions)
	elapsed = time.Since(start)
	fmt.Printf("Finished Parsing Human RNA Bed file in : %s\n", elapsed)

	allRegions = remove0RNA(allRegions)

	start = time.Now()
	allRegions = calculateRPKM(allRegions, lemurSimLineCount, humanSimLineCount, humanXrLineCount, lemurXrLineCount, lemurRNALineCount, humanRNALineCount)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating RPKM in : %s\n", elapsed)

	writeCounts("/cta/users/umitakkose/lemur/random_final_go/res.tsv", allRegions)

	writeQuartersIdentity("/cta/users/umitakkose/lemur/random_final_go/res", allRegions)

	writeQuartersRNALemur("/cta/users/umitakkose/lemur/random_final_go/res", allRegions)

	writeQuartersRNAHuman("/cta/users/umitakkose/lemur/random_final_go/res", allRegions)

}

func remove0RNA(allRegions map[string]Region) map[string]Region {
	for name := range allRegions {
		if allRegions[name].HumanRNACount == 0 || allRegions[name].LemurRNACount == 0 {
			delete(allRegions, name)
		} else if allRegions[name].HumanXrCount == 0 || allRegions[name].LemurXrCount == 0 {
			delete(allRegions, name)
		} else if allRegions[name].HumanSimCount == 0 || allRegions[name].LemurSimCount == 0 {
			delete(allRegions, name)
		}
	}

	return allRegions
}

func writeQuartersIdentity(prefix string, allRegions map[string]Region) {

	allRegionsLength := len(allRegions)
	keys := make([]Region, 0, allRegionsLength)
	for name := range allRegions {
		keys = append(keys, allRegions[name])
	}
	sort.Slice(keys, func(i, j int) bool {
		return keys[i].Identity < keys[j].Identity
	})

	quarter1 := prefix + "_identity_quarter1.tsv"
	fileQuarter1, err := os.Create(quarter1)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter1.Close()
	writer1 := bufio.NewWriter(fileQuarter1)
	quarterEnd1 := allRegionsLength / 4
	for i := 0; i < quarterEnd1; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer1.WriteString(line)
		writer1.Flush()
	}

	quarter2 := prefix + "_identity_quarter2.tsv"
	fileQuarter2, err := os.Create(quarter2)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter2.Close()
	writer2 := bufio.NewWriter(fileQuarter2)
	quarterEnd2 := allRegionsLength / 2
	for i := quarterEnd1; i < quarterEnd2; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer2.WriteString(line)
		writer2.Flush()
	}

	quarter3 := prefix + "_identity_quarter3.tsv"
	fileQuarter3, err := os.Create(quarter3)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter3.Close()
	writer3 := bufio.NewWriter(fileQuarter3)
	quarterEnd3 := (allRegionsLength / 4) * 3
	for i := quarterEnd2; i < quarterEnd3; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer3.WriteString(line)
		writer3.Flush()
	}

	quarter4 := prefix + "_identity_quarter4.tsv"
	fileQuarter4, err := os.Create(quarter4)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter4.Close()
	writer4 := bufio.NewWriter(fileQuarter4)
	quarterEnd4 := allRegionsLength
	for i := quarterEnd3; i < quarterEnd4; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer4.WriteString(line)
		writer4.Flush()
	}

}

func writeQuartersRNALemur(prefix string, allRegions map[string]Region) {

	allRegionsLength := len(allRegions)
	keys := make([]Region, 0, allRegionsLength)
	for name := range allRegions {
		keys = append(keys, allRegions[name])
	}
	sort.Slice(keys, func(i, j int) bool {
		return keys[i].LemurRNARPKM < keys[j].LemurRNARPKM
	})

	quarter1 := prefix + "_RNALemur_quarter1.tsv"
	fileQuarter1, err := os.Create(quarter1)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter1.Close()
	writer1 := bufio.NewWriter(fileQuarter1)
	quarterEnd1 := allRegionsLength / 4
	for i := 0; i < quarterEnd1; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer1.WriteString(line)
		writer1.Flush()
	}

	quarter2 := prefix + "_RNALemur_quarter2.tsv"
	fileQuarter2, err := os.Create(quarter2)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter2.Close()
	writer2 := bufio.NewWriter(fileQuarter2)
	quarterEnd2 := allRegionsLength / 2
	for i := quarterEnd1; i < quarterEnd2; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer2.WriteString(line)
		writer2.Flush()
	}

	quarter3 := prefix + "_RNALemur_quarter3.tsv"
	fileQuarter3, err := os.Create(quarter3)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter3.Close()
	writer3 := bufio.NewWriter(fileQuarter3)
	quarterEnd3 := (allRegionsLength / 4) * 3
	for i := quarterEnd2; i < quarterEnd3; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer3.WriteString(line)
		writer3.Flush()
	}

	quarter4 := prefix + "_RNALemur_quarter4.tsv"
	fileQuarter4, err := os.Create(quarter4)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter4.Close()
	writer4 := bufio.NewWriter(fileQuarter4)
	quarterEnd4 := allRegionsLength
	for i := quarterEnd3; i < quarterEnd4; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer4.WriteString(line)
		writer4.Flush()
	}

}

// change file writes
func writeQuartersRNAHuman(prefix string, allRegions map[string]Region) {

	allRegionsLength := len(allRegions)
	keys := make([]Region, 0, allRegionsLength)
	for name := range allRegions {
		keys = append(keys, allRegions[name])
	}
	sort.Slice(keys, func(i, j int) bool {
		return keys[i].HumanRNARPKM < keys[j].HumanRNARPKM
	})

	quarter1 := prefix + "_RNAHuman_quarter1.tsv"
	fileQuarter1, err := os.Create(quarter1)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter1.Close()
	writer1 := bufio.NewWriter(fileQuarter1)
	quarterEnd1 := allRegionsLength / 4
	for i := 0; i < quarterEnd1; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer1.WriteString(line)
		writer1.Flush()
	}

	quarter2 := prefix + "_RNAHuman_quarter2.tsv"
	fileQuarter2, err := os.Create(quarter2)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter2.Close()
	writer2 := bufio.NewWriter(fileQuarter2)
	quarterEnd2 := allRegionsLength / 2
	for i := quarterEnd1; i < quarterEnd2; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer2.WriteString(line)
		writer2.Flush()
	}

	quarter3 := prefix + "_RNAHuman_quarter3.tsv"
	fileQuarter3, err := os.Create(quarter3)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter3.Close()
	writer3 := bufio.NewWriter(fileQuarter3)
	quarterEnd3 := (allRegionsLength / 4) * 3
	for i := quarterEnd2; i < quarterEnd3; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer3.WriteString(line)
		writer3.Flush()
	}

	quarter4 := prefix + "_RNAHuman_quarter4.tsv"
	fileQuarter4, err := os.Create(quarter4)
	if err != nil {
		log.Fatal(err)
	}
	defer fileQuarter4.Close()
	writer4 := bufio.NewWriter(fileQuarter4)
	quarterEnd4 := allRegionsLength
	for i := quarterEnd3; i < quarterEnd4; i++ {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t\n", keys[i].HumanXrRPKM, keys[i].HumanXrCount, keys[i].LemurXrRPKM, keys[i].LemurXrCount, keys[i].HumanSimRPKM, keys[i].HumanSimCount, keys[i].LemurSimRPKM, keys[i].LemurSimCount, keys[i].HumanRNARPKM, keys[i].HumanRNACount, keys[i].LemurRNARPKM, keys[i].LemurRNACount, keys[i].Genic, keys[i].Identity)
		writer4.WriteString(line)
		writer4.Flush()
	}

}

func writeCounts(outfile string, allRegions map[string]Region) {
	file, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	writer := bufio.NewWriter(file)

	for _, el := range allRegions {
		line := fmt.Sprintf("%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%f\t%d\t%t\t%f\t%s\t%d\t%d\t%s\t%d\t%d\n", el.HumanXrRPKM, el.HumanXrCount, el.LemurXrRPKM, el.LemurXrCount, el.HumanSimRPKM, el.HumanSimCount, el.LemurSimRPKM, el.LemurSimCount, el.HumanRNARPKM, el.HumanRNACount, el.LemurRNARPKM, el.LemurRNACount, el.Genic, el.Identity, el.LemurChr, el.LemurStart, el.LemurEnd, el.HumanChr, el.HumanStart, el.HumanEnd)
		writer.WriteString(line)
		writer.Flush()
	}

}

func calculateRPKM(allRegions map[string]Region, lemurSimLineCount int, humanSimLineCount int, lemurXrLineCount int, humanXrLineCount int, lemurRNALineCount int, humanRNALineCount int) map[string]Region {

	lemurSimLineCount = 0
	humanSimLineCount = 0
	lemurXrLineCount = 0
	humanXrLineCount = 0
	lemurRNALineCount = 0
	humanRNALineCount = 0

	for _, el := range allRegions {
		lemurSimLineCount += el.LemurSimCount
		humanSimLineCount += el.HumanSimCount
		lemurXrLineCount += el.LemurXrCount
		humanXrLineCount += el.HumanXrCount
		lemurRNALineCount += el.LemurRNACount
		humanRNALineCount += el.HumanRNACount
	}

	pmsflemurSim := float64(lemurSimLineCount) / 1000000
	pmsfhumanSim := float64(humanSimLineCount) / 1000000
	pmsflemurXr := float64(lemurXrLineCount) / 1000000
	pmsfhumanXr := float64(humanXrLineCount) / 1000000
	pmsflemurRNA := float64(lemurRNALineCount) / 1000000
	pmsfhumanRNA := float64(humanRNALineCount) / 1000000
	fmt.Println(pmsflemurSim, pmsfhumanSim, pmsflemurXr, pmsfhumanXr, pmsflemurRNA, pmsfhumanRNA)

	for i, el := range allRegions {
		el.LemurSimRPKM = (float64(el.LemurSimCount) / pmsflemurSim) / float64(el.LemurLength)
		el.HumanSimRPKM = (float64(el.HumanSimCount) / pmsfhumanSim) / float64(el.HumanLength)
		el.LemurXrRPKM = (float64(el.LemurXrCount) / pmsflemurXr) / float64(el.LemurLength)
		el.HumanXrRPKM = (float64(el.HumanXrCount) / pmsfhumanXr) / float64(el.HumanLength)
		el.LemurRNARPKM = (float64(el.LemurRNACount) / pmsflemurRNA) / float64(el.LemurLength)
		el.HumanRNARPKM = (float64(el.HumanRNACount) / pmsfhumanRNA) / float64(el.HumanLength)

		allRegions[i] = el
	}

	return allRegions
}

func parseHumanRNABed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.HumanRNACount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount
}

func parseLemurRNABed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := "chr" + fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.LemurRNACount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount
}

func parseLemurSimBed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := "chr" + fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.LemurSimCount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount

}

func parseHumanSimBed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.HumanSimCount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount

}

func parseLemurBed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := "chr" + fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.LemurXrCount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount

}

func parseHumanBed(infile string, allRegions map[string]Region, regions map[string]Ranges) (map[string]Region, int) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	lineCount := 0

	for i := range regions {
		regions[i].Sort()
	}

	for scanner.Scan() {

		line := scanner.Text()
		fields := strings.Split(line, "\t")

		chr := fields[0]
		start, _ := strconv.Atoi(fields[1])
		end, _ := strconv.Atoi(fields[2])
		length := end - start
		middle := start + length/2

		region := regions[chr].Search(middle)
		if region != "" {
			lineCount++
			reg := allRegions[region]
			reg.HumanXrCount++
			allRegions[region] = reg
		}
	}

	return allRegions, lineCount

}

func parseRegions(infile string) (map[string]Region, map[string]Ranges, map[string]Ranges) {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	regions := make(map[string]Region)
	humanRegions := make(map[string]Ranges)
	lemurRegions := make(map[string]Ranges)

	nameHash := 0

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Split(line, "\t")
		human := strings.Split(fields[3], ":")

		identity, _ := strconv.ParseFloat(fields[6], 64)
		var genic bool
		if fields[7] == "." {
			genic = false
		} else {
			genic = true
		}

		name := "range" + strconv.Itoa(nameHash)

		lemurChr := fields[0]
		lemurStart, _ := strconv.Atoi(fields[1])
		lemurEnd, _ := strconv.Atoi(fields[2])
		lemurLength := lemurEnd - lemurStart

		humanChr := human[0][4:]
		humanStart, _ := strconv.Atoi(human[1])
		humanEnd, _ := strconv.Atoi(human[2])
		humanLength := humanEnd - humanStart

		lemurRegions[lemurChr] = append(lemurRegions[lemurChr],
			Range{Min: lemurStart,
				Max:   lemurEnd,
				Value: name})

		humanRegions[humanChr] = append(humanRegions[humanChr],
			Range{Min: humanStart,
				Max:   humanEnd,
				Value: name})

		region := Region{
			Identity:      identity,
			LemurSimRPKM:  0,
			HumanSimRPKM:  0,
			LemurSimCount: 0,
			HumanSimCount: 0,
			LemurXrRPKM:   0,
			HumanXrRPKM:   0,
			LemurXrCount:  0,
			HumanXrCount:  0,
			LemurRNACount: 0,
			HumanRNACount: 0,
			LemurRNARPKM:  0,
			HumanRNARPKM:  0,
			Genic:         genic,
			LemurChr:      lemurChr,
			LemurStart:    lemurStart,
			LemurEnd:      lemurEnd,
			LemurLength:   lemurLength,
			HumanChr:      humanChr,
			HumanStart:    humanStart,
			HumanEnd:      humanEnd,
			HumanLength:   humanLength}

		regions[name] = region

		nameHash++

	}

	return regions, lemurRegions, humanRegions
}
