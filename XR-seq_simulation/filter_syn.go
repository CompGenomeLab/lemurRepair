package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"log"
	"math/rand"
	"os"
	"sort"
	"time"
)

// Record -- single fastq record
type Record struct {
	Title      string
	Seq        string
	Identifier string
	Qual       string
	Score      float64
}

// Nuc -- frequency of single nucleotide
type Nuc struct {
	A float64
	C float64
	T float64
	G float64
}

func main() {

	xrfastqpath := os.Args[1]
	simulatedfastqpath := os.Args[2]
	filteredfastqpath := os.Args[3]

	fmt.Println("Started reading xr fastq")
	start := time.Now()
	xrfastq := readFastq(xrfastqpath)
	elapsed := time.Since(start)
	fmt.Printf("Finished reading xr fastq in: %s\n", elapsed)

	fmt.Println("Started calculating xr frequency")
	start = time.Now()
	frequency := calculateFrequency(xrfastq)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating xr frequency in: %s\n", elapsed)

	fmt.Println("Started reading simulated fastq")
	start = time.Now()
	fastqWscores := readAndScoreFastq(simulatedfastqpath, frequency)
	elapsed = time.Since(start)
	fmt.Printf("Finished reading simulated fastq in: %s\n", elapsed)

	fmt.Println("Started calculating simulated frequency")
	start = time.Now()
	simulatedFrequency := calculateFrequency(fastqWscores)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating simulated frequency in: %s\n", elapsed)

	fmt.Println("Frequency of simulated reads")
	printFrequency(simulatedFrequency)

	fmt.Println("Started sorting simulated reads")
	start = time.Now()

	sort.Slice(fastqWscores, func(i, j int) bool {
		return fastqWscores[i].Score > fastqWscores[j].Score
	})

	elapsed = time.Since(start)
	fmt.Printf("Finished sorting simulated reads in: %s\n", elapsed)

	fmt.Println("Started calculating filtered frequency")
	start = time.Now()
	filteredFrequency := calculateFrequency(fastqWscores[:10000000])
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating filtered frequency in: %s\n", elapsed)

	fmt.Println("Frequency of filtered reads")
	printFrequency(filteredFrequency)

	fmt.Println("Writing filtered fastq to file")
	start = time.Now()
	writeFastq(filteredfastqpath, fastqWscores[:10000000])
	elapsed = time.Since(start)
	fmt.Printf("Finished writing filtered fastq in: %s\n", elapsed)
}

func printFrequency(freq []Nuc) {
	b, err := json.Marshal(freq)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(string(b))
}

func calculateFrequency(records []Record) []Nuc {

	length := float64(len(records))

	counts := make([][]int, 26)
	for i := 0; i < 26; i++ {
		counts[i] = make([]int, 4)
	}
	for _, r := range records {
		for i, s := range r.Seq {
			if s == 'A' {
				counts[i][0]++
			} else if s == 'C' {
				counts[i][1]++
			} else if s == 'T' {
				counts[i][2]++
			} else if s == 'G' {
				counts[i][3]++
			}
		}
	}
	var f []Nuc
	for _, base := range counts {
		n := Nuc{
			A: float64(base[0]) / length,
			C: float64(base[1]) / length,
			T: float64(base[2]) / length,
			G: float64(base[3]) / length}
		f = append(f, n)
	}

	return f
}

func calculateScore(seq string, freq []Nuc) float64 {
	var score float64 = 1.0
	for i, s := range freq {
		if seq[i] == 'A' {
			score *= s.A
		} else if seq[i] == 'C' {
			score *= s.C
		} else if seq[i] == 'T' {
			score *= s.T
		} else if seq[i] == 'G' {
			score *= s.G
		}
	}

	return score
}

func writeFastq(outfile string, records []Record) {
	file, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	writer := bufio.NewWriter(file)

	for _, record := range records {
		line := fmt.Sprintf("%s\n%s\n%s\n%s\n", record.Title, record.Seq, record.Identifier, record.Qual)
		writer.WriteString(line)
	}

}

func readFastq(infile string) []Record {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	var fastq []Record

	for scanner.Scan() {

		title := scanner.Text()
		scanner.Scan()
		seq := scanner.Text()
		scanner.Scan()
		identifier := scanner.Text()
		scanner.Scan()
		qual := scanner.Text()

		rec := Record{
			Title:      title,
			Seq:        seq,
			Identifier: identifier,
			Qual:       qual}
		if len(rec.Seq) == 26 {
			fastq = append(fastq, rec)
		}
	}

	return fastq
}

func readAndScoreFastq(infile string, frequency []Nuc) []Record {
	rand.Seed(time.Now().UnixNano())
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	var fastq []Record

	for scanner.Scan() {

		title := scanner.Text()
		scanner.Scan()
		seq := scanner.Text()
		scanner.Scan()
		identifier := scanner.Text()
		scanner.Scan()
		qual := scanner.Text()
		score := calculateScore(seq, frequency)

		rec := Record{
			Title:      title,
			Seq:        seq,
			Identifier: identifier,
			Qual:       qual,
			Score:      score}
		if len(rec.Seq) == 26 {
			fastq = append(fastq, rec)
		}
	}

	return fastq
}
