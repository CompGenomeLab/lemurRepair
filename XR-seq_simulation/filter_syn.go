package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"log"
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

	oligomerLengthPtr := flag.Int("length", 26, "length of oligomer")
	numReadsPtr := flag.Int("numReads", 30000000, "number of reads in filtered fastq, must be smaller than inputFile's number of reads")
	modelFastaPtr := flag.String("in", "./model.fa", "name of input fasta file, to extract frequencies")
	synFastqPtr := flag.String("sim", "./simulated.fastq", "name of simulated input fastq file")
	filteredFastqPtr := flag.String("out", "./filtered.fastq", "name of filtered output fastq file")

	flag.Parse()

	fmt.Println("Started reading input fasta")
	start := time.Now()
	modelFasta := readFasta(*modelFastaPtr, *oligomerLengthPtr)
	elapsed := time.Since(start)
	fmt.Printf("Finished reading input fasta in: %s\n", elapsed)

	fmt.Println("Started calculating input frequency")
	start = time.Now()
	frequency := calculateFrequency(modelFasta, *oligomerLengthPtr)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating input frequency in: %s\n", elapsed)

	fmt.Println("Started reading simulated fastq")
	start = time.Now()
	fastqWscores := readAndScoreFastq(*synFastqPtr, frequency, *oligomerLengthPtr)
	elapsed = time.Since(start)
	fmt.Printf("Finished reading simulated fastq in: %s\n", elapsed)

	fmt.Println("Started calculating simulated frequency")
	start = time.Now()
	simulatedFrequency := calculateFrequency(fastqWscores, *oligomerLengthPtr)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating simulated frequency in: %s\n", elapsed)

	fmt.Println("Frequency of simulated reads")
	printFrequency(simulatedFrequency, "./starting_frequency.json")

	fmt.Println("Started sorting simulated reads")
	start = time.Now()

	sort.Slice(fastqWscores, func(i, j int) bool {
		return fastqWscores[i].Score > fastqWscores[j].Score
	})

	elapsed = time.Since(start)
	fmt.Printf("Finished sorting simulated reads in: %s\n", elapsed)

	fmt.Println("Started calculating filtered frequency")
	start = time.Now()
	filteredFrequency := calculateFrequency(fastqWscores[:*numReadsPtr], *oligomerLengthPtr)
	elapsed = time.Since(start)
	fmt.Printf("Finished calculating filtered frequency in: %s\n", elapsed)

	fmt.Println("Frequency of filtered reads")
	printFrequency(filteredFrequency, "./filtered_frequency.json")

	fmt.Println("Writing filtered fastq to file")
	start = time.Now()
	writeFastq(*filteredFastqPtr, fastqWscores[:*numReadsPtr])
	elapsed = time.Since(start)
	fmt.Printf("Finished writing filtered fastq in: %s\n", elapsed)
}

func printFrequency(freq []Nuc, outfile string) {
	b, err := json.Marshal(freq)
	if err != nil {
		log.Fatal(err)
	}
	file, err := os.Create(outfile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	writer := bufio.NewWriter(file)
	writer.WriteString(string(b))
}

func calculateFrequency(records []Record, ol int) []Nuc {

	length := float64(len(records))

	counts := make([][]int, ol)
	for i := 0; i < ol; i++ {
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

func readFastq(infile string, ol int) []Record {
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
		if len(rec.Seq) == ol {
			fastq = append(fastq, rec)
		}
	}

	return fastq
}

func readFasta(infile string, ol int) []Record {
	file, err := os.Open(infile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	var fasta []Record

	for scanner.Scan() {

		title := scanner.Text()
		scanner.Scan()
		seq := scanner.Text()

		rec := Record{
			Title: title,
			Seq:   seq}
		if len(rec.Seq) == ol {
			fasta = append(fasta, rec)
		}
	}

	return fasta
}

func readAndScoreFastq(infile string, frequency []Nuc, ol int) []Record {
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
		if len(rec.Seq) == ol {
			fastq = append(fastq, rec)
		}
	}

	return fastq
}
