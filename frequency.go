package main

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"log"
	"os"
)

// Record -- single fasta record
type Record struct {
	Title string
	Seq   string
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
	fastaPtr := flag.String("in", "./input.fa", "path of input fasta file which we'll be calculating frequencies")
	outPtr := flag.String("out", ".frequency.json", "path of frequency output file")

	flag.Parse()

	fmt.Println("Reading and parsing input fasta")
	fasta := readFasta(*fastaPtr, *oligomerLengthPtr)
	fmt.Println("Calculating frequencies")
	frequency := calculateFrequency(fasta, *oligomerLengthPtr)
	printFrequency(frequency, *outPtr)
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

		if len(seq) == ol {
			fasta = append(fasta, rec)
		}

	}

	return fasta
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
	for _, b := range counts {
		n := Nuc{
			A: float64(b[0]) / length,
			C: float64(b[1]) / length,
			T: float64(b[2]) / length,
			G: float64(b[3]) / length}
		f = append(f, n)
	}

	return f
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
