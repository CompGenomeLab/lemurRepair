import argparse


def parsebiomart(infile, outpath):
    print('Parsing BIOMART file')
    regions = []
    chr_list = ['X', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32']
    with open(infile, 'r') as infile_handler:
        for line in infile_handler:
            clash = False
            line = line.strip().split('\t')
            if line[2] in chr_list:
                if line[6] == 'protein_coding':
                    for i, region in enumerate(regions):
                        r_mid = int(line[3]) + ((int(line[4]) - int(line[3])) // 2)
                        r_len = int(line[4]) - int(line[3])
                        if region['g_name'] == f'{line[0]}':
                            clash = True
                            break
                        elif abs(region['mid'] - r_mid) <= ((region['len']//2)+20000+(r_len//2)):
                            clash = True
                            regions.pop(i)
                            break
                    if not clash:
                        regions.append({
                            'g_name': f'{line[0]}',
                            'name': f'{line[1]}',
                            'chr': f'chr{line[2]}',
                            'tss': int(line[3]),
                            'tes': int(line[4]),
                            'len': int(line[4]) - int(line[3]),
                            'mid': int(line[3]) + ((int(line[4]) - int(line[3])) // 2),
                            'strand': '+' if f'{line[5]}'=='1' else '-',
                        })


    for region in regions:
        with open(f'{outpath}_tss.bed', 'w') as tss_bed:
            with open(f'{outpath}_tes.bed', 'w') as tes_bed:
                for line in infile_handler:
                    line = line.strip().split('\t')
                    tss_bed.write(f'{region['chr']}\t{int(region['tss'])-10050}\t{int(region['tss'])+10050}\t{region['g_name']}\t{region['len']}\t{region['strand']}\n')
                    tes_bed.write(f'{region['chr']}\t{int(region['tes'])-10050}\t{int(region['tes'])+10050}\t{region['g_name']}\t{region['len']}\t{region['strand']}\n')

def main:
    parser = argparse.ArgumentParser()
	parser.add_argument(
		'--biomart',
		action='store',
		dest='biomart',
		help='biomart export with gene start and end locations'
	)

    parser.add_argument(
		'--out',
		action='store',
		dest='out',
		help='output path for tss and tes bed files'
	)

	args = parser.parse_args()

    parsebiomart(args.biomart, args.out)


if __name__ == "__main__":
    main()
