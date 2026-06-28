import fasta
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("fasta")
args = parser.parse_args()



entries = fasta.read_mfasta(args.fasta)

print(len(entries),sum(len(e.sequence) for e in entries))