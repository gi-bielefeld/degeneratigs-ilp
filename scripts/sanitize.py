from fasta import read_mfasta,FastaSequence
from argparse import ArgumentParser
import sys
parser = ArgumentParser()

parser.add_argument("fasta")
parser.add_argument("--only-first",type=int)
args = parser.parse_args()

entries = read_mfasta(args.fasta)

if args.only_first:
    entries=entries[:args.only_first]
mod_entries = []
for entry in entries:
    
    rogue_chars = set(entry.sequence).difference("ACGT")
    if len(rogue_chars)==0:
        mod_entries.append(entry)
        continue
    seq_parts = [entry.sequence]
    for r in rogue_chars:
        new_parts = []
        for part in seq_parts:
            new_parts.extend(part.split(r))
        seq_parts=new_parts
    for i,part in enumerate(seq_parts):
        sq = FastaSequence(entry.id+f" part{i}",part)
        mod_entries.append(sq)
for entry in mod_entries:
    print(entry.to_str())