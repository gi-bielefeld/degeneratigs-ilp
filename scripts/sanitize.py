from fasta import read_mfasta,FastaSequence
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("fasta")
args = parser.parse_args()

entries = read_mfasta(args.fasta)

mod_entries = []
for entry in entries:
    rogue_chars = set(entry.sequence).difference("ACGT")
    if len(rogue_chars)==0:
        mod_entries.append(entry)
        continue
    seq_parts = []
    for r in rogue_chars:
        new_parts = []
        for part in seq_parts:
            new_parts.extend(part.split(r))
        seq_parts=new_parts
    for i,part in enumerate(seq_parts):
        mod_entries.append(FastaSequence(entry.id+f" part{i}"))
for entry in mod_entries:
    print(entry.to_str())