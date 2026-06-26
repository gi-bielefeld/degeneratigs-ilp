from typing import List

class FastaSequence:
    """
    A class to represent a FASTA sequence entry.
    
    Attributes:
        id (str): The identifier (header until the first space)
        header (str): The full header line without the '>' character
        sequence (str): The biological sequence (DNA, RNA, or protein)
    """
    id : str
    header: str
    sequence: str
    def __init__(self,header,sequence):
        """        
        Args:
            header (str): The header line (with or without the '>' character)
            sequence (str): The sequence string (may contain whitespace/newlines)
        
        Example:
            >>> seq = FastaSequence(">sp|P12345|PROT_HUMAN Protein name", "ACGTACGT")
            >>> seq.id
            'sp|P12345|PROT_HUMAN'
            >>> seq.header
            'sp|P12345|PROT_HUMAN Protein name'
        """
        self.header = header.lstrip('>')
        self.id = self.header.split()[0]
        self.sequence = ''.join(sequence.split())
        
    def __len__(self):
        return len(self.sequence)
        
    def to_str(self,max_columns:int=60):
        """       
        This method wraps the sequence (but not the header!) at a specified column 
        width.
        
        Args:
            max_columns (int): Maximum number of sequence characters per line.        
        Returns:
            str: Formatted FASTA string with:
                 - Header line starting with '>'
                 - Sequence wrapped at max_columns width
        
        Example:
            >>> seq = FastaSequence(">test1", "ACGTACGTACGT")
            >>> print(seq.to_str(max_columns=5))
            >test1
            ACGTA
            CGTAC
            GT
        """
        lines = [f'>{self.header}']
        for i in range(0, len(self.sequence), max_columns):
            lines.append(self.sequence[i:i+max_columns])
        return '\n'.join(lines)


def read_mfasta(path : str) -> List[FastaSequence]:
    """
    A multi-FASTA file contains one or more sequences, where each sequence
    starts with a header line beginning with '>' followed by one or more
    lines of sequence data.
    """
    sequences = []
    with open(path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    sequences.append(FastaSequence(header, ''.join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            sequences.append(FastaSequence(header, ''.join(seq_lines)))
    return sequences