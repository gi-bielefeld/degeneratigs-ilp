from typing import Set
import sys
def get_kmerset(s : str, k : int) -> Set[str]:
    return set(s[i:i+k] for i in range(len(s)-k+1))


def build_trie(kmerset : Set[str]):
    root = dict()
    for kmer in kmerset:
        curr_level = root
        for char in kmer:
            if not char in curr_level:
                curr_level[char]=dict()
            curr_level=curr_level[char]
    return root

def get_forbidden_prefixes(trie,alphabet,k):
    stack = [("",trie)]
    forbidden = set()
    while len(stack):
        prefix,subtrie = stack.pop()
        for a in alphabet:
            if a in subtrie:
                stack.append((prefix+a,subtrie[a]))
            else:
                if len(prefix+a) > k:
                    continue 
                forbidden.add(prefix+a)
    return forbidden

def get_hamming_neighbors(kmer,alphabet):
    neighbors = []
    for i, c1 in enumerate(kmer):
        for c2 in alphabet:
            if c2!=c1:
                neighbors.append(kmer[0:i]+c2+kmer[i+1:])
    return neighbors

def get_forbidden_neighbors_dict(kmerset,alphabet):
    res = dict()
    for kmer in kmerset:
        neighbors = get_hamming_neighbors(kmer,alphabet)
        neighbors.difference_update(kmerset)
        res[kmer]=neighbors
    return res 

def get_neighborless(kmerset,alphabet):
    neighborless = set()
    for kmer in kmerset:
        real_neighbors = [neighbor for neighbor in get_hamming_neighbors(kmer,alphabet) if neighbor in kmerset]
        if len(real_neighbors)==0:
            neighborless.add(kmer)
    return neighborless

def get_breakmers(kmerset,alphabet):
    breakmers = set()
    print(f"Determining neighborless k-mers...",file=sys.stderr)
    neighborless=get_neighborless(kmerset,alphabet)
    print(f"Of {len(kmerset)} k-mers {len(neighborless)} are neighborless.",file=sys.stderr)
    print(f"Determining break k-mers...",file=sys.stderr)
    for kmer in neighborless:
        next = [kmer[1:]+a for a in alphabet if kmer[1:]+a in kmerset]
        prev = [a+kmer[:-1] for a in alphabet if a+kmer[:-1] in kmerset]
        if len(next)<=1 and len(prev)<=1:
            nblnext = len(next)==0 or next[0] in neighborless
            nblprev = len(prev)==0 or prev[0] in neighborless
            if nblnext and nblprev:
                breakmers.add(kmer)
    print(f"Of {len(kmerset)} k-mers {len(breakmers)} are breakmers.",file=sys.stderr)
    return breakmers


def get_components(kmerset : Set[str],alphabet):
    components = []
    while len(kmerset) > 0:
        stack = [next(iter(kmerset))]
        component = set()
        while len(stack)>0:
            curr = stack.pop()
            if not curr in kmerset:
                continue
            kmerset.discard(curr)
            component.add(curr)
            for a in alphabet:
                nxt = curr[1:]+a
                prev = a+curr[:-1]
                if nxt in kmerset:
                    stack.append(nxt)
                if prev in kmerset:
                    stack.append(prev)
            for neighbor in get_hamming_neighbors(curr,alphabet):
                if neighbor in kmerset:
                    stack.append(neighbor)
        if len(component)>0:
            components.append(component)
    return components
                
        