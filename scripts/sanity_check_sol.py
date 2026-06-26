objective=-1
with open("test.sol") as f:
    vars = dict()
    for line in f:
        
        line=line.strip()
        if line.startswith("#"):
            objective=int(line.split(" = ")[1])
            continue
        var,val = line.split()
        vars[var]=int(val)


fin_kmers = set()
for var,vale in vars.items():

    if var.startswith("R") and vale==1:
        fin_kmers.add(var.split("_")[1])


degens = dict()
for kmer in fin_kmers:
    degens[kmer]=[[] for _ in kmer]
for var, vale in  vars.items():
    if var.startswith("t") and vale==1:
        _,kmer,i,a = var.split("_")
        i=int(i)
        if kmer in fin_kmers:
            degens[kmer][i].append(a)

print(degens)
    
