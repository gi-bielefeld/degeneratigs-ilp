import util as ut
from fasta import read_mfasta
import sys
from argparse import ArgumentParser
import os
import gurobipy as gp
from collections import Counter
import networkx as nx

def ilp(kmerset,k,alphabet,baselength=float("inf")):
    binaries = []
    generals = []
    constraints = []
    print("Using base string length",baselength,file=sys.stderr)
    baselength = min((k+1)*len(kmerset),baselength)
    obj_terms = []
    trie = ut.build_trie(kmerset)
    forbidden = ut.get_forbidden_prefixes(trie,alphabet,k)
    for i in range(baselength):
        charvars = []
        for c in alphabet:
            charvars.append(f"t_{i}_{c}")
            constraints.append(f"t_{i}_{c} + t_{i}_stop <= 1")
        binaries.extend(charvars)
        binaries.append(f"t_{i}_stop")
        obj_terms.append(f"t_{i}_stop")
    for pat in forbidden:
        for i in range(baselength-len(pat)+1):
            forbidden_terms = []
            for j,c in enumerate(pat):
                forbidden_terms.append(f"t_{i+j}_{c}")
            not_all = len(forbidden_terms)-1
            constraints.append(" + ".join(forbidden_terms) + f" <= {not_all}")
    for kmer in kmerset:
        kmerterms = []
        for i in range(baselength-k+1):
            kmerterms.append(f"r_{kmer}_{i}")
            for j,c in enumerate(kmer):
                constraints.append(f"r_{kmer}_{i} - t_{i+j}_{c} <= 0")
        binaries.extend(kmerterms)
        constraints.append(" + ".join(kmerterms)+" >= 1")
    return"maximize " + " + ".join(obj_terms), constraints, binaries, generals

def ilp_graph(kmerset,k,alphabet):
    binaries = []
    generals = []
    constraints = defaultdict()
    obj_terms = []
    max_l = len(kmerset)
    for kmer in kmerset:
        for i, c in enumerate(kmer):
            for a in alphabet:
                binaries.append(f"t_{kmer}_{i}_{a}")
            constraints.append()

    return obj, constraints,binaries,generals,


def get_variables_graph(kmerset,k,alphabet):
    binaries = []
    generals = []
    for kmer in kmerset:
        binaries.append(f"R_{kmer}")
        binaries.append(f"s_{kmer}")
        binaries.append(f"e_{kmer}")
        generals.append(f"l_{kmer}")
        for i, c in enumerate(kmer):
            for a in alphabet:
                binaries.append(f"t_{kmer}_{i}_{a}")
                
                if a != c:
                    neighbor = kmer[0:i]+a+kmer[i+1:]
                    if a < c and neighbor in kmerset:
                        binaries.append(f"r_{kmer}_{neighbor}")
        #next
        for a in alphabet:
            next = kmer[1:]+a
            if next in kmerset and next!=kmer:
                binaries.append(f"x_{kmer}_{next}")
    return binaries, generals


def add_variables_gurobi(kmerset,k,alphabet, model : gp.Model):
    variables = dict((x,dict()) for x in "Rseltrx")
    for kmer in kmerset:
        variables["R"][kmer]=model.addVar(name=f"R_{kmer}",vtype=gp.GRB.BINARY)
        variables["s"][kmer]=model.addVar(name=f"s_{kmer}",vtype=gp.GRB.BINARY)
        variables["e"][kmer]=model.addVar(name=f"e_{kmer}",vtype=gp.GRB.BINARY)
        variables["l"][kmer]=model.addVar(name=f"l_{kmer}",vtype=gp.GRB.INTEGER)
        for i, c in enumerate(kmer):
            for a in alphabet:
                variables["t"][(kmer,i,a)]=model.addVar(name=f"t_{kmer}_{i}_{a}",vtype=gp.GRB.BINARY)
                
                if a != c:
                    neighbor = kmer[0:i]+a+kmer[i+1:]
                    if a < c and neighbor in kmerset:
                        variables["r"][(kmer,neighbor)]=model.addVar(name=f"r_{kmer}_{neighbor}",vtype=gp.GRB.BINARY)
        #next
        for a in alphabet:
            next = kmer[1:]+a
            if next in kmerset and next!=kmer:
                variables["x"][(kmer,next)]=model.addVar(name=f"x_{kmer}_{next}",vtype=gp.GRB.BINARY)
    model.update()
    return variables
                        

def print_g01(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        for i,c in enumerate(kmer):
            print(f"c01.{n}: t_{kmer}_{i}_{c} = 1",file=file)
            n+=1


def add_g01(kmerset,alphabet,variables,model):
    n=1
    for kmer in kmerset:
        for i,c in enumerate(kmer):
            #print(f"c01.{n}: t_{kmer}_{i}_{c} = 1",file=file)
            t = variables["t"][(kmer,i,c)]
            model.addLConstr(t, gp.GRB.EQUAL, 1,name=f"c01.{n}")
            n+=1

def print_g0203(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        for neighbor in ut.get_hamming_neighbors(kmer,alphabet):
            if neighbor in kmerset and neighbor < kmer:
                for i in range(len(kmer)):
                    for a in alphabet:
                        print(f"c02.{n}: t_{kmer}_{i}_{a} - t_{neighbor}_{i}_{a} + r_{kmer}_{neighbor} <= 1",file=file)
                        print(f"c03.{n}: t_{neighbor}_{i}_{a} - t_{kmer}_{i}_{a} + r_{kmer}_{neighbor} <= 1",file=file)
                        n+=1

def add_g0203(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        for neighbor in ut.get_hamming_neighbors(kmer, alphabet):
            if neighbor in kmerset and neighbor < kmer:
                for i in range(len(kmer)):
                    for a in alphabet:
                        t1 = variables["t"][(kmer, i, a)]
                        t2 = variables["t"][(neighbor, i, a)]
                        r  = variables["r"][(kmer, neighbor)]

                        model.addLConstr(
                            t1 - t2 + r,
                            gp.GRB.LESS_EQUAL,
                            1,
                            name=f"c02.{n}"
                        )
                        model.addLConstr(
                            t2 - t1 + r,
                            gp.GRB.LESS_EQUAL,
                            1,
                            name=f"c03.{n}"
                        )
                        n += 1

def print_g0405(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        for i, c in enumerate(kmer):
            for a in alphabet:
                if a != c:
                    neighbor = kmer[0:i]+a+kmer[i+1:]
                    rvar = f"r_{kmer}_{neighbor}" if kmer > neighbor else f"r_{neighbor}_{kmer}"
                    if neighbor in kmerset:
                        print(f"c0405.{n}: t_{kmer}_{i}_{a} - {rvar} <= 0",file=file)
                        n+=1


def add_g0405(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        for i, c in enumerate(kmer):
            for a in alphabet:
                if a != c:
                    neighbor = kmer[0:i] + a + kmer[i+1:]
                    if neighbor in kmerset:
                        if kmer > neighbor:
                            r = variables["r"][(kmer, neighbor)]
                        else:
                            r = variables["r"][(neighbor, kmer)]

                        t = variables["t"][(kmer, i, a)]

                        model.addLConstr(
                            t - r,
                            gp.GRB.LESS_EQUAL,
                            0,
                            name=f"c0405.{n}"
                        )
                        n += 1

def print_g06(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        for i, c in enumerate(kmer):
            for a in alphabet:
                if a != c:
                    neighbor = kmer[0:i]+a+kmer[i+1:]
                    if not neighbor in kmerset:
                        print(f"c06.{n}: t_{kmer}_{i}_{a} = 0",file=file)
                        n+=1


def add_g06(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        for i, c in enumerate(kmer):
            for a in alphabet:
                if a != c:
                    neighbor = kmer[0:i] + a + kmer[i+1:]
                    if neighbor not in kmerset:
                        t = variables["t"][(kmer, i, a)]
                        model.addLConstr(
                            t,
                            gp.GRB.EQUAL,
                            0,
                            name=f"c06.{n}"
                        )
                        n += 1

def print_g0708(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        for a in alphabet:
            next = kmer[1:]+a
            if next in kmerset:
                for i in range(len(kmer)-1):
                    for c in alphabet:
                        print(f"c07.{n}: t_{kmer}_{i+1}_{c} - t_{next}_{i}_{c} + x_{kmer}_{next} <= 1",file=file)
                        print(f"c08.{n}: t_{next}_{i}_{c} - t_{kmer}_{i+1}_{c} + x_{kmer}_{next} <= 1")
                        n+=1


def add_g0708(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        for a in alphabet:
            nxt = kmer[1:] + a
            if nxt == kmer:
                continue
            if nxt in kmerset:
                for i in range(len(kmer) - 1):
                    for c in alphabet:
                        t1 = variables["t"][(kmer, i+1, c)]
                        t2 = variables["t"][(nxt, i, c)]
                        x  = variables["x"][(kmer, nxt)]

                        model.addLConstr(
                            t1 - t2 + x,
                            gp.GRB.LESS_EQUAL,
                            1,
                            name=f"c07.{n}"
                        )
                        model.addLConstr(
                            t2 - t1 + x,
                            gp.GRB.LESS_EQUAL,
                            1,
                            name=f"c08.{n}"
                        )
                        n += 1

def print_g09(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        terms = []
        for a in alphabet:
            next = kmer[1:]+a
            if next in kmerset:
                terms.append(f"x_{kmer}_{next}")
        xterms = " + ".join(terms)
        print(f"c09.{n}: {xterms} + e_{kmer} - R_{kmer} = 0",file=file)
        n+=1

def add_g09(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        expr = gp.LinExpr()
        for a in alphabet:
            nxt = kmer[1:] + a
            if nxt in kmerset and nxt!=kmer:
                expr += variables["x"][(kmer, nxt)]

        expr += variables["e"][kmer]
        expr -= variables["R"][kmer]

        model.addLConstr(expr, gp.GRB.EQUAL, 0, name=f"c09.{n}")
        n += 1


def print_g10(kmerset,alphabet,file=sys.stdout):
    n=1
    for kmer in kmerset:
        terms = []
        for a in alphabet:
            prev = a+kmer[:-1]
            if prev in kmerset and prev!=kmer:
                terms.append(f"x_{prev}_{kmer}")
        xterms = " + ".join(terms)
        print(f"c10.{n}: {xterms} + s_{kmer} - R_{kmer} = 0",file=file)
        n+=1

def add_g10(kmerset, alphabet, variables, model):
    n = 1
    for kmer in kmerset:
        expr = gp.LinExpr()
        for a in alphabet:
            prev = a + kmer[:-1]
            if prev in kmerset and prev != kmer:
                expr += variables["x"][(prev, kmer)]

        expr += variables["s"][kmer]
        expr -= variables["R"][kmer]

        model.addLConstr(expr, gp.GRB.EQUAL, 0, name=f"c10.{n}")
        n += 1

def print_g11(kmerset,alphabet,file=sys.stdout):
    n=1
    maxlen = len(kmerset)
    for kmer in kmerset:
        print(f"c11.{n}: l_{kmer} + {maxlen} s_{kmer} <= {maxlen}",file=file)

def add_g11(kmerset, alphabet, variables, model):
    n = 1
    maxlen = len(kmerset)

    for kmer in kmerset:
        model.addLConstr(
            variables["l"][kmer] + maxlen * variables["s"][kmer],
            gp.GRB.LESS_EQUAL,
            maxlen,
            name=f"c11.{n}"
        )
        n += 1

def print_g12(kmerset,alphabet,file=sys.stdout):
    n=1
    escape = len(kmerset)+1
    for kmer in kmerset:
        for a in alphabet:
            next = kmer[1:]+a
            if next in kmerset and next!=kmer:
                print(f"c12.{n}: l_{next} - l_{kmer} - {escape} x_{kmer}_{next} >= {1-escape}",file=file)
                n+=1


def add_g12(kmerset, alphabet, variables, model):
    n = 1
    escape = len(kmerset) + 1

    for kmer in kmerset:
        for a in alphabet:
            nxt = kmer[1:] + a
            if nxt in kmerset and nxt != kmer:
                model.addLConstr(
                    variables["l"][nxt]
                    - variables["l"][kmer]
                    - escape * variables["x"][(kmer, nxt)],
                    gp.GRB.GREATER_EQUAL,
                    1 - escape,
                    name=f"c12.{n}"
                )
                n += 1

def print_g1314(kmerset,alphabet,file=sys.stdout):
    g14c = 1
    g13c = 1
    for kmer in kmerset:
        g13terms = []
        for i,c in enumerate(kmer):
            for a in alphabet:
                if a < c:
                    print(f"c14.{g14c}: t_{kmer}_{i}_{a} + R_{kmer} <= 1",file=file)
                    g14c+=1
                    g13terms.append(f"t_{kmer}_{i}_{a}")
        g13terms = " + ".join(g13terms)
        print(f"c13.{g13c}: {g13terms} +  R_{kmer} >= 1")
        g13c+=1


def add_g1314(kmerset, alphabet, variables, model):
    g14c = 1
    g13c = 1

    for kmer in kmerset:
        expr = gp.LinExpr()

        for i, c in enumerate(kmer):
            for a in alphabet:
                if a < c:
                    t = variables["t"][(kmer, i, a)]
                    R = variables["R"][kmer]

                    model.addLConstr(
                        t + R,
                        gp.GRB.LESS_EQUAL,
                        1,
                        name=f"c14.{g14c}"
                    )
                    g14c += 1
                    expr += t

        expr += variables["R"][kmer]

        model.addLConstr(
            expr,
            gp.GRB.GREATER_EQUAL,
            1,
            name=f"c13.{g13c}"
        )
        g13c += 1

def print_objective_g(kmerset,k,file=sys.stdout):
    rvars = []
    svars = []
    for kmer in kmerset:
        rvars.append(f"R_{kmer}")
        svars.append(f"{k-1} s_{kmer}")
    rvars = " + ".join(rvars)
    svars = " + ".join(svars)
    print(f"minimize {rvars} + {svars}",file=file)

def add_objective_g(kmerset, k, variables, model):
    """
    Implements:
        minimize sum(R_kmer) + (k-1)*sum(s_kmer)
    """

    obj = gp.LinExpr()

    for kmer in kmerset:
        obj += variables["R"][kmer]
        obj += (k - 1) * variables["s"][kmer]

    model.setObjective(obj, gp.GRB.MINIMIZE)
    

def fmt_ilp_gurobi(obj,cons,bin,general,file=sys.stdout):
    print(obj,file=file)
    print("subject to",file=file)
    for c in cons:
        print(c,file=file)
    print("binaries")
    for b in bin:
        print(b,file=file)
    print("generals")
    for g in  general:
        print(g,file=file)
    print("end")

def ez_sol(kmerset):
    vars = []
    for kmer in kmerset:
        vars.append(f"s_{kmer} = 1")
        vars.append(f"R_{kmer} = 1")
        for i,c in enumerate(kmer):
            for a in alphabet:
                if a==c:
                    vars.append(f"t_{kmer}_{i}_{a} = 1")
                else:
                    vars.append(f"t_{kmer}_{i}_{a} = 0")
    return vars

def extract_representatives(variables, model, tol=1e-6):
    edges = set()
    representatives = dict()

    for kmer, var in variables["R"].items():
        val = var.X
        if val > tol:
            representatives[kmer]=[]

    for key, var in variables["x"].items():
        val = var.X
        if val > tol:
            edges.add(key)

    # --- t variables ---
    for (kmer,i,a), var in variables["t"].items():
        if not kmer in representatives:
            continue
        val = var.X
        if val > tol:
            representatives[kmer].add((i,a))
    return representatives,edges


def get_degenerate_strings(variables,alphabet, tol=1e-6):
    deg_strs = []
    for kmer, var in variables["s"].items():
        val = var.X
        if val > 1-tol:
            assert(variables["R"][kmer].X > 1-tol)
            deg_str = kmer
            deg_pos = [(i,a) for i,c in enumerate(kmer) for a in alphabet if a!=c and variables["t"][(kmer,i,a)].X>1-tol]
            last=kmer
            while True:
                nxt = [last[1:]+a for a in alphabet]
                nxt = [n for n in nxt if n in variables["R"] and n!=last]
                nxt = [n for n in nxt if variables["x"][(last,n)].X > 1-tol]
                assert(len(nxt)<=1)
                if len(nxt)==0:
                    break
                nxt = nxt[0]
                deg_str=deg_str+(nxt[-1])
                deg_pos.extend([(len(deg_str)-1,a) for a in alphabet if a!=nxt[-1] and variables["t"][(nxt,len(nxt)-1,a)].X>1-tol])
                last=nxt
            deg_strs.append((deg_str,tuple(sorted(deg_pos))))
    return deg_strs

def add_adapters(deg_strs,k,start_adapters,end_adapters):
    for string, degpos in deg_strs:
        #technically k-1-mer
        start_kmer = string[0:k-1]
        start_degpos = [pos for pos,_ in degpos if pos < k]
        end_kmer = string[len(string)-k+1:]
        end_degpos = [pos for pos,_ in degpos if pos >= len(string)-k]
        if len(start_degpos)==0:
            start_adapters[start_kmer]=(string,degpos)
        if len(end_degpos)==0:
            end_adapters[end_kmer]=(string,degpos)

def breakmer_unitigs(breakmers,alphabet):
    unitigs = []
    remaining = set(breakmers.copy())
    while len(remaining)>0:
        start = remaining.pop()
        unitig = start
        last=start
        #extend as far as possible to the right
        while True:
            nxt = [last[1:]+a for a in alphabet]
            nxt = [n for n in nxt if n in remaining]
            assert(len(nxt)<=1)
            if len(nxt)==0:
                break
            else:
                nxt=nxt[0]
                remaining.remove(nxt)
                unitig=unitig+nxt[-1]
                last=nxt
        #extend as far as possible to the left
        last = start
        while True:
            prev = [a+last[:-1] for a in alphabet]
            prev = [p for p in prev if p in remaining]
            assert(len(prev)<=1)
            if len(prev)==0:
                break
            else:
                prev=prev[0]
                remaining.remove(prev)
                unitig=prev[0]+unitig
                last = prev
        unitigs.append(unitig)
        
    return unitigs


def breakmer_unitigs_graph(breakmers,k):
    graph = nx.DiGraph()
    for breakmer in breakmers:
        start = breakmer[0:k-1]
        end = breakmer[1:]
        if not graph.has_node(start):
            graph.add_node(start)
        if not graph.has_node(end):
            graph.add_node(end)
        assert(not graph.has_edge(start,end))
        graph.add_edge(start,end)
    for node in graph.nodes():
        assert(graph.in_degree(node)<=1)
        assert(graph.out_degree(node)<=1)
    walks = extract_simple_paths_and_cycles(graph)
    unitigs = [w[0][0]+"".join(e[-1] for _,e in w) for w in walks]
    return unitigs

def extract_simple_paths_and_cycles(graph):
    walks = []
    ugraph = nx.Graph(graph)
    for comp in nx.connected_components(ugraph):
        start = [v for v in comp if graph.in_degree(v)==0]
        assert(len(start)<=1)
        if len(start)==0:
            start =next(iter(comp))
        else:
            start = start[0]
        
        curr = start
        edges_seen = set()
        walk = []
        while graph.out_degree(curr) >= 1:
            nxt=next(graph.neighbors(curr))
            e=(curr,nxt)
            if e in edges_seen:
                break
            edges_seen.add(e)
            walk.append(e)
            curr=nxt
        walks.append(walk)
    return walks

def degenerate_spectrum(string,deg_pos,k):
    spectrum = Counter()
    for i in range(len(string)-k+1):
        base_kmer=string[i:i+k]
        #spectrum[base_kmer]+=1
        degen_positions = [j for j,_ in deg_pos if j >= i and j < i+k]
        kmers = set([base_kmer])
        for j in degen_positions:
            new_kmers = []
            alt_chars = [c for (x,c) in deg_pos if x==j]
            rel_pos = j-i
            for kmer in kmers:
                for a in alt_chars:
                    new_kmers.append(kmer[0:rel_pos]+a+kmer[rel_pos+1:])
            kmers.update(new_kmers)
        spectrum.update(kmers)
    return spectrum


def tie_unitigs_to_degen_strings(b_unitigs,deg_strs,start_adapters,end_adapters,k):
    deg_strs = set(deg_strs)
    #print(deg_strs)
    #print(start_adapters)
    assert(set(start_adapters.values()).issubset(deg_strs))
    assert(set(end_adapters.values()).issubset(deg_strs))
    unitig_starts = set()
    unitig_ends = set()

    for unitig in b_unitigs:
        assert(set(start_adapters.values()).issubset(deg_strs))
        assert(set(end_adapters.values()).issubset(deg_strs))
        start = unitig[0:k-1]
        end = unitig[len(unitig)-k+1:]
        #unitigs cannot have the same start or end
        assert(start not in unitig_starts)
        assert(end not in unitig_starts)
        assert(end not in unitig_ends)
        assert(start not in unitig_ends)
        unitig_starts.add(start)
        unitig_ends.add(end)
        #print(start,end)
        start_degenerate = False
        end_degenerate = False
        #assert(start!=end)
        if start in end_adapters:
            sstr,sdegpos = end_adapters[start]
            del end_adapters[start]
            deg_strs.remove((sstr,sdegpos))
            start_degenerate=True
        else:
            sstr,sdegpos = start,()
        
        if end in start_adapters:
            estr,edegpos = start_adapters[end]
            del start_adapters[end]
            if (estr,edegpos) != (sstr,sdegpos):
                deg_strs.remove((estr,edegpos))
                end_degenerate=True
            else:
                estr,edegpos = end, ()
        else:
            estr,edegpos = end, ()
        newstr = sstr+unitig[k-1:]+estr[k-1:]
        assert(newstr.startswith(sstr))
        assert(newstr.endswith(estr))
        assert(unitig in newstr)
        #print("Sane???",sstr,unitig,estr,newstr)
        shift = len(sstr)-len(start)+len(unitig)-len(end)
        newdegpos = tuple(list(sdegpos)+[(j+shift,c) for (j,c) in edegpos])
        #update data structures
        deg_strs.add((newstr,newdegpos))
        newstart=newstr[:k-1]
        if len([i for (i,_) in newdegpos if i < k])==0 and start_degenerate:
            start_adapters[newstart]=(newstr,newdegpos)
        newend = newstr[len(newstr)-k+1:]
        if len([i for (i,_) in newdegpos if i >= len(newstr)-k])==0 and end_degenerate:
            end_adapters[newend]=(newstr,newdegpos)
        
    return deg_strs

UNITIG = "u"
DEGENERATE = "d"
FAKESTART = "f"
FAKEEND = "e"
def degenerate_string_graph(b_unitigs,deg_strs,k):
    graph = nx.DiGraph()
    for unitig in b_unitigs:
        start = unitig[0:k-1]
        end = unitig[len(unitig)-k+1:]
        assert(len(start)==k-1)
        assert(len(end)==k-1)
        if not graph.has_node(start):
            graph.add_node(start)
        if not graph.has_node(end):
            graph.add_node(end)
        graph.add_edge(start,end,type=UNITIG,string=(unitig,[]))
    for string, degpos in deg_strs:
        start = string[0:k-1]
        end = string[len(string)-k+1:]
        assert(len(start)==k-1)
        assert(len(end)==k-1)
        if len([i for i,_ in degpos if i < k-1]) > 0 or (graph.has_node(start) and graph.out_degree(start)>0):
            start = (FAKESTART,(string,degpos))
        if len([i for i,_ in degpos if i >= len(string) -k+1]) > 0 or (graph.has_node(end) and graph.in_degree(end)>0):
            end = (FAKEEND,(string,degpos))
        if not graph.has_node(start):
            graph.add_node(start)
        if not graph.has_node(end):
            graph.add_node(end)
        if graph.in_degree(end)>0 or graph.out_degree(start)>0:
            #do not add if alrea
            continue 
        graph.add_edge(start,end,type=DEGENERATE,string=(string,degpos))
    for x in graph.nodes:
        #print(f"x='{x}'",graph[x])
        assert(graph.in_degree(x)<=1)
        assert(graph.out_degree(x)<=1)
    degeneratigs = []
    walks = extract_simple_paths_and_cycles(graph)
    for walk in walks:
        prev_string = ""
        prev_pos = []
        for (u,v) in walk:
            string,degpos = graph[u][v]["string"]

            if prev_string=="":
                prev_string+=string
                shift = 0
            else:
                shift = len(prev_string)-k+1
                prev_string+=string[k-1:]
            for (j,c) in degpos:
                prev_pos.append((j+shift,c))
        degeneratigs.append((prev_string,prev_pos))
                
    return degeneratigs

    
    

    
    
    




if __name__=="__main__":
    #print(degenerate_spectrum("ACC",[(0,"T"),(1,"T"),(2,"T")],3))
    #sys.exit(0)
    parser = ArgumentParser()
    parser.add_argument("fasta",help="Multiple fasta file")
    parser.add_argument("output",help="Output")
    parser.add_argument("--alphabet",help="Provide an alphabet. Otherwise will be inferred from strings.")
    parser.add_argument("-k",type=int,default=31)
    parser.add_argument("--warm-start")
    parser.add_argument("--work-dir",default="degeneratigs_workdir")
    parser.add_argument("--sanity-check-spectra",action="store_true",help="Slow, but validates the spectrum solution")
    parser.add_argument("--timelim",type=int, help="time limit for solving each instance",default=300)
    parser.add_argument("--show-gurobi-logs",action="store_true")
    args = parser.parse_args()
    f = read_mfasta(args.fasta)
    strings = [entry.sequence for entry in  f] #if set(entry.sequence).issubset(set("ACGT"))]
    if not args.alphabet:
        alphabet = set("".join(strings))
    else:
        alphabet = set(args.alphabet)
    print("Using alphabet",alphabet,file=sys.stderr)
    print("Determining kmer-set using k =",args.k,file=sys.stderr)
    all_kmers = set()
    for s in strings:
        kmers = ut.get_kmerset(s,k=args.k)
        all_kmers.update(kmers)
    
    if args.alphabet:
        print("Filtering non alphabet k-mers")
        all_kmers = set(kmer for kmer in all_kmers if set(kmer).issubset(alphabet))
    breakmers = ut.get_breakmers(all_kmers,alphabet)
    print(f"Removing all breakmers from solution",file=sys.stderr)
    all_kmers.difference_update(breakmers)
    print(f"Calculating breakmer unitigs")
    b_unitigs = breakmer_unitigs_graph(breakmers,args.k)
    print(f"Determined {len(b_unitigs)} breakmer unitigs.")
    print(f"Breaking remaining graph ({len(all_kmers)} kmers) into components",file=sys.stderr)
    components=ut.get_components(all_kmers,alphabet)
    start_adapters = dict()
    end_adapters = dict()
    all_deg_strs = []
    n_comp = len(components)
    print(f"Split the problem into {n_comp} components.")
    for compnum,comp in enumerate(components,start=1):
        print(f"Building a model for component {compnum}/{n_comp} with {len(comp)} k-mers")
        model = gp.Model()
        variables = add_variables_gurobi(comp,args.k,alphabet,model)
        for c in [add_g01,add_g0203,add_g0405,add_g06,add_g0708,add_g09,add_g10,add_g11,add_g12,add_g1314]:
            c(comp,alphabet,variables,model)
        add_objective_g(comp,args.k,variables,model)
        model.update()
        if not args.show_gurobi_logs:
            model.setParam("OutputFlag", 0)
        model.setParam("TimeLimit",args.timelim)
        print("Solving...")
        model.optimize()
        print(f"Solved (gap: {model.MIPGap}) {model.objVal}")

        #for _,vs in variables.items():
        #    for _,v in vs.items():
        #        if v.X >0.99:
        #            print(v) 
        deg_strs = get_degenerate_strings(variables,alphabet)
        all_deg_strs.extend(deg_strs)
        if args.sanity_check_spectra:
            spectrum = Counter()
            for string,degpos in deg_strs:
                part_spectrum = degenerate_spectrum(string,degpos,args.k)
                spectrum.update(part_spectrum)
            for kmer,i in spectrum.items():
                assert(i==1)
            assert(set(spectrum)==comp)
        
        add_adapters(deg_strs,args.k,start_adapters,end_adapters)
        model.close()
    print(f"Before chaining: {len(b_unitigs)} breakmer-unitigs and {len(all_deg_strs)} degenerate strings (total {len(b_unitigs) + len(all_deg_strs)})")
    final_deg_strs = degenerate_string_graph(b_unitigs,all_deg_strs,args.k)
    print(f"Determined final {len(final_deg_strs)} degenerate strings of cumulative length {sum(len(s[0]) for s in final_deg_strs)}")
    with open(args.output,"w") as outf:
        for i,(string,degpos) in enumerate(final_deg_strs):
            dgpsstr = ",".join(str(i)+c for i,c in degpos)
            print(f">D{i}: {dgpsstr}",file=outf)
            print(string,file=outf)

    #print(f"Re-attaching breakmer unitigs...")
    #deg_strs = tie_unitigs_to_degen_strings(b_unitigs,all_deg_strs,start_adapters,end_adapters,args.k)
    
    
        
        #print(comp)
        #sys.exit()

    


    
