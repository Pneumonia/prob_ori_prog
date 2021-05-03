#Import

#lokal aligment

def swiss_waterman(seq1,seq2): #input: ("str","str")
    gap_cost, miss_match,match = -1,-1,1
    trace_matrix = [[0]*(len(seq2)+1)] #[  [inhalt zeile,x],[] zeilen,y ] #erste zeile nur 0, deswegen y+1
    for y in range(len(seq1)):#zeilen
        trace_matrix += [[0]]#erster eintrag jeder  zeile = 0,(erste spalte = 0),deshalb x +1
        for x in range(len(seq2)):#spalte  in zeile
            maximum =max(trace_matrix[y+1][x]+gap_cost,trace_matrix[y][x+1]+gap_cost,trace_matrix[y][x]+match if seq1[y]==seq2[x]else trace_matrix[y][x]+miss_match)
            trace_matrix[y+1] = trace_matrix[y+1] + [0 if maximum <= 0 else maximum]
    lokal_aligment = [(i,j) for i,y in enumerate(trace_matrix) for j,x in enumerate(y) if x == max(x for y in trace_matrix for x in y)]#coordinaten
    def tracer(trace_matrix,lokal_aligment):#start bei max -> ende 0#überarbeiten
        y_c, x_c = lokal_aligment[0][0], lokal_aligment[0][1]
        scoreXY = [trace_matrix[y_c][x_c] - gap_cost, trace_matrix[y_c][x_c] - miss_match,trace_matrix[y_c][x_c] - match]  # rückverfolgend
        coordinate_score = [trace_matrix[y_c][x_c-1],trace_matrix[y_c-1][x_c-1],trace_matrix[y_c-1][x_c]]
        if any(coordinate_score)==0:#why not all? konterintuitive oder ich bin blöd, nachlesen!
            return lokal_aligment
        coordinate = list(zip(coordinate_score,[(y_c,x_c-1),(y_c-1,x_c-1),(y_c-1,x_c)]))
        coordinate = [x[1] for x in coordinate if x[0] == max(x[0] for x in coordinate if x[0] in scoreXY)]#legaler move
        for x in coordinate:
            lokal_aligment = [x] + lokal_aligment # muss ans ende damit [0][0]
            return tracer(trace_matrix,lokal_aligment)#ka aber funktioniert, glaub ich zumindest

    lokal_aligment = [tracer(trace_matrix,[align])[1:] for align in lokal_aligment]
    for seq in lokal_aligment:
        gap_check1, gap_check2,seq_safe1,seq_safe2 = None, None,[],[]
        for base in seq:
            seq_safe1 += [seq1[base[0]-1] if gap_check1 != base[0] else "-"]
            seq_safe2 += [seq2[base[1]-1] if gap_check2 != base[1] else "-"]
            gap_check1,gap_check2 = base[0],base[1]#lücken erkennen 
        lokal_aligment[lokal_aligment.index(seq)] = list(zip(seq_safe1,seq_safe2))
    return lokal_aligment


#test_lokal = swiss_waterman("atcgatcg","gggeee")
#print("lokal: ", test_lokal)
####################################################################################################################################
#global aligment

def needlemann_wunsch(seq1,seq2): #input: ("[str1_1,str1_2]","[str2,]")
    gap_cost, miss_match,match = -1,-1,1 #gap und missmatch müssen 0>sein
    trace_matrix = [list(range(0,(gap_cost*len(seq2[0]))-1,gap_cost))] #pro spalte  um gapcosten verringern
    for y in range(len(seq1[0])): #zeilen == len str1_1 == str1_2....
        trace_matrix += [[trace_matrix[y][0] +gap_cost]]#pro zeile
        for x in range(len(seq2[0])):#spalte  in zeile len str2_1 == len str2_2....
            gap_seq1 = trace_matrix[y+1][x]+gap_cost#same y
            gap_seq2 = trace_matrix[y][x+1]+gap_cost#same x
            base_y_x = [[s1[y] for s1 in seq1]] + [[s2[x] for s2 in seq2]] #[[,,],[,,]]
            score_diagonal = sum([trace_matrix[y][x] + match if s1 == s2 and s1 != "-" else trace_matrix[y][x]+miss_match for s1 in base_y_x[0] for s2 in base_y_x[1]])/(len(seq1)*len(seq2))
            maximum = max(gap_seq1,gap_seq2,score_diagonal)
            trace_matrix[y+1] = trace_matrix[y+1] + [maximum]
    print(trace_matrix)
    def tracer(trace_matrix):#start bei max trace punkt -> ende beginn#überarbeiten#brauche maß für legale schritte..
        global_aligment = [(len(trace_matrix)-1,len(trace_matrix[0])-1)]
        y_c, x_c = global_aligment[0][0], global_aligment[0][1]
        if y_c<= 0 and x_c <=0:
            return [(0,0)]
        elif x_c <=0:
            coordinate = list(zip([trace_matrix[y_c-1][x_c]],[(y_c-1,x_c)]))
        elif y_c <=0:
            coordinate = list(zip([trace_matrix[y_c][x_c-1]],[(y_c,x_c-1)]))
        else:
            coordinate = list(zip([trace_matrix[y_c][x_c-1],trace_matrix[y_c-1][x_c-1],trace_matrix[y_c-1][x_c]],[(y_c,x_c-1),(y_c-1,x_c-1),(y_c-1,x_c)]))
        coordinate = [x[1] for x in coordinate if x[0] == max(x[0] for x in coordinate)]
        print("coor: ",coordinate)
        for x in coordinate:
            print(x)

            new_trace=[r[:x[1]+1] for r in trace_matrix[:x[0]+1]]

            #return min([tracer(new_trace) + [x]],key=len)  #score  additing?
    global_aligment = tracer(trace_matrix)[1:] + [(len(seq1[0]),len(seq2[0]))]

    seq_safe =  [[seq[pos[0]-1] if pos[0] != global_aligment[i-1][0] else "-" for i,pos in enumerate(global_aligment)]for seq in seq1]
    seq_safe += [[seq[pos[1]-1] if pos[1] != global_aligment[i-1][1] else "-" for i,pos in enumerate(global_aligment)]for seq in seq2]
    seq_safe = ["".join(_[1:]) for _ in seq_safe]

    global_aligment = seq_safe

    return global_aligment #["seq1","seq2"]

#test_seq = ["aaaeee","eeefff","fffggg","gggeee","ggeeff"]

test_global = needlemann_wunsch(["atcg"],["aagg"])
print("global: ",test_global)
###########################################################################################################################################
#distance  Matrix
def distance_matrix(seq_array):#["seq1","seq2","seq3"] #euclidean  von julius besser
    def distance(seq1,seq2):#ohne sub_matrix-> einfügen
        seq_d = 0
        global_seq = needlemann_wunsch([seq1],[seq2])
        seq_d =  sum([0 if global_seq[0][i] == global_seq[1][i] else 1 for i,y in enumerate(global_seq[0])])
        return seq_d
    d_matrix = [[0 if y == x else distance(y,x) for x in seq_array] for y in seq_array]
    return d_matrix

#test_distance = (distance_matrix(["atcg","atce","atee","aeee","eeee"]))
#for line in test_distance:
#    print(*line)

####################################################################################################################################
#guide Tree
def guide_tree(d_matrix):#[[],[],[]]
    if len(d_matrix) == len(d_matrix[1]):
        tree_verlauf = [x for x,_ in enumerate(d_matrix)]
    else:
        tree_verlauf = d_matrix[0]
        d_matrix = d_matrix[1:]
    if len(d_matrix) <= 1:
        return []
    min_coor  = [(y,x) for y,_ in enumerate(d_matrix) for x,_ in enumerate(d_matrix) if d_matrix[y][x] == min(t for i,d in enumerate(d_matrix) for j,t in enumerate(d) if j != i)]
    #erstmal nür fürs erste min_coor -> verbessern!!
    tree = sorted(min_coor[0])#(,)-> für beginn #sorted wichtig sonst index probleme
    new_matrix = []
    for y in range(len(d_matrix)):
        new_matrix += [[]]
        for x in range(len(d_matrix[y])):
            if y not in tree and x not in tree:
                new_matrix[y] += [d_matrix[y][x]]
                continue
            elif y in tree:
                new_matrix[y] += [(d_matrix[tree[0]][x] + d_matrix[tree[1]][x]) / 2]
            elif x in tree:
                new_matrix[y] += [(d_matrix[y][tree[0]] + d_matrix[y][tree[1]]) / 2]
    d_matrix = [y[:tree[0]]+y[tree[0]+1:] for y in new_matrix[:tree[0]]+ new_matrix[tree[0]+1:]]
    d_matrix = [[d_matrix[y][x] if x != y else 0 for y,_ in enumerate(d_matrix) ]for x,_ in enumerate(d_matrix)]

    tree_verlauf[tree[0]] = [tree_verlauf[tree[0]]] + [tree_verlauf[tree[1]]]
    tree_verlauf[tree[1]:tree[1] + 1] = []


    d_matrix = [tree_verlauf] + d_matrix
    return [tree_verlauf[tree[0]]] + guide_tree(d_matrix)

#print("guide: ",guide_tree(test_distance))
##################################################################################################################################
#msa
test_seq = ["aaaeeeaaa","eeefffaaa","fffgggaaa","gggeeeaaa","ggeefffaa"]
def msa(seq_array):#["seq1","","","seqn"]
    distanceMatrix = distance_matrix(seq_array)
    guideTree = guide_tree(distanceMatrix)
    print(distanceMatrix)
    seq_array = [[x]for x in seq_array]

    named_array = [str(x) for x in range(len(seq_array))]

    for x in guideTree:

        seq1 =seq_array[named_array.index(str(x[0]))]
        seq2 =seq_array[named_array.index(str(x[1]))]
        seq_array[named_array.index(str(x[0]))] = needlemann_wunsch(seq1,seq2)
        del(seq_array[named_array.index(str(x[1]))])
        named_array[named_array.index(str(x[0]))] = str(x)
        del (named_array[named_array.index(str(x[1]))])

    return seq_array

#test = msa(test_seq)
#print(test)
#class
