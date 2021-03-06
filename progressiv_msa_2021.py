import numpy as np
import copy
import random
import os as os
import sys as sys

#safe location, variablen input anpassen!
curr_path = os.getcwd()
in_script = sys.argv[0]
#global pfad input
global input_file
input_file = sys.argv[1]
input_file_pfad = curr_path + "/" + input_file # python3 [msa.py,file.sto]

#read fasta/sto
with open(input_file_pfad, "r") as datei:
    test = datei.readlines()
    #print("end: ", test[-1])
    test = [x for x in test if x != ""]
    test_anno = [x.replace("\n","") for x in test if x.startswith("#")]
    test_namen=[x.replace("\n","").split()[0] for x in test if (x.startswith("#")==False) and (x.startswith("/")==False)]
    test = [x.replace("\n","").replace(".","").split()[1] for x in test if x.startswith("#")==False and x.startswith("/")==False]
    test_namen = dict(zip(test,test_namen))
    datei.close()

#head = annotation,
global head
head = [input_file,input_file_pfad,test_anno,test_namen,test]

#
check_path = curr_path+"/results/" + "msa_check_" + input_file
check = open(check_path,"w")

[check.write(f"#={x}\n") for x in head]
#gap cost in bh62 matrix
global gap_c
gap_c = -4
global var_gap_c
var_gap_c = -1

#redumdante aligments
global align_dic
align_dic = {}

# * gap symbol durch - ersetzt -4 durch var_gap cost erseten
b62h={a: i for i,a in enumerate(' A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  - '.split())}
b62 = np.array([
[ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
[-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
[-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
[-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
[ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
[-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
[-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
[ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
[-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
[-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
[-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
[-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
[-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
[-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
[-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
[ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
[ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
[-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
[-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
[ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
[-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
[-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
[ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
[-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1]])
#------------------------------------------------------------------------------------------
#opt globales aligment
def path_finder_function(align1, align2): #["str1","str2"],["str1","str2"]
    #redumdante aligments
    if str(align1)+str(align2) in align_dic:
        return align_dic[str(align1)+str(align2)]

    #gap und multiple gap cost in eine funktion zusammenbringen
    #erstellt vektor mit alles var_gap und findet opt
    def gap_cost_function(score_matrix):
        x_axe = [value + gap_c + var_gap_c * (y - i - 1) for i, value in enumerate(score_matrix[x][:y])]
        y_axe = [value + gap_c + var_gap_c * (x - i - 1) for i, value in enumerate(np.take(score_matrix, y, axis=1)[:x])]
        if max(x_axe) >= max(y_axe):
            score = [max(x_axe),"x_axe", np.argmax(x_axe)+1]
        elif max(x_axe) < max(y_axe):
            score = [max(y_axe),"y_axe", np.argmax(y_axe)+1]
        return score#(score,weg,leange)
    #score=durschnitt aller comb
    def match_cost_function(align1,align2,sub_matrix=b62,sub_head=b62h):
        x_basen,y_basen = [seq[x - 1] for seq in align1],[seq[y - 1] for seq in align2]
        score = sum([sub_matrix[sub_head[x.upper()],sub_head[y.upper()]] if y.upper() in sub_head and x.upper() in sub_head else gap_c for x in x_basen for y in y_basen]) / (len(x_basen) * len(y_basen))
        return score
    #matrix f??r aligment pfad score
    pfad_matrix = np.zeros((len(align1[0])+1,len(align2[0])+1))#zeros((shape),dtype)
    #merkt sich ersten, min pfad durch die matrix, [x,y,[gap oder match,leange gap fuer variable gap cost]]
    move_matrix = np.empty((len(align1[0])+1,len(align2[0])+1,2),dtype=object)
    for x in range(np.shape(pfad_matrix)[0]):
        #edge_case
        for y in range(np.shape(pfad_matrix)[1]):
            if x == 0 and y == 0:
                pfad_matrix[0][0] = 0
                move_matrix[0][0] = ["match",0]
            elif x==0:
                pfad_matrix[x][y] = gap_c + var_gap_c * (y - 1)
                move_matrix[0][y] = ["x_axe",1]
            elif y ==0:
                pfad_matrix[x][y] = gap_c + var_gap_c * (x - 1)
                move_matrix[x][0] = ["y_axe",1]
            else:
                #find best score
                score_match = pfad_matrix[x-1][y-1] + match_cost_function(align1,align2)
                gap = gap_cost_function(pfad_matrix)
                #legalmove_matrix, schaltstelle_pref = gap > match
                if score_match > gap[0]:
                    move_matrix[x][y] = ["match", 0]
                    pfad_matrix[x][y] = score_match
                elif gap[0] >= score_match:
                    pfad_matrix[x][y] = gap[0]
                    move_matrix[x][y] = gap[1:]
    #matrix wird in jedem schritt verkleinert bis edge
    def traceback_function(move_matrix):  # input=np.array 3d

        position = [] #[(0,1),(1,1),(2,2)...]
        while (np.shape(move_matrix)[0] > 0 and np.shape(move_matrix)[1] > 0):
            if move_matrix[-1][-1][0] == "match":
                position += [np.shape(move_matrix)[:2]]
                move_matrix = np.delete(move_matrix, (-1), axis=0)
                move_matrix = np.delete(move_matrix, (-1), axis=1)
            elif move_matrix[-1][-1][0] == "x_axe":#links,spalte
                for _ in range(len(move_matrix[-1]) - move_matrix[-1][-1][1]):
                    position += [np.shape(move_matrix)[:2]]
                    move_matrix = np.delete(move_matrix, (-1), axis=1)
            elif move_matrix[-1][-1][0] == "y_axe":#zeile,oben
                for _ in range(len(move_matrix) - move_matrix[-1][-1][1]):
                    position += [np.shape(move_matrix)[:2]]
                    move_matrix = np.delete(move_matrix, (-1), axis=0)
        return position  # [[()],[(),()]]
    #return score_matrix #2d numpy array
    position = traceback_function(move_matrix)
    #redumdanten align
    align_dic[str(align1)+str(align2)] = position
    #check
    check.write(f"#align1: {align1}, align2: {align2} : {position}\n")
    return position#[(),()]

#min score, luecken sind zaelen als +1
#uneanlichkeitsscore
def score_function(seq1,seq2):
    score = sum([0 if x == y and x !="-" else 1 for x,y in zip(seq1,seq2)])
    check.write(f"#seq1,seq2: {seq1} : {seq2}, score: {score}\n")
    return score
#--------------------------------------------------------------------------------------------------------------
#macht aus pos seq mit luecen
def sequence_function(seq1,seq2,pfad):#-2 weil shape verwendet und index #[[""],[""],[""]]
    #print("start: seq1,seq2,aligment: " , seq1," ",seq2)
    #pfad in enum ist 1 k??rzer am anfang als der der die last pos trackt, -2 bei pos da shape bei (1,1) anf??ngt
    new_seq1 = ["".join(["-" if pfad[::-1][i][0] == pos[0] else s1[pos[0]-2] for i,pos in enumerate(pfad[-2::-1])]) for s1 in seq1]
    new_seq2 = ["".join(["-" if pfad[::-1][i][1] == pos[1] else s2[pos[1]-2] for i,pos in enumerate(pfad[-2::-1])]) for s2 in seq2]
    check.write(f"#seq1,seq2: {seq1} {seq2}, pfad: {pfad}, aligment: {new_seq1} {new_seq2}\n")
    return new_seq1 + new_seq2 #["","",""]
#-------------------------------------------------------------------------------------------------------------------------
#erzeugt guide tree
def guide_tree_function(score_matrix):
    score_matrix = score_matrix.astype(float)
    guide_t = []
    ind = [i for i, seq in enumerate(score_matrix)]
    while(np.shape(score_matrix)!=(1,1)):
        min_d = sorted([(i,j) for i,y in enumerate(score_matrix) for j,x in enumerate(y[:i]) if x == min([x for i,y in enumerate(score_matrix) for x in y[:i]]) and x != 0][0])
        guide_t += [[ind[min_d[0]],ind[min_d[1]]]]
        del ind[min_d[1]]
        score_matrix = [[y if j != min_d[0] and i != min_d[0] else (y+score_matrix[min_d[1]][j])/2 if j != min_d[0] else (y + score_matrix[i][min_d[1]])/2 for j,y in enumerate(x)] for i,x in enumerate(score_matrix)]
        score_matrix = np.delete(score_matrix,min_d[1],axis=0)
        score_matrix = np.delete(score_matrix,min_d[1],axis=1)
        score_matrix[min_d[0]][min_d[0]] = 0
    return guide_t #output [(paar1),[(),[()]]], seq_liste die mit verkleinert wird?
#---------------------------------------------------------------------------------------------------------------
#findet fuer seq die mit min uneanlichkeit und zu dieser die uneanlichste seq
def close_bransh_finder(score_matrix,offset=0):
    score_matrix = score_matrix.astype(float)
    star = np.argmin([sum(x) for x in score_matrix[:len(score_matrix)-offset]])
    guide = [star]
    guide += [np.argmin(score_matrix[star][len(score_matrix)-(len(score_matrix)+offset):])]
    return guide

#macht aus liste von seq eine quad. matrix mit uneanlichkeits scores
def score_matrix_function(seq_list):
    pos_list = [(s1,s2,path_finder_function(s1,s2)) for i,s1 in enumerate(seq_list) for j,s2 in enumerate(seq_list[:i+1])]
    seq_list = [sequence_function(*args) for args in pos_list]
    score_list = [score_function(x,y) if x!=y else 0 for x,y in seq_list]
    score_list = [[x for x in score_list[j+1:i+1]] for j,i in zip([-1] + [i for i,j in enumerate(score_list) if j == 0],[i for i,j in enumerate(score_list) if j == 0])]
    score_list = [x + [y[i] for y in score_list[i+1:]] for i,x in enumerate(score_list)]
    score_list = np.array(score_list)
    return score_list
#----------------------------------------------------------------------------------------------------------------
#findet richtige seq fuer output
def annotion_find(seq,anotations=test_namen):
    return anotations[seq.replace("-","").replace(".","")]

#bald main
def msa(seq_list,head=head):#["","","",]
    seq_list = [[x] for x in seq_list] #[[""],[""],[""]
    score_matrix = score_matrix_function(seq_list)
    guide_tree = guide_tree_function(score_matrix)

    head+= [f"seq_list: {seq_list}"]
    #fehlerhaft
    #head+= [f"score_matrix: {score_matrix}"]
    #[result.write(f"#{x}\n") for x in score_matrix]
    head+=[f"guide_tree: {guide_tree}"]

    [check.write(f"#={x}\n") for x in head]

    #msa baum methode----------------------------------------------------------
    #seq werden anhand guide tree aligniert, bei aligments mehrerer seq wird durschnitt der jeweiligen pos als weert genommen
    tree_list = copy.deepcopy(seq_list)

    tree_msa = open(os.getcwd() + "/results/tree_msa_"+ input_file,"w")
    tree_msa.write("# STOCKHOLM 1.0\n")
    #[tree_msa.write(f"#={x}\n") for x in head]

    for pos in guide_tree:
        matrix = path_finder_function(tree_list[pos[0]], tree_list[pos[1]])
        seq = sequence_function(tree_list[pos[0]],tree_list[pos[1]],matrix)
        tree_list[pos[0]] = seq
    #tree_msa.write("#=tree_aligment\n#=laengen_check: "+str(set([len(x) for x in tree_list[0]]))+"\n#=result:\n")
    [tree_msa.write(f"{annotion_find(w)}\t{w}\n") for w in tree_list[0]]
    tree_msa.write("//\n")
    tree_msa.close()
    [check.write(f"{w}\n") for w in tree_list[0]]
    check.write("//\n")
    #tree_star----------------------------------------------------
    #guide tree -> aligmentreinfolge
    #bei mehr als 2 seq gegeneinander wird ermittelt welche seq den geringsten average uneanlichkeit hat und diese mit
    #der aenlichsten seq aus alig2 aligniert
    #luecken werden nach vorblid des musteraligment in andere aligmen eingefuegt
    tree_star_msa = open(os.getcwd() + "/results/tree_star_msa_"+ input_file,"w")
    tree_star_msa.write("# STOCKHOLM 1.0\n")
    #[tree_star_msa.write(f"#={x}\n") for x in head]

    star_list = copy.deepcopy(seq_list)
    for pos in guide_tree:
        seq = [[x] for x in star_list[pos[0]] + star_list[pos[1]]]
        score_matrix = score_matrix_function(seq)
        star_guide = close_bransh_finder(score_matrix,len(star_list[pos[1]]))
        matrix = path_finder_function([star_list[pos[0]][star_guide[0]]], [star_list[pos[1]][star_guide[1]]])
        seq = sequence_function(star_list[pos[0]],star_list[pos[1]],matrix)
        star_list[pos[0]] = seq
    #tree_star_msa.write("#=ergebniss_star-tree_aligment\n#=laengen_check: "+str(set([len(x) for x in star_list[0]]))+"\n#=result:\n")
    [tree_star_msa.write(f"{annotion_find(w)}\t{w}\n") for w in star_list[0]]
    tree_star_msa.write("//\n")
    tree_star_msa.close()
    [check.write(f"{w}\n") for w in star_list[0]]
    check.write("//\n")

    return None
"""
to do :
    -in blosum62 die gap costen variable machen
    refacctor code -> class, executable 
    output zwischenschritte in file 
    track time
    np.where(matrix==matrix.max/min()) -> bessere variante coor von val zu finden
"""
check.close()
#print(align_dic.keys())