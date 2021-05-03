import numpy as np

seq1 = ["atcg","atcg"]
seq2 = ["atcg","atce"]
global score_gap
global score_match
global score_miss
score_gap,score_match,score_miss = -1,1,-1


def score_matrix(seq1,seq2):#seq1 = y, seq2 = x, [[seq1,seq2],[seq1,seq2]]
    matrix = np.zeros((len(seq1[0])+1,len(seq2[0])+1))
    for y in range(np.shape(matrix)[0]): #shape(y,x)
        for x in range(np.shape(matrix)[1]):
            if y == 0 and x == 0:
                matrix[0][0] == 0
            elif y==0 and x > 0:
                matrix[0][x] = matrix[0][x-1] + score_gap
            elif y>0 and x == 0:
                matrix[y][0] = matrix[y-1][0] + score_gap
            else:
                gap = max(matrix[y][x-1]+score_gap,matrix[y-1][x]+ score_gap)
                symbol = [a[y-1] for a in seq1] + [a[x-1] for a in seq2]
                match_miss = matrix[y-1][x-1] + score_match if all([True if s == seq1[0][y-1] and s != "-" else False for s in symbol]) else matrix[y-1][x-1] + score_miss
                score = max(gap,match_miss)
                matrix[y][x] = score
    return matrix#np matrix

print("global aligment: " ,"\n",score_matrix(seq1,seq2))
#-----------------------------------------------------------------------------------------------------------------

def traceback(matrix):#input matrix output: pos der symbole im aligment zueinander
    seq_code = (len(matrix), len(matrix[-1]))
    #precursor - was davor
    links = (len(matrix),len(matrix[-1])-1) if len(matrix[0])>1 else None
    diagonal = (len(matrix)-1,len(matrix[-2])-1) if len(matrix[0]) > 1 and len(matrix) > 1 else None
    unten = (len(matrix)-1,len(matrix[-2])) if len(matrix) > 1 else None
    precursor = [links,diagonal,unten]
    #end_bedingung
    if any(precursor) == False:#alle sind None
        return [seq_code] #+ [(0,0)]
    #max score und valide moves  davon
    start = matrix[-1][-1]  # ende matrix
    valid_moves = [start-score_gap,start-score_match,start-score_miss] #gap, match, miss
    precursor = [x for x in precursor if matrix[x[0]-1][x[1]-1] in valid_moves]
    precursor = [x for x in precursor if matrix[x[0]-1][x[1]-1] == max([matrix[x[0]-1][x[1]-1] for x in precursor])]
    precursor = precursor[0] #(y,x)
    #rekursion
    new_matrix = [t[:precursor[1]] for t in matrix[:precursor[0]] ]
    return [seq_code] + traceback(new_matrix)

print("seq: ", traceback(score_matrix(seq1, seq2)))

#------------------------------------------------------------------------------------------------------------------
#def seq_aligment():#input [(y,x),(y,x),(y,x)]
