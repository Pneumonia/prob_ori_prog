pfad = "/home/lary/Documents/Studium/4.Semester/problemorientierte_programmierung/beleg2021/PF00545_seed_130_seq.fasta"
tree = "/home/lary/workwork/PycharmProjects/test00/brojekt_biot/tree_result.txt"
with open(pfad,"r") as datei:
    original = datei.readlines()
    datei.close()

ano = [x for x in original if "#" in x]
seq = [x.replace(".","").replace("\n","") for x in original if "#" not in x][:-2]
print("ano: ", ano)
print("seq: ",len(seq) ,seq)
seq_name = [x.split()[0] for x in seq]
seq = [x.split()[1] for x in seq]

print("seq_name: ", seq_name,seq)

with open(tree,"r") as datei:
    tree_msa = datei.readlines()
    datei.close()

tree_msa = [x.replace("-","") + " " + x for x in tree_msa]

print("tree_msa: ", len(tree_msa),tree_msa)
for i,x in enumerate(tree_msa):
    s = x.split()[0]
    ind = seq.index(s)
    tree_msa[i] = seq_name[ind] + " " + x.split()[1].replace("\n","")

print(len(tree_msa))
print(tree_msa)
print(set([len(x.split()[1]) for x in tree_msa]))


datei = open(tree,"w")

for x in ano[:-1]:
    datei.write(x)

for x in tree_msa:
    datei.write(x+ "\n")
datei.write("\\")
datei.close()
