
pfad = "/home/lary/workwork/PycharmProjects/test00/brojekt_biot"
with open("/home/lary/workwork/PycharmProjects/test00/brojekt_biot/PF00545_seed_130_seq.fasta","r") as datei:
    sto = datei.readlines()
    datei.close()

sto = [x.replace("\n","") if x.startswith("#")==True else x.replace(".","").replace("\n","") for x in sto]
print(sto)

new = open(pfad+"/PF00545_seed_130_seq.sto","w")

[new.write(x+"\n") for x in sto]