
import numpy as np
perfTable = [["dataset", \
        "Monotonicity", \
            "Criticality", \
            "Monotonic Paths", \
            "Critical Pairs", \
            "\\rowcolor{white}Total"]]
sec_colonne = [""]*6
# ,"{\\small\\autoref{sec_cpRegularity}:$\\quad$\\textbf{1)+2)}}", \
#         "{\\small\\autoref{sec_cpRegularity}:$\\quad$\\textbf{3)+4)}}", \
#         "{\\small\\autoref{sec_nonProgressive_diagram}:$\\quad$\\textbf{1)}}", \
#         "{\\small\\autoref{sec_nonProgressive_diagram}:$\\quad$\\textbf{2)+3)}}",\
#         ""]
dataset_list=['minMax','isabel','aneurism', 'random']
for i,dataset in enumerate(dataset_list):
    dataset=dataset.split(".")[0]
    if(dataset=="ctBones"):
        perfTable.append(["foot"])
        perfTable.append(["foot"])
    elif(dataset=="at"):
        perfTable.append(["AT"])
        perfTable.append(["AT"])
    elif(dataset=="ethanediol"):
        perfTable.append(["ethaneDiol"])
        perfTable.append(["ethaneDiol"])
    else:
        perfTable.append([dataset])
        perfTable.append([dataset])
    # perfTableCP.append([dataset])


    # nonprog 
    cpL=[]
    propL=[]
    pairsL=[]
    for itry in range(0,12):
        timeNP=0
        mono=0
        cp=0 
        prop=0
        pairs=0
        with open("./outputs_pd/nonProg_"+dataset+"_try"+str(itry), "r") as f:
            for line in f.readlines():
                line=line.split(" ")
                if(len(line)>1):
                    line2 = line
                    if(line2[0]=='initial'):
                        cp+=float(line2[-2])
                    if(line2[0]=='FIRSTPROPAGATION'):
                        prop+=float(line2[-1])
                    if(line2[0]=='TRIPLETS'):
                        pairs+=float(line2[-1])
                    if(line2[0]=='PAIRS'):
                        pairs+=float(line2[-1])
            f.close()
            cpL.append(cp)
            pairsL.append(pairs)
            propL.append(prop)
    cpL.sort()
    pairsL.sort()
    propL.sort()
    cpL=cpL[1:-1]
    pairsL=pairsL[1:-1]
    propL=propL[1:-1]
    tot=0
    perfTable[2*i+1].append("\multicolumn{1}{c}{$\emptyset$}")
    mean = sum(cpL)/len(cpL)
    tot+=mean
    perfTable[2*i+1].append(format(mean, '.2f'))
    mean = sum(propL)/len(propL)
    tot+=mean
    perfTable[2*i+1].append(format(mean, '.2f'))
    mean = sum(pairsL)/len(pairsL)
    tot+=mean
    perfTable[2*i+1].append(format(mean, '.2f'))
    perfTable[2*i+1].append(format(tot,'.2f'))
    # Prog
    monoL=[]
    cpL=[]
    propL=[]
    pairsL=[]
    for itry in range(0,12):
        timeNP=0
        mono=0
        cp=0
        pairs=0
        prop=0
        with open("./outputs_pd/prog_"+dataset+"_try"+str(itry),"r") as f:
            dl=-1
            for line in f.readlines():
                line=line.split(" ")
                if(len(line)>1):
                    line2 = line
                    if(line2[0]=='decimation'):
                        dl=int(line2[-1])
                    if(1):
                        if(line2[0]=='MONOTONY'):
                            mono+=float(line[-2])
                        if(line2[0]=='CRITICAL'):
                            cp+=float(line[-1])
                        if(line2[0]=='PROPAGATION'):
                            prop+=float(line[-1])
                        if(line2[0]=='TRIPLETS'):
                            pairs+=float(line[-1])
                        if(line2[0]=='PAIRS'):
                            pairs+=float(line[-1])
            f.close()
            monoL.append(mono)
            cpL.append(cp)
            pairsL.append(pairs)
            propL.append(prop)
    monoL.sort()
    cpL.sort()
    pairsL.sort()
    propL.sort()
    monoL=monoL[1:-1]
    cpL=cpL[1:-1]
    pairsL=pairsL[1:-1]
    propL=propL[1:-1]
    tot=0
    mean = sum(monoL)/len(monoL)
    tot+=mean
    perfTable[2*i+2].append(format(mean, '.2f'))
    mean = sum(cpL)/len(cpL)
    tot+=mean
    perfTable[2*i+2].append(format(mean, '.2f'))
    mean = sum(propL)/len(propL)
    tot+=mean
    perfTable[2*i+2].append(format(mean, '.2f'))
    mean = sum(pairsL)/len(pairsL)
    tot+=mean
    perfTable[2*i+2].append(format(mean, '.2f'))
    perfTable[2*i+2].append(format(tot, '.2f'))

perfTable=np.array(perfTable, dtype=object)
sec_colonne = np.array(sec_colonne, dtype=object)
perfTable=np.insert(perfTable, [0]*perfTable.shape[1], "x x")
perfTable=perfTable.reshape(2*len(dataset_list)+2,6)
perfTable[0,:]=perfTable[1,:]
perfTable[1,:]=sec_colonne
perfTable=np.array(perfTable).T

#print(perfTable[0,2::2])
#with open("./table4.tex","w") as f:
    #f.write("\\begin{tabular}{|ll"+"".join("|rr" for x in dataset_list)+"|}\n")
    ## f.write("%s\\\\\n" % " & ".join(str(col).title() for col in table[0]))
    #f.write("\hline\n")
    #f.write("&"+"".join("&\multicolumn{2}{c|}{"+x[0].upper()+x[1:]+"}" for x in perfTable[0,2::2])+"\\\\\n")
    #f.write(" \multicolumn{2}{|c|}{Step} "+"".join(" & NP & Prog" for x in dataset_list) + "\\\\\n")
    ## f.write("{Data set} & 1th & 8th & 1th & 8th& 1th & 8th\\\\\n")
    #f.write("\hline\n")
    #for row in perfTable[1:-1]:
        #f.write("%s\\\\\n" % " & ".join(col for col in row))
    #f.write("\hline\n")
    #for row in perfTable[-1:]:
        #f.write("%s\\\\\n" % " & ".join(col for col in row))
    #f.write("\hline\n")
    #f.write("\\end{tabular}\n")
    #f.close()

for i in range(2,2+2*len(dataset_list)):
    for j in range(1,6):
        if(perfTable[j,i] != '\\multicolumn{1}{c}{$\\emptyset$}'):
            perfTable[j,i] = format((float(perfTable[j,i])/float(perfTable[-1,i])*100),'.1f')
print(perfTable)

print(perfTable[0,2::2])
with open("./table3.tex","w") as f:
    f.write("\\begin{tabular}{|ll"+"".join("|rr" for x in dataset_list)+"|}\n")
    # f.write("%s\\\\\n" % " & ".join(str(col).title() for col in table[0]))
    f.write("\hline\n")
    f.write("&"+"".join("&\multicolumn{2}{c|}{"+x[0].upper()+x[1:]+"}" for x in perfTable[0,2::2])+"\\\\\n")
    f.write(" \multicolumn{2}{|c|}{Step} "+"".join(" & NP & Prog" for x in dataset_list) + "\\\\\n")
    f.write("\hline\n")
    for row in perfTable[-1:]:
        f.write("%s\\\\\n" % " & ".join(col for col in row))
    f.write("\hline\n")
    for row in perfTable[1:-1]:
        f.write("%s\\\\\n" % " & ".join(col for col in row))
    f.write("\hline\n")
    f.write("\\end{tabular}\n")
    f.close()
